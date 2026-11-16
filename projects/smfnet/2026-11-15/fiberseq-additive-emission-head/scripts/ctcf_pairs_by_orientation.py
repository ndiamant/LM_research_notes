"""Observed CTCF-CTCF pair co-accessibility split by CTCF motif orientation.

Pairs and their correlations come from element_pairs_efficiency.py (active cCREs that are
CTCF-bound in candidates.bed, 2-15 kb apart, >= 50 shared molecules; partial r of the two
elements' accessibility given molecule labelling: window level, b_m, both). Each element's
motif: the best MA0139.2 (smf_hmm references, as the HMM work) log-odds hit, uniform
background, either strand, within +/-150 bp of the element centre. Its strand is called if
the hit beats the 95th percentile of best hits in 5,000 random 300 bp windows of the same
chromosomes (so ~5% of calls could be chance hits) and the other strand's best hit is >= 2
bits weaker; otherwise the element is ambiguous. Orientation of a pair, JASPAR's forward
strand as in the insertion experiment: convergent (left +, right -, pointing at each
other), divergent (left -, right +), tandem (same strand).
"""
import json
import os
from collections import defaultdict

import numpy as np
from pyfaidx import Fasta
from scipy import stats

OUT = os.environ["SCRATCH"] + "/tm_inspect/coaccessibility"
FASTA = "/scratch/groups/btrippe/ndiamant/reference/hg38.fa"
MEME = "/home/users/diamant/repos/smf_hmm/evaluation/tf_chip/references/MA0139.2.meme"
HALF = 150
rng = np.random.default_rng(0)

rows, on = [], False
for line in open(MEME):
    if line.startswith("letter-probability"):
        on = True
        continue
    v = line.split()
    if on and len(v) == 4:
        rows.append([float(x) for x in v])
    elif on and rows:
        break
pwm = np.log2((np.array(rows) + 1e-3) / 0.25)  # w x 4
W = len(pwm)
rc_pwm = pwm[::-1, ::-1]
CODE = np.full(256, 4, np.int8)
for i, b in enumerate(b"ACGT"):
    CODE[b] = i
    CODE[b + 32] = i
fa = Fasta(FASTA, as_raw=True)


def best(seq):
    """Best score per strand over a sequence (N scores -inf)."""
    s = CODE[np.frombuffer(seq.encode(), np.uint8)]
    n = len(s) - W + 1
    idx = np.arange(n)[:, None] + np.arange(W)[None, :]
    win = s[idx]
    ok = (win < 4).all(1)
    w = np.where(win < 4, win, 0)
    fwd = np.where(ok, pwm[np.arange(W), w].sum(1), -np.inf)
    rev = np.where(ok, rc_pwm[np.arange(W), w].sum(1), -np.inf)
    return fwd.max(), rev.max()


pairs = [json.loads(l) for l in open(f"{OUT}/element_pairs_efficiency.jsonl")]
pairs = [p for p in pairs if p["kind"] == "CTCF-CTCF"]
chroms = sorted({p["chrom"] for p in pairs})
null = []
for _ in range(5000):
    c = chroms[rng.integers(len(chroms))]
    L = len(fa[c])
    x = int(rng.integers(1_000_000, L - 1_000_000))
    f, r = best(fa[c][x : x + 2 * HALF])
    if np.isfinite(max(f, r)):
        null.append(max(f, r))
thr = float(np.percentile(null, 95))

strand = {}
for p in pairs:
    for pos in (p["p"], p["q"]):
        key = (p["chrom"], pos)
        if key in strand:
            continue
        f, r = best(fa[p["chrom"]][pos - HALF : pos + HALF + W])
        top = max(f, r)
        if top < thr or abs(f - r) < 2:
            strand[key] = None
        else:
            strand[key] = "+" if f > r else "-"

groups = defaultdict(list)
for p in pairs:
    a, b = strand[(p["chrom"], p["p"])], strand[(p["chrom"], p["q"])]
    if a is None or b is None:
        kind = "ambiguous"
    elif a == "+" and b == "-":
        kind = "convergent"
    elif a == "-" and b == "+":
        kind = "divergent"
    else:
        kind = "tandem"
    p["orientation"] = kind
    groups[kind].append(p)
calls = [v for v in strand.values()]
print(f"{len(pairs)} CTCF-CTCF pairs, {len(strand)} elements; motif threshold {thr:.1f} bits "
      f"(95th pct of random 300 bp windows); strand called for {sum(c is not None for c in calls)} "
      f"(+ {calls.count('+')}, - {calls.count('-')})")
print(f"{'orientation':11s} {'pairs':>5s} {'raw r':>14s} {'partial r (level)':>18s} {'partial r (both)':>17s} {'p<0.001':>8s} {'median kb':>9s}")
for kind in ("convergent", "divergent", "tandem", "ambiguous"):
    g = groups[kind]
    if not g:
        continue
    cells = []
    for key in ("raw", "r_level", "r_both"):
        v = np.array([x[key] for x in g], float)
        v = v[~np.isnan(v)]
        cells.append(f"{v.mean():+.3f}+/-{v.std(ddof=1) / np.sqrt(len(v)):.3f}")
    sig = np.mean([x["pvalue_both"] < 0.001 for x in g])
    print(f"{kind:11s} {len(g):5d} {cells[0]:>14s} {cells[1]:>18s} {cells[2]:>17s} {sig:8.1%} {np.median([x['sep'] for x in g]) / 1000:9.1f}")
conv = np.array([x["r_both"] for x in groups["convergent"]])
for other in ("divergent", "tandem"):
    o = np.array([x["r_both"] for x in groups[other]])
    if len(o) > 2 and len(conv) > 2:
        t = stats.ttest_ind(conv, o, equal_var=False)
        u = stats.mannwhitneyu(conv, o)
        print(f"convergent - {other}: {conv.mean() - o.mean():+.3f} (Welch p {t.pvalue:.2g}, Mann-Whitney p {u.pvalue:.2g})")
print("\nby distance, partial r given both covariates (n):")
for lo, hi in ((2000, 5000), (5000, 10000), (10000, 15001)):
    cells = []
    for kind in ("convergent", "divergent", "tandem"):
        v = [x["r_both"] for x in groups[kind] if lo <= x["sep"] < hi]
        cells.append(f"{kind} {np.mean(v):+.3f} ({len(v)})" if v else f"{kind} - (0)")
    print(f"  {lo // 1000}-{hi // 1000} kb: " + ", ".join(cells))
with open(f"{OUT}/ctcf_pairs_by_orientation.jsonl", "w") as h:
    for p in pairs:
        h.write(json.dumps(p) + "\n")
