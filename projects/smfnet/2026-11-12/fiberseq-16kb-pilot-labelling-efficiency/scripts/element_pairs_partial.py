"""Are there individual CRE pairs that co-vary on the same molecules beyond labelling level?

Elements: active cCREs of candidates.bed (overlapping an IDR ATAC peak), split by whether
both, one or neither are CTCF-bound. Every pair 2-15 kb apart, genome-wide. For each pair,
the molecules covering both elements +/-100 bp; each molecule's accessible fraction at
each element (+/-100 bp), and its labelling level: its accessible fraction over
[p - 4 kb, q + 4 kb] excluding +/-300 bp around both elements. The partial correlation of
the two elements' accessibility given labelling level (both regressed on it, residuals
correlated) is tested against 0 (Fisher z, n - 4 degrees of freedom) and against random
position pairs at the same distances (each real pair shifted 3-6 kb), an empirical null.
"""

import json
import os
from collections import defaultdict

import numpy as np
from scipy import stats

from smf_net.figwig_data import open_lanes

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
OUT = os.environ["SCRATCH"] + "/tm_inspect/coaccessibility"
HALF, EXCL, PAD, MIN_N = 100, 300, 4000, 50
rng = np.random.default_rng(0)
lanes = open_lanes(S + "/figwig")

elements = defaultdict(list)
for line in open("/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/candidates.bed"):
    f = line.rstrip("\n").split("\t")
    if f[4] == "active":
        elements[f[0]].append(((int(f[1]) + int(f[2])) // 2, "CTCF-bound" in f[7]))


def fraction(calls, lo, hi):
    c = calls[:, lo:hi]
    n = (c != -1).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, (c == 0).sum(1) / n, np.nan), n


def partial_r(chrom, p, q):
    l = lanes[chrom]
    a, b = p - PAD, q + PAD
    if a < 0 or b > l.chrom_length:
        return None
    # molecules covering both elements +/- HALF
    idx = l.overlapping(p - HALF, q + HALF, (q + HALF) - (p - HALF))
    if len(idx) < MIN_N:
        return None
    calls = l.calls(idx, a, b - a).numpy()
    x, nx = fraction(calls, p - a - HALF, p - a + HALF)
    y, ny = fraction(calls, q - a - HALF, q - a + HALF)
    rest = calls.copy()
    for c in (p - a, q - a):
        rest[:, max(c - EXCL, 0) : c + EXCL] = -1
    level, nl = fraction(rest, 0, b - a)
    ok = (nx >= 10) & (ny >= 10) & (nl >= 500)
    n = int(ok.sum())
    if n < MIN_N:
        return None
    x, y, m = x[ok], y[ok], level[ok]
    design = np.column_stack([np.ones(n), m])
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if rx.std() == 0 or ry.std() == 0:
        return None
    r = float(np.corrcoef(rx, ry)[0, 1])
    raw = float(np.corrcoef(x, y)[0, 1]) if x.std() and y.std() else np.nan
    z = np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(n - 4)
    return {"n": n, "r": r, "raw": raw, "pvalue": float(2 * stats.norm.sf(abs(z))),
            "x_mean": float(x.mean()), "y_mean": float(y.mean())}


rows = []
for chrom in sorted(lanes):
    el = sorted(elements[chrom])
    pos = np.array([e[0] for e in el])
    for i, (p, c1) in enumerate(el):
        for j in range(i + 1, np.searchsorted(pos, p + 15000, side="right")):
            q, c2 = el[j]
            if q - p < 2000:
                continue
            kind = "CTCF-CTCF" if c1 and c2 else "CTCF-other" if c1 or c2 else "other-other"
            res = partial_r(chrom, p, q)
            if res is not None:
                rows.append({"chrom": chrom, "p": p, "q": q, "sep": q - p, "kind": kind, **res})
            if rng.random() < 0.25:  # random-position null for a quarter of the pairs
                shift = int(rng.choice([-1, 1]) * rng.integers(3000, 6000))
                res = partial_r(chrom, p + shift, q + shift)
                if res is not None:
                    rows.append({"chrom": chrom, "p": p + shift, "q": q + shift, "sep": q - p,
                                 "kind": "random", **res})
    print(chrom, len(rows), flush=True)

with open(f"{OUT}/element_pairs_partial.jsonl", "w") as handle:
    for r in rows:
        handle.write(json.dumps(r) + "\n")

# Benjamini-Hochberg over all real pairs together
real = [r for r in rows if r["kind"] != "random"]
pv = np.array([r["pvalue"] for r in real])
order = np.argsort(pv)
q_values = np.empty(len(pv))
q_values[order] = np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
for r, qv in zip(real, q_values):
    r["qvalue"] = float(qv)

print(f"\npartial r given molecule labelling level; pairs 2-15 kb apart, >= {MIN_N} shared molecules")
print(f"{'group':12s} {'min n':>5s} {'pairs':>6s} {'mean r':>7s} {'p95':>5s} {'p99':>5s} "
      f"{'p<0.001':>8s} {'r>0.3':>6s}")
for min_n in (50, 100, 150):
    for kind in ("CTCF-CTCF", "CTCF-other", "other-other", "random"):
        sel = [r for r in rows if r["kind"] == kind and r["n"] >= min_n]
        if len(sel) < 10:
            continue
        rr = np.array([r["r"] for r in sel])
        pp = np.array([r["pvalue"] for r in sel])
        print(f"{kind:12s} {min_n:5d} {len(sel):6d} {rr.mean():+7.3f} {np.percentile(rr, 95):5.2f} "
              f"{np.percentile(rr, 99):5.2f} {(pp < 0.001).mean():8.2%} {(rr > 0.3).mean():6.2%}")
sig = sorted([r for r in real if r["qvalue"] < 0.05], key=lambda r: -r["r"])
print(f"\n{len(sig)} real pairs at FDR < 0.05 (of {len(real)}); top 15:")
for r in sig[:15]:
    print(f"  {r['kind']:11s} {r['chrom']}:{r['p']:,}-{r['q']:,} ({r['sep'] / 1000:.1f} kb) n={r['n']} "
          f"r={r['r']:+.2f} (raw {r['raw']:+.2f}) access {r['x_mean']:.2f}/{r['y_mean']:.2f} q={r['qvalue']:.1e}")
