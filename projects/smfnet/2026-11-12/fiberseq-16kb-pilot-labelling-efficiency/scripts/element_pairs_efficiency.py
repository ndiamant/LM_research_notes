"""CRE pairs, subsampled, with per-molecule labelling efficiency b_m as a covariate.

As element_pairs_partial.py, but: all CTCF-CTCF pairs, and random subsets of the other
groups (CTCF-other 2,000, other-other 4,000, random positions 3,000); and the partial
correlation given (a) the window-based labelling level, (b) b_m from
molecule_efficiency.py, (c) both.

Original description:

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
efficiency = {c: np.load(f"{S}/efficiency/{c}.efficiency.npz")["efficiency"] for c in lanes}
LIMIT = {"CTCF-other": 2000, "other-other": 4000, "random": 3000}

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
    b_m = efficiency[chrom][idx]
    ok = (nx >= 10) & (ny >= 10) & (nl >= 500) & ~np.isnan(b_m)
    n = int(ok.sum())
    if n < MIN_N:
        return None
    x, y, m, e = x[ok], y[ok], level[ok], b_m[ok]
    out = {"n": n, "raw": float(np.corrcoef(x, y)[0, 1]) if x.std() and y.std() else np.nan,
           "x_mean": float(x.mean()), "y_mean": float(y.mean())}
    for name, covariates in (("level", [m]), ("efficiency", [e]), ("both", [m, e])):
        design = np.column_stack([np.ones(n), *covariates])
        rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
        ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
        if rx.std() == 0 or ry.std() == 0:
            return None
        r = float(np.corrcoef(rx, ry)[0, 1])
        z = np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(n - 3 - len(covariates))
        out[f"r_{name}"] = r
        out[f"pvalue_{name}"] = float(2 * stats.norm.sf(abs(z)))
    return out


candidates = {"CTCF-CTCF": [], "CTCF-other": [], "other-other": []}
for chrom in sorted(lanes):
    el = sorted(elements[chrom])
    pos = np.array([e[0] for e in el])
    for i, (p, c1) in enumerate(el):
        for j in range(i + 1, np.searchsorted(pos, p + 15000, side="right")):
            q, c2 = el[j]
            if q - p >= 2000:
                kind = "CTCF-CTCF" if c1 and c2 else "CTCF-other" if c1 or c2 else "other-other"
                candidates[kind].append((chrom, p, q))
random_pairs = []
for kind in list(candidates):
    pool = candidates[kind]
    if kind in LIMIT and len(pool) > LIMIT[kind] * 3:  # oversample: some lack 50 molecules
        pool = [pool[i] for i in rng.choice(len(pool), LIMIT[kind] * 3, replace=False)]
    candidates[kind] = pool
    random_pairs += pool
random_pairs = [random_pairs[i] for i in rng.choice(len(random_pairs), LIMIT["random"] * 3, replace=False)]
rows = []
for kind, pool in list(candidates.items()) + [("random", random_pairs)]:
    kept = 0
    for chrom, p, q in pool:
        if kind in LIMIT and kept >= LIMIT[kind]:
            break
        shift = int(rng.choice([-1, 1]) * rng.integers(3000, 6000)) if kind == "random" else 0
        res = partial_r(chrom, p + shift, q + shift)
        if res is not None:
            rows.append({"chrom": chrom, "p": p + shift, "q": q + shift, "sep": q - p, "kind": kind, **res})
            kept += 1
    print(kind, kept, flush=True)

with open(f"{OUT}/element_pairs_efficiency.jsonl", "w") as handle:
    for r in rows:
        handle.write(json.dumps(r) + "\n")

real = [r for r in rows if r["kind"] != "random"]
for cov in ("level", "efficiency", "both"):
    pv = np.array([r[f"pvalue_{cov}"] for r in real])
    order = np.argsort(pv)
    qv = np.empty(len(pv))
    qv[order] = np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    for r, v in zip(real, qv):
        r[f"qvalue_{cov}"] = float(v)

print(f"\npairs 2-15 kb apart, >= {MIN_N} shared molecules; partial r given each covariate")
print(f"{'group':12s} {'pairs':>6s} | {'raw':>6s} | level: mean  p<.001  r>.3 | efficiency: mean  p<.001  r>.3 | both: mean  p<.001  r>.3")
for kind in ("CTCF-CTCF", "CTCF-other", "other-other", "random"):
    sel = [r for r in rows if r["kind"] == kind]
    raw = np.nanmean([r["raw"] for r in sel])
    cells = []
    for cov in ("level", "efficiency", "both"):
        rr = np.array([r[f"r_{cov}"] for r in sel])
        pp = np.array([r[f"pvalue_{cov}"] for r in sel])
        cells.append(f"{rr.mean():+.3f} {(pp < 0.001).mean():6.1%} {(rr > 0.3).mean():5.1%}")
    print(f"{kind:12s} {len(sel):6d} | {raw:+.3f} | " + " | ".join(cells))
sig = sorted([r for r in real if r["qvalue_both"] < 0.05], key=lambda r: -r["r_both"])
print(f"\n{len(sig)} of {len(real)} real pairs at FDR < 0.05 given both covariates; top 15:")
for r in sig[:15]:
    print(f"  {r['kind']:11s} {r['chrom']}:{r['p']:,}-{r['q']:,} ({r['sep'] / 1000:.1f} kb) n={r['n']} "
          f"r={r['r_both']:+.2f} (raw {r['raw']:+.2f}, level {r['r_level']:+.2f}, efficiency {r['r_efficiency']:+.2f}) "
          f"access {r['x_mean']:.2f}/{r['y_mean']:.2f} q={r['qvalue_both']:.1e}")
