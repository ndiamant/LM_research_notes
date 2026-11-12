"""Co-accessibility of element pairs 2-15 kb apart on the same Fiber-seq molecules.

Elements are the cCRE centres of candidates.bed (ENCODE GM12878 classes). For each pair,
a 16,384 bp window centred on the pair's midpoint, the molecules covering >= 75% of it,
and each molecule's accessible fraction of called A/T sites within +/-100 bp of each
element. r is the correlation across molecules (>= 20 covering both), raw and with each
molecule's mean over its window removed (100 bp bins, as in coaccessibility_by_distance).

Groups: CTCF-CTCF (both CTCF-bound), other active pairs (both ATAC-active, neither
CTCF-bound), random positions (each CTCF pair shifted 3-6 kb, distance kept), and the
CTCF pairs with molecules shuffled (noise).
"""

import json
import os
from collections import defaultdict

import numpy as np

from smf_net.figwig_data import open_lanes, read_span

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
CANDIDATES = "/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/candidates.bed"
OUT = os.environ["SCRATCH"] + "/tm_inspect/coaccessibility"
WINDOW, HALF, BIN, MIN_N = 16384, 100, 100, 20
MIN_SEP, MAX_SEP = 2000, 15000
rng = np.random.default_rng(0)
lanes = open_lanes(S + "/figwig")

elements = defaultdict(list)
with open(CANDIDATES) as handle:
    for line in handle:
        f = line.rstrip("\n").split("\t")
        centre = (int(f[1]) + int(f[2])) // 2
        ctcf = "CTCF-bound" in f[7]
        elements[f[0]].append((centre, ctcf, f[4] == "active"))


def pairs(chrom, kind):
    rows = sorted(elements[chrom])
    pos = np.array([r[0] for r in rows])
    out = []
    for i, (p, ctcf, active) in enumerate(rows):
        for j in range(i + 1, np.searchsorted(pos, p + MAX_SEP, side="right")):
            q, ctcf2, active2 = rows[j]
            if q - p < MIN_SEP:
                continue
            if kind == "ctcf" and ctcf and ctcf2:
                out.append((p, q))
            elif kind == "active" and active and active2 and not (ctcf or ctcf2):
                out.append((p, q))
    return out


def pair_r(chrom, p, q, shuffle=False):
    l = lanes[chrom]
    a = (p + q) // 2 - WINDOW // 2
    if a < 0 or a + WINDOW > l.chrom_length:
        return None
    idx = l.overlapping(*read_span(a, WINDOW, 500, 0, 0.75))
    if len(idx) < MIN_N:
        return None
    calls = l.calls(idx, a, WINDOW).numpy()

    def access(lo, hi):
        c = calls[:, lo:hi]
        n = (c != -1).sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(n > 0, (c == 0).sum(1) / n, np.nan)

    x = access(p - a - HALF, p - a + HALF)
    y = access(q - a - HALF, q - a + HALF)
    n_bins = WINDOW // BIN
    bins = calls[:, : n_bins * BIN].reshape(len(calls), n_bins, BIN)
    nb = (bins != -1).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mol_mean = np.nanmean(np.where(nb > 0, (bins == 0).sum(-1) / nb, np.nan), axis=1)
    if shuffle:
        y = rng.permutation(y)
    ok = ~np.isnan(x) & ~np.isnan(y)
    if ok.sum() < MIN_N:
        return None
    out = []
    for u, v in ((x[ok], y[ok]), (x[ok] - mol_mean[ok], y[ok] - mol_mean[ok])):
        out.append(np.nan if u.std() == 0 or v.std() == 0 else float(np.corrcoef(u, v)[0, 1]))
    return out


results = {k: [] for k in ("ctcf", "active", "random", "shuffled")}
for chrom in sorted(lanes):
    ctcf_pairs = pairs(chrom, "ctcf")
    active_pairs = pairs(chrom, "active")
    if len(active_pairs) > 2 * len(ctcf_pairs):  # keep the groups comparable in size
        active_pairs = [active_pairs[i] for i in rng.choice(len(active_pairs), 2 * len(ctcf_pairs), replace=False)]
    for p, q in ctcf_pairs:
        for kind, args in (("ctcf", (p, q)), ("shuffled", (p, q))):
            r = pair_r(chrom, *args, shuffle=kind == "shuffled")
            if r is not None:
                results[kind].append((q - p, *r))
        shift = int(rng.choice([-1, 1]) * rng.integers(3000, 6000))
        r = pair_r(chrom, p + shift, q + shift)
        if r is not None:
            results["random"].append((q - p, *r))
    for p, q in active_pairs:
        r = pair_r(chrom, p, q)
        if r is not None:
            results["active"].append((q - p, *r))
    print(chrom, {k: len(v) for k, v in results.items()}, flush=True)

with open(f"{OUT}/element_pairs.json", "w") as handle:
    json.dump(results, handle)

print(f"\npairs {MIN_SEP // 1000}-{MAX_SEP // 1000} kb apart; r = correlation across molecules")
print(f"{'group':10s} {'n':>6s} | raw: mean   p95   p99  >0.3 | mean removed: mean   p95   p99  >0.3")
for kind, rows in results.items():
    a = np.array(rows, dtype=float)
    raw, res = a[:, 1], a[:, 2]
    raw, res = raw[~np.isnan(raw)], res[~np.isnan(res)]
    print(f"{kind:10s} {len(a):6d} | {raw.mean():9.3f} {np.percentile(raw, 95):5.2f} "
          f"{np.percentile(raw, 99):5.2f} {(raw > 0.3).mean():5.1%} | {res.mean():18.3f} "
          f"{np.percentile(res, 95):5.2f} {np.percentile(res, 99):5.2f} {(res > 0.3).mean():5.1%}")
