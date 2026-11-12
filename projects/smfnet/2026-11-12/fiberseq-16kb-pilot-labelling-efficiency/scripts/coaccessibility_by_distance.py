"""Observed co-accessibility of Fiber-seq molecules by distance, out to 16 kb.

Regions and reads as the 16 kb run draws them: the kept-window centres
(kept_windows.bed), a 16,384 bp SMF window each, and the molecules covering at least 75%
of it. Each molecule's calls are summarised in 100 bp bins as the fraction of its called
A/T sites that are accessible (m6A); bins correlate across a region's molecules
(pairwise complete, >= 20 molecules), and correlations are averaged by bin separation.

- raw: the correlation as is.
- molecule mean removed: each molecule's mean over its covered bins subtracted first,
  which removes molecule-wide differences (labelling efficiency, a fibre that is open
  throughout) that correlate bins at any distance.
- shuffled: molecules permuted independently in each bin; the level noise alone gives.
"""

import json
import os
import sys

import numpy as np

from smf_net.figwig_data import open_lanes, read_bed_regions, read_span

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
OUT = os.environ["SCRATCH"] + "/tm_inspect/coaccessibility"
CHROMS = sys.argv[1].split(",") if len(sys.argv) > 1 else ["chr21"]
WINDOW, BIN, MIN_PAIRS = 16384, 100, 20
N_BINS = WINDOW // BIN
rng = np.random.default_rng(0)

lanes = open_lanes(S + "/figwig")
classes = {}
with open(S + "/kept_windows.bed") as handle:
    for line in handle:
        f = line.split("\t")
        classes[(f[0], int(f[1]))] = f[4]


def binned(calls):
    """[molecules, bins] accessible fraction, NaN where a bin has no call."""
    c = calls[:, : N_BINS * BIN].reshape(len(calls), N_BINS, BIN)
    called = (c != -1).sum(-1)
    accessible = (c == 0).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(called > 0, accessible / called, np.nan)


def correlations(x):
    """[bins, bins] pairwise-complete Pearson correlation across molecules."""
    m = ~np.isnan(x)
    v = np.where(m, x, 0.0)
    mf = m.astype(float)
    n = mf.T @ mf
    s1 = v.T @ mf  # sum of x_i over molecules with both i and j
    s2 = (v * v).T @ mf
    sxy = v.T @ v
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = sxy - s1 * s1.T / n
        var_i = s2 - s1 * s1 / n
        r = cov / np.sqrt(var_i * var_i.T)
    r[n < MIN_PAIRS] = np.nan
    return r


sep = np.abs(np.subtract.outer(np.arange(N_BINS), np.arange(N_BINS)))
sums = {k: np.zeros(N_BINS) for k in ("raw", "residual", "shuffled")}
counts = {k: np.zeros(N_BINS) for k in sums}
by_class = {c: {k: np.zeros(N_BINS) for k in ("raw", "residual")} for c in ("active", "background")}
class_counts = {c: {k: np.zeros(N_BINS) for k in ("raw", "residual")} for c in by_class}
n_regions = n_molecules = 0

for chrom in CHROMS:
    l = lanes[chrom]
    for c, start in read_bed_regions(S + "/kept_windows.bed", 2114):
        if c != chrom:
            continue
        cls = classes[(c, start)]
        a = start + 1057 - WINDOW // 2
        if a < 0 or a + WINDOW > l.chrom_length:
            continue
        idx = l.overlapping(*read_span(a, WINDOW, 500, 0, 0.75))
        if len(idx) < 32:
            continue
        x = binned(l.calls(idx, a, WINDOW).numpy())
        residual = x - np.nanmean(x, axis=1, keepdims=True)
        shuffled = np.stack([rng.permutation(col) for col in x.T], axis=1)
        for name, data in (("raw", x), ("residual", residual), ("shuffled", shuffled)):
            r = correlations(data)
            ok = ~np.isnan(r) & (sep > 0)
            np.add.at(sums[name], sep[ok], r[ok])
            np.add.at(counts[name], sep[ok], 1)
            if name in ("raw", "residual"):
                np.add.at(by_class[cls][name], sep[ok], r[ok])
                np.add.at(class_counts[cls][name], sep[ok], 1)
        n_regions += 1
        n_molecules += len(idx)

result = {
    "chroms": CHROMS, "window": WINDOW, "bin": BIN, "n_regions": n_regions,
    "mean_molecules": n_molecules / max(n_regions, 1),
    "separation_bp": (np.arange(N_BINS) * BIN).tolist(),
}
for name in sums:
    with np.errstate(invalid="ignore"):
        result[name] = (sums[name] / counts[name]).tolist()
for cls in by_class:
    for name in ("raw", "residual"):
        with np.errstate(invalid="ignore"):
            result[f"{cls}_{name}"] = (by_class[cls][name] / class_counts[cls][name]).tolist()
with open(f"{OUT}/coaccessibility_{'_'.join(CHROMS)}.json", "w") as handle:
    json.dump(result, handle)

print(f"{n_regions} regions on {','.join(CHROMS)}, {result['mean_molecules']:.0f} molecules each")
for d in (100, 200, 500, 1000, 1500, 2000, 4000, 8000, 12000, 16000):
    k = min(d // BIN, N_BINS - 1)
    print(f"{k * BIN:6d} bp: raw {result['raw'][k]:.4f}  molecule mean removed "
          f"{result['residual'][k]:.4f}  shuffled {result['shuffled'][k]:+.4f}  "
          f"(active {result['active_residual'][k]:.4f}, background "
          f"{result['background_residual'][k]:.4f})")
