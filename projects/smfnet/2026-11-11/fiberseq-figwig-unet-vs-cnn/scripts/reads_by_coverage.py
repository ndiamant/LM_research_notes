"""Reads available per window by the fraction of the window each read must cover.

Windows are the GM12878 H5's kept windows (active and GC-matched background),
re-centred at each window size, kept only if they fully contain an IDR ATAC peak
(ENCFF614SMH), the peaks that defined the H5's active set.
"""

import gzip
import json
import os
from collections import defaultdict

import numpy as np

from smf_net.figwig_data import open_lanes, read_bed_regions

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
OUT = os.environ["SCRATCH"] + "/tm_inspect/coverage_plot"
WIDTHS = (2048, 4096, 8192, 16384)
FRACTIONS = np.round(np.linspace(0, 1, 21), 2)
MIN_OVERLAP = 500

lanes = open_lanes(S + "/figwig")

peaks = defaultdict(list)
with gzip.open(S + "/inputs/ENCFF614SMH.bed.gz", "rt") as handle:
    for line in handle:
        f = line.split("\t")
        peaks[f[0]].append((int(f[1]), int(f[2])))
peak_start, peak_end = {}, {}
for chrom, rows in peaks.items():
    rows.sort()
    peak_start[chrom] = np.array([r[0] for r in rows])
    peak_end[chrom] = np.array([r[1] for r in rows])

classes = {}
with open(S + "/kept_windows.bed") as handle:
    for line in handle:
        f = line.split("\t")
        classes[(f[0], int(f[1]))] = f[4]
windows = read_bed_regions(S + "/kept_windows.bed", 2114)
centres = [(chrom, start + 1057, classes[(chrom, start)]) for chrom, start in windows]

results = {}
for width in WIDTHS:
    counts, kept_class = [], []
    for chrom, centre, cls in centres:
        a, b = centre - width // 2, centre + width - width // 2
        if chrom not in lanes or a < 0 or b > lanes[chrom].chrom_length:
            continue
        ps, pe = peak_start.get(chrom), peak_end.get(chrom)
        if ps is None:
            continue
        lo, hi = np.searchsorted(ps, a), np.searchsorted(ps, b)
        if not (pe[lo:hi] <= b).any():
            continue
        chrom_lanes = lanes[chrom]
        idx = chrom_lanes.overlapping(a, b, MIN_OVERLAP)
        covered = np.minimum(chrom_lanes.end[idx], b) - np.maximum(chrom_lanes.start[idx], a)
        need = np.maximum(MIN_OVERLAP, np.ceil(FRACTIONS * width))
        counts.append((covered[None, :] >= need[:, None]).sum(axis=1))
        kept_class.append(cls)
    counts = np.array(counts)  # [windows, fractions]
    kept_class = np.array(kept_class)
    results[width] = {
        "n_windows": len(counts),
        "n_active": int((kept_class == "active").sum()),
        "n_background": int((kept_class == "background").sum()),
        "median": np.median(counts, axis=0).tolist(),
        "q25": np.percentile(counts, 25, axis=0).tolist(),
        "q75": np.percentile(counts, 75, axis=0).tolist(),
        "share_ge32": (counts >= 32).mean(axis=0).tolist(),
    }
    r = results[width]
    print(
        f"{width}: {r['n_windows']} windows contain a peak "
        f"({r['n_active']} active, {r['n_background']} background centres); median reads "
        + ", ".join(f"{f:.1f}:{m:.0f}" for f, m in zip(FRACTIONS, r["median"]) if f in (0, 0.5, 0.8, 1))
    )

with open(OUT + "/reads_by_coverage.json", "w") as handle:
    json.dump({"fractions": FRACTIONS.tolist(), "n_centres": len(centres), "widths": results}, handle)
