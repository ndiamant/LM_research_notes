"""Kept H5 windows, re-centred at each size, that still have 32 reads covering >= 75%."""

import os

import numpy as np

from smf_net.figwig_data import open_lanes, read_bed_regions

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
lanes = open_lanes(S + "/figwig")
classes = {}
with open(S + "/kept_windows.bed") as handle:
    for line in handle:
        f = line.split("\t")
        classes[(f[0], int(f[1]))] = f[4]
centres = [(c, s + 1057, classes[(c, s)]) for c, s in read_bed_regions(S + "/kept_windows.bed", 2114)]
for width in (2048, 4096, 8192, 16384):
    kept = {"active": [0, 0], "background": [0, 0]}
    reads = []
    for chrom, centre, cls in centres:
        a = centre - width // 2
        b = a + width
        if a < 0 or b > lanes[chrom].chrom_length:
            continue
        idx = lanes[chrom].overlapping(a, b, 500)
        covered = np.minimum(lanes[chrom].end[idx], b) - np.maximum(lanes[chrom].start[idx], a)
        n = int((covered >= 0.75 * width).sum())
        kept[cls][1] += 1
        kept[cls][0] += n >= 32
        reads.append(n)
    print(
        f"{width}: active {kept['active'][0]:,}/{kept['active'][1]:,}, "
        f"background {kept['background'][0]:,}/{kept['background'][1]:,} with >= 32 reads; "
        f"median {np.median(reads):.0f} reads covering >= 75%"
    )
