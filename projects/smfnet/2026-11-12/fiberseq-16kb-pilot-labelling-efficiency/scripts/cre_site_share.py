"""Share of observed calls inside IDR ATAC peaks (ENCFF614SMH), 1 kb centres vs 16 kb windows."""
import gzip
import os
from collections import defaultdict

import numpy as np

from smf_net.figwig_data import open_lanes, read_bed_regions, read_span

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
lanes = open_lanes(S + "/figwig")
peaks = defaultdict(list)
for line in gzip.open(S + "/inputs/ENCFF614SMH.bed.gz", "rt"):
    f = line.split("\t")
    peaks[f[0]].append((int(f[1]), int(f[2])))
classes = {}
for line in open(S + "/kept_windows.bed"):
    f = line.split("\t")
    classes[(f[0], int(f[1]))] = f[4]
rng = np.random.default_rng(0)
regions = [(c, s) for c, s in read_bed_regions(S + "/kept_windows.bed", 2114) if c in ("chr8", "chr20")]
regions = [regions[i] for i in rng.choice(len(regions), 300, replace=False)]
tot = defaultdict(lambda: np.zeros(2))
for chrom, start in regions:
    centre = start + 1057
    l = lanes[chrom]
    a = centre - 8192
    if a < 0 or a + 16384 > l.chrom_length:
        continue
    idx = l.overlapping(*read_span(a, 16384, 500, 0, 0.75))
    calls = l.calls(idx, a, 16384).numpy()
    inside = np.zeros(16384, bool)
    for p0, p1 in peaks[chrom]:
        if p1 > a and p0 < a + 16384:
            inside[max(p0 - a, 0) : min(p1 - a, 16384)] = True
    called = (calls != -1).sum(0)
    mid = slice(8192 - 512, 8192 + 512)
    cls = classes[(chrom, start)]
    for key, sl in ((f"centre 1 kb, {cls}", mid), ("centre 1 kb, all", mid), ("16 kb window, all", slice(None))):
        tot[key] += (called[sl][inside[sl]].sum(), called[sl].sum())
for k, (inside_calls, calls) in sorted(tot.items()):
    print(f"{k:28s} {inside_calls / calls:6.1%} of calls in IDR peaks")
