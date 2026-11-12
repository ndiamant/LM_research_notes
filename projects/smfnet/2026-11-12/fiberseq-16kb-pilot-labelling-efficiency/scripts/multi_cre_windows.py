"""How many 16 kb windows (kept-window centres, training chromosomes) hold several CREs?

Elements: candidates.bed cCRE centres (thinned to >= 1,057 bp apart, so counts are lower
bounds), active = overlaps an IDR ATAC peak, CTCF = CTCF-bound class. Counted inside the
16,384 bp SMF window, centre element included.
"""
import os
from collections import defaultdict

import numpy as np

from smf_net.figwig_data import read_bed_regions

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
TRAIN = {f"chr{c}" for c in (2, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, "X")}
active, ctcf = defaultdict(list), defaultdict(list)
for line in open("/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/candidates.bed"):
    f = line.rstrip("\n").split("\t")
    c = (int(f[1]) + int(f[2])) // 2
    if f[4] == "active":
        active[f[0]].append(c)
    if "CTCF-bound" in f[7]:
        ctcf[f[0]].append(c)
active = {k: np.sort(v) for k, v in active.items()}
ctcf = {k: np.sort(v) for k, v in ctcf.items()}
classes = {}
for line in open(S + "/kept_windows.bed"):
    f = line.split("\t")
    classes[(f[0], int(f[1]))] = f[4]
rows = []
for chrom, start in read_bed_regions(S + "/kept_windows.bed", 2114):
    if chrom not in TRAIN:
        continue
    a = start + 1057 - 8192
    count = lambda arr: int(np.searchsorted(arr, a + 16384) - np.searchsorted(arr, a)) if arr is not None else 0
    rows.append((classes[(chrom, start)], count(active.get(chrom)), count(ctcf.get(chrom))))
cls = np.array([r[0] for r in rows]); n_active = np.array([r[1] for r in rows]); n_ctcf = np.array([r[2] for r in rows])
print(f"{len(rows):,} training windows (16 kb SMF window around each kept-window centre)")
for k in (1, 2, 3, 5, 8):
    print(f"  >= {k} active cCREs: {(n_active >= k).mean():6.1%}  (active centres {(n_active[cls == 'active'] >= k).mean():6.1%}, background centres {(n_active[cls == 'background'] >= k).mean():6.1%})")
for k in (1, 2, 3):
    print(f"  >= {k} CTCF-bound:  {(n_ctcf >= k).mean():6.1%}")
print(f"  mean active cCREs per window {n_active.mean():.2f} (median {np.median(n_active):.0f}); mean CTCF {n_ctcf.mean():.2f}")
print(f"  >= 3 active and >= 2 CTCF: {((n_active >= 3) & (n_ctcf >= 2)).mean():.1%} ({((n_active >= 3) & (n_ctcf >= 2)).sum():,} windows)")
