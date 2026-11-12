"""Split-half check of b_m's region independence: per val kept window, the mean b_m of a random
half of its molecules against the centre accessibility measured on the other half (and, for
contrast, accessibility on the same half)."""
import os
import numpy as np
from scipy import stats
from smf_net.figwig_data import open_lanes, read_bed_regions

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
rng = np.random.default_rng(1)
lanes = open_lanes(S + "/figwig")
eff = {c: np.load(f"{S}/efficiency/{c}.efficiency.npz")["efficiency"] for c in ("chr8", "chr20")}
same, other = [], []
for chrom, start in read_bed_regions(S + "/kept_windows.bed", 2114):
    if chrom not in eff:
        continue
    l = lanes[chrom]
    idx = l.overlapping(start, start + 2114, 500)
    idx = idx[~np.isnan(eff[chrom][idx])]
    if len(idx) < 40:
        continue
    calls = l.calls(idx, start, 2114).numpy()
    acc = (calls == 0).sum(1) / np.maximum((calls != -1).sum(1), 1)
    half = rng.permutation(len(idx)) < len(idx) // 2
    b = eff[chrom][idx]
    same.append((b[half].mean(), acc[half].mean()))
    other.append((b[half].mean(), acc[~half].mean()))
same, other = np.array(same), np.array(other)
print(f"{len(same):,} val windows")
print(f"  mean b_m (half A) vs accessibility of the same half A:  Spearman {stats.spearmanr(same[:, 0], same[:, 1])[0]:+.3f}")
print(f"  mean b_m (half A) vs accessibility of the other half B: Spearman {stats.spearmanr(other[:, 0], other[:, 1])[0]:+.3f}")
