"""Does b_m measure labelling and not region activity? (chr8, chr20)

1. Labelling: b_m vs each molecule's raw rate (accessible / called over its whole span) and its
   background-site rate; protected-stretch lengths (runs of protected calls between m6A calls,
   from the first to the last base of the run) by b_m quintile.
2. Region independence: per val kept window, the mean b_m and the mean raw rate of molecules
   covering its 2,114 bp centre, against the centre's observed accessible fraction, and by
   window class (active / background).
"""

import os

import numpy as np
from scipy import stats

from smf_net.figwig_data import open_lanes, read_bed_regions

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
rng = np.random.default_rng(0)
lanes = open_lanes(S + "/figwig")
eff = {c: np.load(f"{S}/efficiency/{c}.efficiency.npz") for c in ("chr8", "chr20")}
classes = {}
for line in open(S + "/kept_windows.bed"):
    f = line.split("\t")
    classes[(f[0], int(f[1]))] = f[4]

# 1. per-molecule: raw rate and protected stretches on a subsample
b_all, raw_all, bg_all, stretch_rows = [], [], [], []
for chrom in ("chr8", "chr20"):
    l, e = lanes[chrom], eff[chrom]
    for i in rng.choice(len(l), 4000, replace=False):
        if np.isnan(e["efficiency"][i]):
            continue
        s0, s1 = int(l.start[i]), int(l.end[i])
        calls = l.calls(np.array([i]), s0, s1 - s0)[0].numpy()
        called = calls != -1
        raw = (calls == 0).sum() / max(called.sum(), 1)
        b_all.append(e["efficiency"][i])
        raw_all.append(raw)
        bg_all.append(e["background_rate"][i])
        # protected stretches: between consecutive m6A calls, if all calls in between are protected
        m6a = np.flatnonzero(calls == 0)
        lengths = np.diff(m6a) - 1
        stretch_rows.append((e["efficiency"][i], lengths[lengths >= 20]))
b_all, raw_all, bg_all = map(np.array, (b_all, raw_all, bg_all))
logit = lambda p: np.log(np.clip(p, 1e-4, 1 - 1e-4) / (1 - np.clip(p, 1e-4, 1 - 1e-4)))
print(f"{len(b_all):,} molecules (chr8, chr20 subsample)")
print(f"  Spearman(b_m, raw molecule rate)        {stats.spearmanr(b_all, raw_all)[0]:+.3f}")
print(f"  Spearman(b_m, background-site rate)     {stats.spearmanr(b_all, bg_all)[0]:+.3f}")
print(f"  Pearson(b_m, logit raw rate)            {np.corrcoef(b_all, logit(raw_all))[0, 1]:+.3f}")
edges = np.quantile(b_all, [0, 0.2, 0.4, 0.6, 0.8, 1])
print("  protected stretches between m6A calls (>= 20 bp), by b_m quintile:")
for k in range(5):
    sel = [ls for b, ls in stretch_rows if edges[k] <= b <= edges[k + 1]]
    ls = np.concatenate(sel)
    hist, bins = np.histogram(ls[(ls >= 60) & (ls <= 260)], bins=np.arange(60, 262, 5))
    mode = bins[np.argmax(hist)] + 2.5
    raw_q = raw_all[(b_all >= edges[k]) & (b_all <= edges[k + 1])]
    print(f"    b {edges[k]:+.2f} to {edges[k + 1]:+.2f}: raw rate median {np.median(raw_q):.3f}, "
          f"mode of 60-260 bp stretches {mode:.0f} bp, median stretch {np.median(ls):.0f} bp, "
          f"share > 250 bp {(ls > 250).mean():.1%}")

# 2. per window: mean b_m vs centre accessibility
rows = []
for chrom, start in read_bed_regions(S + "/kept_windows.bed", 2114):
    if chrom not in eff:
        continue
    l, e = lanes[chrom], eff[chrom]
    idx = l.overlapping(start, start + 2114, 500)
    b = e["efficiency"][idx]
    ok = ~np.isnan(b)
    if ok.sum() < 20:
        continue
    calls = l.calls(idx, start, 2114).numpy()
    access = (calls == 0).sum() / max((calls != -1).sum(), 1)
    rate = e["background_rate"][idx]
    rows.append((classes[(chrom, start)], access, b[ok].mean(), np.nanmean(rate)))
cls = np.array([r[0] for r in rows])
acc, mb, mr = (np.array([r[k] for r in rows]) for k in (1, 2, 3))
print(f"\n{len(rows):,} val windows")
print(f"  Spearman(window accessibility, mean b_m)                {stats.spearmanr(acc, mb)[0]:+.3f}")
print(f"  Spearman(window accessibility, mean background-site rate) {stats.spearmanr(acc, mr)[0]:+.3f}")
for c in ("active", "background"):
    print(f"  {c:10s} windows: mean b_m {mb[cls == c].mean():+.3f} (sd {mb[cls == c].std():.3f}), "
          f"centre accessibility {acc[cls == c].mean():.3f}")
