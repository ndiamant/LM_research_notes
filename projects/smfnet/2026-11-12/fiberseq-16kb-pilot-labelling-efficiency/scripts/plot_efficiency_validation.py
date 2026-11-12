"""Figure: does the per-molecule efficiency b_m measure labelling and not region activity?

(a) b_m over all molecules genome-wide; (b) b_m vs the molecule's raw m6A rate (chr8/chr20
subsample, 2D histogram); (c) protected-stretch lengths (runs between consecutive m6A calls)
by b_m quintile; (d) per val window, mean b_m against centre accessibility measured on the same
molecules and on the other half of the window's molecules, with the raw background-site rate
for contrast (Spearman in the panel titles).
"""

import os
from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from smf_net.figwig_data import open_lanes, read_bed_regions

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
OUT = os.environ["SCRATCH"] + "/fiberseq/efficiency"
rng = np.random.default_rng(0)
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]  # validated ordinal ramp
plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})

all_b = np.concatenate([np.load(f)["efficiency"] for f in sorted(glob(f"{S}/efficiency/*.efficiency.npz"))])
all_b = all_b[~np.isnan(all_b)]

lanes = open_lanes(S + "/figwig")
eff = {c: np.load(f"{S}/efficiency/{c}.efficiency.npz") for c in ("chr8", "chr20")}
b_s, raw_s, stretches = [], [], []
for chrom in ("chr8", "chr20"):
    l, e = lanes[chrom], eff[chrom]
    for i in rng.choice(len(l), 4000, replace=False):
        if np.isnan(e["efficiency"][i]):
            continue
        calls = l.calls(np.array([i]), int(l.start[i]), int(l.end[i] - l.start[i]))[0].numpy()
        b_s.append(e["efficiency"][i])
        raw_s.append((calls == 0).sum() / max((calls != -1).sum(), 1))
        lengths = np.diff(np.flatnonzero(calls == 0)) - 1
        stretches.append(lengths[lengths >= 20])
b_s, raw_s = np.array(b_s), np.array(raw_s)

windows = []
for chrom, start in read_bed_regions(S + "/kept_windows.bed", 2114):
    if chrom not in eff:
        continue
    l, e = lanes[chrom], eff[chrom]
    idx = l.overlapping(start, start + 2114, 500)
    idx = idx[~np.isnan(e["efficiency"][idx])]
    if len(idx) < 40:
        continue
    calls = l.calls(idx, start, 2114).numpy()
    acc = (calls == 0).sum(1) / np.maximum((calls != -1).sum(1), 1)
    half = rng.permutation(len(idx)) < len(idx) // 2
    b, rate = e["efficiency"][idx], e["background_rate"][idx]
    windows.append((b[half].mean(), acc[half].mean(), acc[~half].mean(), rate[half].mean()))
w = np.array(windows)

fig, axes = plt.subplots(1, 4, figsize=(17, 4.2), facecolor=SURFACE)
ax = axes[0]
ax.hist(all_b, bins=np.linspace(-4, 3, 141), color=BLUE, linewidth=0)
ax.set_xlabel("Per-molecule efficiency b (logit offset)")
ax.set_ylabel("Molecules")
ax.set_title(f"(a) {len(all_b) / 1e6:.1f}M molecules: median {np.median(all_b):+.2f},\n"
             f"IQR {np.percentile(all_b, 25):+.2f} to {np.percentile(all_b, 75):+.2f}", loc="left", fontsize=9.5)

ax = axes[1]
h = ax.hist2d(b_s, raw_s, bins=[np.linspace(-3, 2.5, 60), np.linspace(0, 0.35, 60)], cmap="Blues", cmin=1)
ax.set_xlabel("Per-molecule efficiency b")
ax.set_ylabel("Raw m6A rate of the molecule")
ax.set_title(f"(b) tracks labelling: Spearman {stats.spearmanr(b_s, raw_s)[0]:.2f}\n"
             f"({len(b_s):,} molecules, chr8 + chr20)", loc="left", fontsize=9.5)
fig.colorbar(h[3], ax=ax, shrink=0.8, label="Molecules")

ax = axes[2]
edges = np.quantile(b_s, [0, 0.2, 0.4, 0.6, 0.8, 1])
bins = np.arange(20, 402, 5)
for k in range(5):
    sel = np.concatenate([ls for b, ls in zip(b_s, stretches) if edges[k] <= b <= edges[k + 1]])
    dens, _ = np.histogram(sel, bins=bins, density=True)
    peak_bins = np.arange(60, 262, 5)
    mode = peak_bins[np.argmax(np.histogram(sel[(sel >= 60) & (sel <= 260)], bins=peak_bins)[0])] + 2.5
    ax.plot(bins[:-1] + 2.5, dens, color=RAMP[k], linewidth=2,
            label=f"Q{k + 1} (b {edges[k]:+.1f} to {edges[k + 1]:+.1f}): peak {mode:.0f} bp")
ax.set_yscale("log")
ax.set_xlabel("Protected stretch between m6A calls (bp)")
ax.set_ylabel("Density")
ax.set_title("(c) nucleosome footprints shrink as labelling rises", loc="left", fontsize=9.5)
ax.legend(frameon=False, fontsize=7.5, title="Efficiency quintile", title_fontsize=8)

ax = axes[3]
r_same = stats.spearmanr(w[:, 0], w[:, 1])[0]
r_other = stats.spearmanr(w[:, 0], w[:, 2])[0]
r_rate = stats.spearmanr(w[:, 3], w[:, 2])[0]
ax.scatter(w[:, 2], w[:, 0], s=3, color=BLUE, alpha=0.25, linewidths=0)
ax.set_xlabel("Window centre accessibility (other half of its molecules)")
ax.set_ylabel("Mean b of half its molecules")
ax.set_title(f"(d) no region activity: Spearman {r_other:+.2f}\n"
             f"(same molecules {r_same:+.2f}; raw background rate {r_rate:+.2f})", loc="left", fontsize=9.5)
for a in axes:
    a.set_facecolor(SURFACE)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)
fig.text(0.01, 0.0, "b: logit P(m6A at site s on molecule m) = a_s + b_m, a_s from the pileup, b_m fitted on each "
         "molecule's sites outside +/-500 bp of IDR peaks and active cCREs. (d): 12,303 val windows.",
         color=MUTED, fontsize=7.5)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(f"{OUT}/efficiency_validation.png", dpi=140, facecolor=SURFACE)
print(f"{OUT}/efficiency_validation.png", r_same, r_other, r_rate)
