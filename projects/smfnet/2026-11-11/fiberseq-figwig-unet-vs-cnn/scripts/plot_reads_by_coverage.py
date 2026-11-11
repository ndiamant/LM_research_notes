import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

OUT = os.environ["SCRATCH"] + "/tm_inspect/coverage_plot"
data = json.load(open(OUT + "/reads_by_coverage.json"))
x = np.array(data["fractions"]) * 100

# Ordinal blue ramp (validated with --ordinal), light to dark with window size.
COLORS = {"2048": "#86b6ef", "4096": "#3987e5", "8192": "#1c5cab", "16384": "#0d366b"}
SURFACE, INK, INK2, MUTED, GRID, AXIS = (
    "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
)
LABEL = {"2048": "2 kb", "4096": "4 kb", "8192": "8 kb", "16384": "16 kb"}

plt.rcParams.update({"font.size": 10, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})
fig, (ax_reads, ax_share) = plt.subplots(1, 2, figsize=(11.5, 4.6), facecolor=SURFACE)

for width, r in data["widths"].items():
    color = COLORS[width]
    label = f"{LABEL[width]} ({r['n_windows']:,} windows)"
    ax_reads.fill_between(x, r["q25"], r["q75"], color=color, alpha=0.18, linewidth=0)
    ax_reads.plot(x, r["median"], color=color, linewidth=2, label=label, solid_capstyle="round")
    ax_reads.annotate(LABEL[width], (x[-1], r["median"][-1]), xytext=(6, 0),
                      textcoords="offset points", va="center", color=INK2, fontsize=9)
    share = np.array(r["share_ge32"]) * 100
    ax_share.plot(x, share, color=color, linewidth=2, label=label, solid_capstyle="round")
    if width == "16384":  # the others overlap at the top; the legend names them
        ax_share.annotate(LABEL[width], (x[-1], share[-1]), xytext=(6, 0),
                          textcoords="offset points", va="center", color=INK2, fontsize=9)

ax_reads.axhline(32, color=MUTED, linewidth=1, linestyle=(0, (4, 3)))
ax_reads.annotate("32 reads drawn per item", (2, 32), xytext=(0, 4), textcoords="offset points",
                  color=MUTED, fontsize=8.5)
ax_reads.set_ylabel("Reads per window (median, band: 25th–75th percentile)")
ax_reads.set_title("Reads that cover at least the given fraction of the window",
                   loc="left", fontsize=10.5, color=INK)
ax_reads.set_ylim(0, None)

ax_share.set_ylabel("Windows with at least 32 such reads")
ax_share.yaxis.set_major_formatter(PercentFormatter(decimals=0))
ax_share.set_ylim(0, 102)
ax_share.set_title("Windows that can still fill an item with distinct reads",
                   loc="left", fontsize=10.5, color=INK)

for ax in (ax_reads, ax_share):
    ax.set_facecolor(SURFACE)
    ax.set_xlabel("Minimum fraction of the window a read must cover")
    ax.xaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.set_xlim(0, 100)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.margins(x=0)
ax_share.legend(title="Window", frameon=False, loc="lower left", fontsize=8.5, title_fontsize=9)

fig.text(0.01, 0.005,
         "GM12878 Fiber-seq, all chromosomes. Windows: the H5's kept windows re-centred at each size, kept if "
         "they fully contain an IDR ATAC peak (ENCFF614SMH). Reads also overlap by at least 500 bp.",
         color=MUTED, fontsize=8)
fig.tight_layout(rect=(0, 0.04, 0.97, 1))
fig.savefig(OUT + "/reads_by_coverage.png", dpi=200, facecolor=SURFACE)
print(OUT + "/reads_by_coverage.png")
