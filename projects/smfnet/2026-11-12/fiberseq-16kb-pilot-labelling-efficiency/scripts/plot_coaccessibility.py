import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

OUT = os.environ["SCRATCH"] + "/tm_inspect/coaccessibility"
name = sys.argv[1] if len(sys.argv) > 1 else "chr21"
r = json.load(open(f"{OUT}/coaccessibility_{name}.json"))
d = np.array(r["separation_bp"])
keep = d > 0
d = d[keep]
raw, res, sh = (np.array(r[k])[keep] for k in ("raw", "residual", "shuffled"))

BLUE, ORANGE = "#2a78d6", "#eb6834"  # categorical slots 1-2 (validated)
SURFACE, INK, INK2, MUTED, GRID, AXIS = (
    "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
)
REACH = 1541  # measured hard reach of the U-Net's convolutions
plt.rcParams.update({"font.size": 10, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})
kb = FuncFormatter(lambda v, _: f"{v / 1000:g} kb" if v >= 1000 else f"{v:g} bp")
fig, (ax, zoom) = plt.subplots(1, 2, figsize=(11.5, 4.4), facecolor=SURFACE)

for a in (ax, zoom):
    a.axhline(0, color=AXIS, linewidth=0.8)
    a.axvline(REACH, color=MUTED, linewidth=1, linestyle=(0, (4, 3)))
    a.plot(d, sh, color=MUTED, linewidth=1.2, linestyle=(0, (2, 2)), label="Molecules shuffled (noise)")
    a.plot(d, raw, color=BLUE, linewidth=2, label="Observed")
    a.plot(d, res, color=ORANGE, linewidth=2, label="Each molecule's mean removed")
    a.set_xscale("log")
    a.set_xlim(100, d[-1])
    a.xaxis.set_major_formatter(kb)
    a.set_xlabel("Distance between 100 bp bins")
    a.set_facecolor(SURFACE)
    a.grid(axis="y", color=GRID, linewidth=0.6)
    a.set_axisbelow(True)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)

ax.set_ylabel("Correlation of accessibility across molecules")
ax.set_title("Co-accessibility along Fiber-seq molecules", loc="left", fontsize=10.5)
ax.annotate("Observed", (12000, raw[d == 12000][0]), xytext=(0, 6), textcoords="offset points",
            color=INK2, fontsize=9, ha="center")
ax.annotate("Molecule mean removed", (6000, 0.0), xytext=(0, 8), textcoords="offset points",
            color=INK2, fontsize=9, ha="center")
ax.annotate("Convolutions reach ±1.5 kb", (REACH, -0.17), xytext=(5, 0), textcoords="offset points",
            color=MUTED, fontsize=8.5, va="center")
zoom.legend(frameon=False, loc="upper right", fontsize=8.5)

zoom.set_xlim(1000, d[-1])
zoom.set_ylim(-0.02, 0.07)
zoom.set_ylabel("Correlation (zoomed)")
zoom.set_title("With each molecule's mean removed, below 0.01 past about 2 kb", loc="left", fontsize=10.5)
zoom.annotate("Molecule mean removed", (1500, res[d == 1500][0]), xytext=(8, 6),
              textcoords="offset points", color=INK2, fontsize=9)
zoom.annotate("Noise", (8000, sh[d == 8000][0]), xytext=(0, -12), textcoords="offset points",
              color=MUTED, fontsize=8.5, ha="center")

fig.text(0.01, 0.005,
         f"GM12878 {', '.join(r['chroms'])}: {r['n_regions']:,} kept-window centres, a 16,384 bp window each, "
         f"{r['mean_molecules']:.0f} molecules covering >= 75% of it on average.\nBins: fraction of called "
         "A/T sites with m6A. Subtracting each molecule's mean biases distant bins slightly negative (about -0.01).",
         color=MUTED, fontsize=7.5)
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(f"{OUT}/coaccessibility_{name}.png", dpi=200, facecolor=SURFACE)
print(f"{OUT}/coaccessibility_{name}.png")
