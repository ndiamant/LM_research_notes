"""(a) Same-sites loss relative to the 16 kb pilot by read efficiency b_m quintile (spec run and
2 kb model, central 1 kb; spec run, 16 kb window), from same_sites.npz. (b) Predicted vs observed
m6A at sites the spec model calls protected (q < 0.05), by b_m quintile: the numbers printed by
emission_calibration.py (160 val regions, fully masked reads), copied below."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = os.environ["SCRATCH"] + "/smfnet_figwig/spec_eval"
r = np.load(f"{D}/same_sites.npz")
edges = np.quantile(r["b"], [0, .2, .4, .6, .8, 1])
rel = {"spec_c": [], "short_c": [], "spec_w": []}
for k in range(5):
    s = (r["b"] >= edges[k]) & (r["b"] <= edges[k + 1])
    pc = r["pilot_c"][s].sum() / r["calls_c"][s].sum()
    pw = r["pilot_w"][s].sum() / r["calls_w"][s].sum()
    rel["spec_c"].append(100 * (r["spec_c"][s].sum() / r["calls_c"][s].sum() / pc - 1))
    rel["short_c"].append(100 * (r["short_c"][s].sum() / r["calls_c"][s].sum() / pc - 1))
    rel["spec_w"].append(100 * (r["spec_w"][s].sum() / r["calls_w"][s].sum() / pw - 1))
pred = [0.0310, 0.0325, 0.0328, 0.0333, 0.0345]
obs = [0.0179, 0.0205, 0.0211, 0.0220, 0.0249]

SURFACE, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED})
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8), facecolor=SURFACE)
x = np.arange(1, 6)
ax = axes[0]
ax.plot(x, rel["spec_w"], color=BLUE, lw=2, marker="o", ms=6, label="Spec run, 16 kb window")
ax.plot(x, rel["spec_c"], color=BLUE, lw=2, marker="o", ms=6, ls="--", label="Spec run, central 1 kb")
ax.plot(x, rel["short_c"], color=GREEN, lw=2, marker="o", ms=6, label="2 kb model, central 1 kb")
ax.axhline(0, color=AXIS, lw=0.8)
ax.set_xticks(x, [f"Q{k}\n{edges[k-1]:+.1f} to {edges[k]:+.1f}" for k in x])
ax.set_xlabel("Read labelling efficiency b_m quintile")
ax.set_ylabel("Loss vs the 16 kb pilot (%)")
ax.set_title("(a) The spec run is worse than the pilot, most for low-labelling reads", loc="left", fontsize=9.5)
ax.legend(frameon=False, fontsize=8)
ax = axes[1]
ax.plot(x, [100 * p for p in pred], color=ORANGE, lw=2, marker="o", ms=6, label="Predicted by the spec model")
ax.plot(x, [100 * o for o in obs], color=INK2, lw=2, marker="o", ms=6, label="Observed")
ax.set_xticks(x, [f"Q{k}" for k in x])
ax.set_ylim(0, 4)
ax.set_xlabel("Read labelling efficiency b_m quintile")
ax.set_ylabel("m6A at sites called protected (q < 0.05), %")
ax.set_title("(b) The fixed epsilon floor overpredicts m6A on protected DNA", loc="left", fontsize=9.5)
ax.legend(frameon=False, fontsize=8, loc="lower right")
for a in axes:
    a.set_facecolor(SURFACE); a.grid(axis="y", color=GRID, lw=0.6)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)
fig.tight_layout()
fig.savefig(f"{D}/spec_eval.png", dpi=120, facecolor=SURFACE)
print(f"{D}/spec_eval.png", rel)
