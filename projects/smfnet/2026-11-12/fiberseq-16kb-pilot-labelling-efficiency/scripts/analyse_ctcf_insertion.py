"""Summarise ctcf_insertion.py: element accessibility and the two elements' co-accessibility by
orientation and distance, paired contrasts across backgrounds, and a figure."""

import json
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

run = sys.argv[1]
rows = [json.loads(line) for line in open(f"{run}/results.jsonl")]
DIST = sorted({r["distance"] for r in rows if r["distance"]})
ORIENT = ("convergent", "divergent", "tandem", "left_only", "right_only")
by = defaultdict(dict)  # (orientation, d) -> background -> row
for r in rows:
    if r["distance"]:
        by[(r["condition"].rsplit("_", 1)[0], r["distance"])][r["background"]] = r


def stat(o, d, key):
    v = np.array([r[key] for r in by[(o, d)].values()], float)
    v = v[~np.isnan(v)]
    return v.mean(), v.std(ddof=1) / np.sqrt(len(v)), len(v)


print(f"{len(rows)} sequences; mean +/- SE over backgrounds")
print(f"{'orientation':11s} {'d (kb)':>6s} {'left acc':>14s} {'right acc':>14s} {'raw r':>14s} {'partial r':>14s}")
for d in DIST:
    for o in ORIENT:
        cells = [stat(o, d, k) for k in ("left_access", "right_access", "raw_r", "partial_r")]
        print(f"{o:11s} {d / 1000:6.0f} " + " ".join(f"{m:+.3f}+/-{s:.3f}" for m, s, _ in cells))

print("\npaired contrasts of partial r (same backgrounds), mean +/- SE")
for d in DIST:
    for a, b in (("convergent", "divergent"), ("convergent", "tandem"), ("convergent", "left_only")):
        common = sorted(set(by[(a, d)]) & set(by[(b, d)]))
        diff = np.array([by[(a, d)][k]["partial_r"] - by[(b, d)][k]["partial_r"] for k in common])
        diff = diff[~np.isnan(diff)]
        print(f"  {d / 1000:4.0f} kb {a} - {b}: {diff.mean():+.3f} +/- {diff.std(ddof=1) / np.sqrt(len(diff)):.3f} (n={len(diff)})")

print("\npartner effect on accessibility (pair minus element alone), mean +/- SE")
for d in DIST:
    for o in ("convergent", "divergent", "tandem"):
        for side, single in (("left_access", "left_only"), ("right_access", "right_only")):
            if o == "tandem" and side == "right_access":
                continue  # tandem's right element is forward; right_only is reverse
            common = sorted(set(by[(o, d)]) & set(by[(single, d)]))
            diff = np.array([by[(o, d)][k][side] - by[(single, d)][k][side] for k in common])
            print(f"  {d / 1000:4.0f} kb {o:10s} {side.split('_')[0]:5s}: {diff.mean():+.4f} +/- {diff.std(ddof=1) / np.sqrt(len(diff)):.4f}")

# ---- figure -----------------------------------------------------------------------
SURFACE, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
COLORS = {"convergent": "#2a78d6", "divergent": "#eb6834", "tandem": "#1baf7a", "left_only": MUTED}
LABELS = {"convergent": "Convergent", "divergent": "Divergent", "tandem": "Tandem",
          "left_only": "Left element alone (vs the empty right site)"}
CORR = LinearSegmentedColormap.from_list("corr", ["#1c5cab", "#f0efec", "#e34948"])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED})
maps = np.load(f"{run}/maps.npz")
fig = plt.figure(figsize=(14, 4.6), facecolor=SURFACE)
grid = fig.add_gridspec(1, 4, width_ratios=[1.3, 1, 1, 1], wspace=0.35)
ax = fig.add_subplot(grid[0])
for o in ("convergent", "divergent", "tandem", "left_only"):
    m = np.array([stat(o, d, "partial_r")[0] for d in DIST])
    se = np.array([stat(o, d, "partial_r")[1] for d in DIST])
    ax.errorbar(np.array(DIST) / 1000, m, yerr=se, color=COLORS[o], linewidth=2, marker="o",
                markersize=5, capsize=3, label=LABELS[o])
ax.axhline(0, color=AXIS, linewidth=0.8)
ax.set_xlabel("Distance between the two elements (kb)")
ax.set_ylabel("Co-accessibility of the two elements\n(partial r given read labelling level)")
ax.set_title("Two inserted 3x CTCF elements", loc="left", fontsize=10)
ax.legend(frameon=False, fontsize=7.5)
ax.grid(axis="y", color=GRID, linewidth=0.6)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
d = 8000 if 8000 in DIST else DIST[-1]
L = maps[f"convergent_{d}"].shape[0] * 100
for k, name in enumerate((f"convergent_{d}", f"divergent_{d}", f"left_only_{d}")):
    a = fig.add_subplot(grid[k + 1])
    im = a.imshow(maps[name], cmap=CORR, vmin=-0.3, vmax=0.3, extent=(0, L / 1000, L / 1000, 0),
                  interpolation="nearest")
    for p in (L / 2 - d / 2, L / 2 + d / 2):
        a.axvline(p / 1000, color=INK2, linewidth=0.5, linestyle=":")
        a.axhline(p / 1000, color=INK2, linewidth=0.5, linestyle=":")
    a.set_title(f"{name.replace('_', ' ').replace(str(d), f'{d // 1000} kb')}, mean map", loc="left", fontsize=9)
    a.set_xlabel("Position (kb)")
fig.colorbar(im, ax=fig.axes[1:], shrink=0.8, label="Correlation given read labelling level")
fig.text(0.01, 0.01, f"{len({r['background'] for r in rows})} background windows (chr8/chr20), 256 generated reads each; "
         "error bars: SE over backgrounds; dotted lines: element positions.", color=MUTED, fontsize=7.5)
fig.savefig(f"{run}/ctcf_pairs.png", dpi=130, facecolor=SURFACE, bbox_inches="tight")
print("wrote", f"{run}/ctcf_pairs.png")
