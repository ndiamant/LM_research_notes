"""CRE pair metrics over training, from score_cre_pairs.py runs at periodic checkpoints.

Usage: plot_checkpoint_scan.py OUT.png PAIRS.tsv LABEL=DIR_PREFIX ... where each series'
runs are DIR_PREFIX<epoch index> (e.g. .../scan/pilot16_epoch -> pilot16_epoch4 ...).
Panels: Spearman of the score with observed partial r (all pairs, and >= 75 shared
molecules; dotted: their split-half ceilings), the same given raw r, AUROC, and AUPRC
(dotted: positive rate). Error bars: 95% intervals over 300 pair bootstraps.
"""
import sys
from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from smf_net.evaluation.cre_pairs import _partial_spearman, ranking_metrics, split_half_ceiling

out_png, pairs_path, series = sys.argv[1], sys.argv[2], sys.argv[3:]
table = pd.read_csv(pairs_path, sep="\t")
table["positive"] = (table["partial_r"] > 0) & (table["pvalue"] < 0.01)
rng = np.random.default_rng(0)
boot = rng.integers(0, len(table), (300, len(table)))
METRICS = {
    "spearman": lambda d, r, raw, pos: stats.spearmanr(d, r)[0],
    "spearman_given_raw": lambda d, r, raw, pos: _partial_spearman(d, r, raw),
    "auroc": lambda d, r, raw, pos: ranking_metrics(d, pos)["auroc"],
    "auprc": lambda d, r, raw, pos: ranking_metrics(d, pos)["auprc"],
}


def evaluate(scores, subset=None):
    m = table.merge(scores[["chrom", "p", "q", "delta"]], on=["chrom", "p", "q"], how="left")
    sel = np.ones(len(m), bool) if subset is None else subset
    cols = [m["delta"].to_numpy(), m["partial_r"].to_numpy(), m["raw_r"].to_numpy(), m["positive"].to_numpy()]
    out = {}
    for name, f in METRICS.items():
        point = f(*(c[sel] for c in cols))
        draws = []
        for row in boot[:50 if name.startswith("au") else 300]:
            row = row[sel[row]]
            draws.append(f(*(c[row] for c in cols)))
        out[name] = (point, *np.nanpercentile(draws, [2.5, 97.5]))
    return out


min75 = (table["n_molecules"] >= 75).to_numpy()
results = {}
for s in series:
    label, prefix = s.split("=", 1)
    rows = []
    for path in sorted(glob(prefix + "*"), key=lambda p: int(p[len(prefix):])):
        scores = pd.read_csv(f"{path}/pair_scores.tsv", sep="\t")
        epoch = int(path[len(prefix):]) + 1
        rows.append((epoch, evaluate(scores), evaluate(scores, min75)))
        print(label, epoch, {k: round(v[0], 3) for k, v in rows[-1][1].items()},
              "min75", {k: round(v[0], 3) for k, v in rows[-1][2].items()}, flush=True)
    results[label] = rows

SURFACE, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6"]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED})
fig, axes = plt.subplots(1, 4, figsize=(15, 3.8), facecolor=SURFACE)
ceil_all, ceil_75 = split_half_ceiling(table), split_half_ceiling(table[min75])
panels = [("spearman", "Spearman(score, partial r)"), ("spearman_given_raw", "... given raw r"),
          ("auroc", "AUROC, p < 0.01 positives"), ("auprc", "AUPRC")]
for ax, (key, title) in zip(axes, panels, strict=True):
    for k, (label, rows) in enumerate(results.items()):
        e = np.array([r[0] for r in rows])
        for j, (style, name) in enumerate((("-", "all pairs"), ("--", ">= 75 molecules"))):
            v = np.array([r[1 + j][key] for r in rows])
            ax.errorbar(e + (k - 0.5) * 0.4, v[:, 0], yerr=[v[:, 0] - v[:, 1], v[:, 2] - v[:, 0]],
                        color=COLORS[k], ls=style, lw=1.8, marker="o", ms=4, capsize=2,
                        label=f"{label}, {name}")
    if key == "spearman":
        ax.axhline(ceil_all["spearman_ceiling"], color=INK2, ls=":", lw=1)
        ax.axhline(ceil_75["spearman_ceiling"], color=MUTED, ls=":", lw=1)
        ax.text(ax.get_xlim()[0], ceil_all["spearman_ceiling"], " ceiling, all", color=INK2, fontsize=7, va="bottom")
        ax.text(ax.get_xlim()[0], ceil_75["spearman_ceiling"], " ceiling, >= 75", color=MUTED, fontsize=7, va="bottom")
    if key == "auprc":
        ax.axhline(table["positive"].mean(), color=INK2, ls=":", lw=1)
        ax.text(ax.get_xlim()[0], table["positive"].mean(), " positive rate", color=INK2, fontsize=7, va="bottom")
    if key == "auroc":
        ax.axhline(0.5, color=AXIS, lw=0.8)
    ax.set_title(title, loc="left", fontsize=9.5)
    ax.set_xlabel("Epoch")
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRID, lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
axes[0].legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0.0, 0.78))
fig.text(0.01, 0.0, f"{len(table):,} val CRE pairs (chr8, chr20), {int(table['positive'].sum())} positives; "
         f"{int(min75.sum()):,} with >= 75 shared molecules. Score: loss increase on an element's calls when the "
         "partner's calls are hidden (training noise masking). Error bars: 95% pair bootstrap.",
         color=MUTED, fontsize=7.5)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(out_png, dpi=130, facecolor=SURFACE)
print("wrote", out_png)
