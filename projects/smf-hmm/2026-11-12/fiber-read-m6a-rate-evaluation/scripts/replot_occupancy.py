"""Redraw occupancy_by_rate.png with one fixed color per state group (analyze_read_rate.py colored
by column order, which differs between models), and print the CTCF/NRF1 nucleosome-flank contrast."""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
D = sys.argv[1]
MODELS = ["phase_type", "ag_kl10", "ag_ol001", "ag_ol001_s05", "ag_kl10_s0"]
COLORS = {"open": "#8a8a86", "TF": "#2a78d6", "other_long": "#eb6834", "histone": "#4a3aa7"}
reads = pd.concat([pd.read_csv(f"{D}/reads_{c}.tsv.gz", sep="\t") for c in ("chr8", "chr20")], ignore_index=True)
bins = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 1.01])
reads["bin"] = pd.cut(reads.mol_rate, bins, right=False)
centers = (bins[:-1] + np.minimum(bins[1:], 0.3)) / 2
fig, axes = plt.subplots(1, len(MODELS), figsize=(3.3 * len(MODELS), 3.2), sharey=True)
for ax, m in zip(axes, MODELS):
    for g, color in COLORS.items():
        col = f"{m}_{g}"
        if col in reads:
            ax.plot(centers, reads.groupby("bin", observed=False)[col].mean(), marker="o", ms=3, color=color, label=g)
    ax.axvspan(0.02, 0.08, color="0.92", zorder=0)
    ax.set_title(m, fontsize=9); ax.set_xlabel("molecule m6A rate")
axes[0].set_ylabel("mean posterior occupancy (512 bp)")
axes[1].legend(fontsize=7, frameon=False)
fig.suptitle("State occupancy by molecule m6A rate, val benchmark reads (gray: Bubb et al. 2-8%)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{D}/occupancy_by_rate.png", dpi=100)
acc = [np.load(f"{D}/accumulators_{c}.npz") for c in ("chr8", "chr20")]
groups = {"low <2%": [0], "2-8%": [1, 2, 3], "8-12%": [4, 5], ">=12%": [6, 7, 8]}
rows = []
for tf in ("ctcf", "nrf1"):
    m6a = sum(a[f"prof_m6a_{tf}_bound"] for a in acc); called = sum(a[f"prof_called_{tf}_bound"] for a in acc)
    x = np.arange(m6a.shape[1]) - m6a.shape[1] // 2
    for g, idx in groups.items():
        r = m6a[idx].sum(0) / np.maximum(called[idx].sum(0), 1)
        linker = r[(np.abs(x) >= 30) & (np.abs(x) <= 50)].mean()
        outer = r[np.abs(x) >= 105].mean()
        rows.append(dict(tf=tf, cls=g, linker_rate=linker, outer_rate=outer, outer_over_linker=outer / linker))
print(pd.DataFrame(rows).round(3).to_string(index=False))
