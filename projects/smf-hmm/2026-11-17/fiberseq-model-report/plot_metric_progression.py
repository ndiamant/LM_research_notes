"""Metric progression along the main Fiber-seq model lineage (val chr8 + chr20, prior x1), for the report.

Usage: python plot_metric_progression.py OUT.png
Reads the evaluation tables under /scratch/users/diamant/fiber_eval_work (idr_summary_*_val.tsv,
kernel_ref/specificity*/{auprc_matched,off_motif}.tsv, prior_log_odds/val_*priors.tsv); the first file
that has a model is used (values are identical across files). Writes the figure and a TSV of the plotted
values next to it.
"""
import glob
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

W = Path("/scratch/users/diamant/fiber_eval_work")
LINEAGE = [  # (model key, label)
    ("ag_kl10", "AG KL 10\n(11-10)"), ("ag_ol001", "other_long 1%\n(11-11)"), ("ag_sptf", "spacers +\nTF/3 (11-12)"),
    ("abl_B", "no chaining\n(11-14)"), ("gelu_s0", "GELU head\n(11-15)"), ("fix_kl300", "KL 300\n(11-16)"),
    ("r2_full", "tie + TF 15-40\n+ flank (11-17)"), ("r3_hist100_190", "histone\n100-190"), ("r4_hist190_tf60", "TF 15-60"), ("r5_tf60_ol61", "other_long\n61-100 (11-18)"),
]
BASELINE = "phase_type"


def first(frames, key, col, value):
    for f in frames:
        hit = f[f[key] == value]
        if len(hit):
            return hit
    return None


idr = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(W / "idr_summary_*val*.tsv")))]
idr = [d[d.split == "val"] for d in idr]
auprc = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(W / "kernel_ref/specificity*/auprc_matched.tsv")))]
off = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(W / "kernel_ref/specificity*/off_motif.tsv")))]
prior = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(W / "prior_log_odds/val_*priors.tsv")))]
rows = []
for m in [k for k, _ in LINEAGE] + [BASELINE]:
    r = {"model": m}
    d = first(idr, "model", None, m)
    if d is not None:
        r.update(auroc=d.auroc.mean(), auroc_access=d.auroc_access.mean(), rho_bound=d.rho_bound_access.mean())
    a = first(auprc, "score", None, m)
    r["auprc"] = a.auprc.mean() if a is not None else np.nan
    o = first(off, "model", None, m)
    r["off_motif_kb"] = float(o.tf_off_per_read_kb.iloc[0]) if o is not None else np.nan
    p = first(prior, "prior", None, f"{m}_prior")
    r["prior_log_odds"] = p.log_odds.mean() if p is not None else np.nan
    rows.append(r)
T = pd.DataFrame(rows).set_index("model")
out = Path(sys.argv[1])
T.round(4).to_csv(out.with_suffix(".tsv"), sep="\t")

INK, INK2, GRID, SURF, BLUE = "#0b0b0b", "#52514e", "#e6e6e3", "#fcfcfb", "#2a78d6"
PANELS = [("auroc", "IDR AUROC"), ("auroc_access", "AUROC, accessible reads"), ("rho_bound", "ρ within bound sites"),
          ("auprc", "accessibility-matched AUPRC"), ("off_motif_kb", "off-motif TF per read-kb (lower is better)"),
          ("prior_log_odds", "prior TF log-odds (9-TF mean)")]
x = np.arange(len(LINEAGE))
fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2), sharex=True, facecolor=SURF)
for ax, (col, title) in zip(axes.ravel(), PANELS):
    ax.set_facecolor(SURF)
    y = T.loc[[k for k, _ in LINEAGE], col].to_numpy(float)
    ok = ~np.isnan(y)
    ax.plot(x[ok], y[ok], color=BLUE, lw=2, zorder=3)
    ax.scatter(x[ok], y[ok], s=40, color=BLUE, edgecolor=SURF, linewidth=2, zorder=4)
    ax.scatter(x[-1:], y[-1:], s=90, color=BLUE, edgecolor=SURF, linewidth=2, zorder=5)
    ax.annotate(f"{y[-1]:.3f}" if col != "off_motif_kb" else f"{y[-1]:.1f}", (x[-1], y[-1]), xytext=(6, 6),
                textcoords="offset points", fontsize=9, color=INK)
    b = T.loc[BASELINE, col]
    if col == "off_motif_kb":  # phase-type (38.6) would squash the axis: give its value in the title
        ax.set_ylim(0, 16)
        title += f"\nphase-type {b:.1f}, off scale"
    elif np.isfinite(b):
        ax.axhline(b, color=INK2, lw=1.2, ls=(0, (4, 3)), zorder=2)
        title += f"\ndashed: phase-type {b:.3f}"
    else:
        title += "\nphase-type has no sequence prior"
    if col in ("auprc", "off_motif_kb"):
        ax.text(0.02, 0.96, "not computed before 11-13", transform=ax.transAxes, va="top", fontsize=8, color=INK2)
    ax.set_title(title, loc="left", fontsize=10.5, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5)
for ax in axes[1]:
    ax.set_xticks(x, [l for _, l in LINEAGE], rotation=45, ha="right", fontsize=8)
fig.suptitle("GM12878 Fiber-seq HMM: benchmark metrics along the model lineage (val chr8 + chr20, prior ×1)", x=0.01, ha="left", fontsize=12, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(out, dpi=120, facecolor=SURF)
print(T.round(3).to_string())
