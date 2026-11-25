"""Relative accessibility A of held-out calls inside MAP TF calls by held-out position efficiency quintile (2026-11-24)."""
import pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = pd.read_csv("heldout_inside_tf.tsv", sep="\t")
cls = {"TF strong motif": ("#3a9a5b", "-"), "TF off, low-eff call": ("#d1495b", "-"), "TF off, high-eff call": ("#e0a030", "-"),
       "histone": ("#555555", ":")}
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
for ax, m in zip(axes, ["AlphaGenome", "sequence-free"]):
    for c, (color, ls) in cls.items():
        t = R[(R.model == m) & (R.cls == c)].sort_values("qbin")
        ax.plot(t.qbin, t.A, ls, color=color, marker="o", ms=3, label=c)
        ax.fill_between(t.qbin, t.A_lo, t.A_hi, color=color, alpha=0.15)
    ax.axhline(1, color="k", lw=0.6); ax.set_title(f"{m}: held-out calls inside MAP calls", fontsize=9)
    ax.set_xlabel("held-out position context-efficiency quintile (0 = least efficient)"); ax.set_ylim(0, 1.05)
axes[0].set_ylabel("accessibility relative to open calls"); axes[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig("heldout_accessibility_inside_tf.png", dpi=130)
