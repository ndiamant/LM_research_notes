"""TF share by context-efficiency decile (outside the strong-footprint mask): prior vs posterior vs MAP, and
the posterior / prior ratio (what the reads add), for both option-B models (2026-11-24)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
c = pd.read_csv("context_deciles.tsv", sep="\t")
c = c[(c.exclusion == "strong") & (c.element == "TF")]
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
col = {"AlphaGenome": "#1f6fb4", "sequence-free": "#d1495b"}; ls = {"prior": ":", "posterior": "-", "MAP": "--"}
for m in col:
    for lvl in ls:
        t = c[(c.model == m) & (c.level == lvl)].sort_values("bin")
        ax = axes[0] if m == "AlphaGenome" else axes[1]
        ax.plot(t.bin, t.share * 1000, ls[lvl], color=col[m], marker="o", ms=3, label=lvl)
        ax.fill_between(t.bin, t.lo * 1000, t.hi * 1000, color=col[m], alpha=0.12)
    p = c[(c.model == m) & (c.level == "posterior")].sort_values("bin").share.to_numpy()
    q = c[(c.model == m) & (c.level == "prior")].sort_values("bin").share.to_numpy()
    axes[2].plot(range(10), p / q, "-o", ms=3, color=col[m], label=m)
for ax, m in zip(axes[:2], col):
    ax.set_title(f"{m}: TF per 1,000 positions (outside strong motifs)", fontsize=9); ax.legend(fontsize=7)
    ax.set_xlabel("context-efficiency decile (0 = AT-rich, least efficient)"); ax.set_ylim(0, 18)
axes[2].axhline(1, color="k", lw=0.5); axes[2].set_title("posterior / prior TF (what the reads add)", fontsize=9)
axes[2].set_xlabel("context-efficiency decile"); axes[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig("context_tf_prior_vs_posterior.png", dpi=130)
