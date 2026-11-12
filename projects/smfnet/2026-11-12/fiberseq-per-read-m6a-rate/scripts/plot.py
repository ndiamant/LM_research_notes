import glob, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = [np.load(f) for f in sorted(glob.glob("/scratch/users/diamant/fiber_read_m6a/chr*.npz"))]
acc = np.concatenate([x["acc"] for x in d]).astype(float); pro = np.concatenate([x["pro"] for x in d]); span = np.concatenate([x["span"] for x in d])
n_at = acc + pro; ok = n_at > 0; rate = acc[ok] / n_at[ok]
q = np.percentile(rate, [1, 5, 25, 50, 75, 95, 99])
print(f"molecules {len(acc):,}, no A/T calls {(~ok).sum():,}; median span {np.median(span):.0f} bp, median A/T calls {np.median(n_at):.0f}")
print("rate percentiles 1/5/25/50/75/95/99:", np.round(q, 4), "mean", rate.mean().round(4))
for t in [0.01, 0.02, 0.05]: print(f"rate < {t}: {(rate < t).mean():.2%}")
for t in [0.3, 0.5, 0.95]: print(f"rate > {t}: {(rate > t).mean():.2%}  (n={(rate > t).sum():,})")
pooled = acc.sum() / n_at.sum(); print("pooled base-level rate", round(pooled, 4))

INK, INK2, BLUE, SURF = "#0b0b0b", "#52514e", "#2a78d6", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})
fig, axes = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True, facecolor=SURF, height_ratios=[3, 2])
bins = np.linspace(0, 1, 201)
for ax, log in zip(axes, [False, True]):
    ax.set_facecolor(SURF)
    ax.hist(rate, bins=bins, color=BLUE, edgecolor=SURF, linewidth=0.4, log=log)
    for t in [0.05, 0.95]: ax.axvline(t, color=INK2, lw=1, ls="--")
    ax.axvline(np.median(rate), color=INK, lw=1)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", color="#e4e3df", lw=0.6); ax.set_axisbelow(True)
axes[0].set_ylabel("Molecules"); axes[1].set_ylabel("Molecules (log)"); axes[1].set_xlabel("Per-molecule m6A rate  (accessible / (accessible + protected) A/T calls)")
axes[0].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v/1e3:.0f}k"))
ymax = axes[0].get_ylim()[1]
axes[0].text(0.93, ymax*0.97, f"ONT filter cutoffs (0.05, 0.95)\n{(rate<0.05).mean():.1%} below 0.05, {(rate>0.95).sum()} above 0.95", va="top", ha="right", fontsize=9, color=INK2)
axes[0].text(np.median(rate), ymax*0.62, f" median {np.median(rate):.3f}", fontsize=9, color=INK)
fig.suptitle(f"GM12878 Fiber-seq: per-molecule m6A rate ({ok.sum()/1e6:.1f}M molecules, all chromosomes)", x=0.02, ha="left", fontsize=12, color=INK)
fig.tight_layout(); fig.savefig("/scratch/users/diamant/fiber_read_m6a/read_m6a_rate_hist.png", dpi=150, facecolor=SURF)
