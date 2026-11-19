"""Figures for the CTCF footprint-position analysis (footprint_position.py outputs; chr20).

Usage: plot_position.py OUT_DIR [TAG=_base] [NULL_TAG=_jit05]
1. examples.png: for 3 high-ChIP sites, each footprint read's kernel P(bound) in motif coordinates (rows sorted by
   footprint centre), next to the same reads re-simulated once under the null (fixed consensus footprint).
2. spread.png: pooled distribution of each footprint's centre minus its site's median, observed vs null; and
   per-site SD of the centre, observed vs null.
"""
import os, sys
import h5py, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, "/home/users/diamant/repos/smf_hmm")
from evaluation.tf_chip.kernel_pbound import call_llr, kernel_pbound, kmer_codes, load_reference
from evaluation.tf_chip.plot_model_comparison import DEFAULT_KERNEL_REFERENCE, read_efficiencies

HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = sys.argv[1]; os.makedirs(OUTD, exist_ok=True)
TAG = os.environ.get("TAG", "_base"); NULL_TAG = os.environ.get("NULL_TAG", "_jit05")
V = "/scratch/users/diamant/gm12878_validation/covariation_windows_v2/chr20-512-valid90-reads10-flank128"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e6e3", "#fcfcfb"
OBS, NUL, NUL2 = "#2a78d6", "#eb6834", "#1baf7a"
KC = LinearSegmentedColormap.from_list("k", ["#eef7f2", "#138a5f", "#06402b"]); KC.set_bad("#ffffff")
BETA, F = load_reference(DEFAULT_KERNEL_REFERENCE)
EFF = h5py.File("/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/active_efficiency.h5")
O = pd.read_csv(f"{HERE}/reads{TAG}.tsv.gz", sep="\t")
N = pd.read_csv(f"{HERE}/null{TAG}.tsv.gz", sep="\t")
N2 = pd.read_csv(f"{HERE}/null{NULL_TAG}.tsv.gz", sep="\t")
site = pd.read_csv(f"{HERE}/sites{TAG}.tsv.gz", sep="\t").dropna(subset=["sd_null"])


def style(ax):
    ax.set_facecolor(SURF); ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)


# --- 2. spread
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), facecolor=SURF)
ax = axes[0]; style(ax)
bins = np.arange(-40.5, 41.5, 1)
dev = lambda D: (D.centre - D.groupby("site").centre.transform("median"))[D.site.isin(site.site)]
for D, lab, col in ((O, "observed reads", OBS), (N, "null: fixed footprint", NUL), (N2, "null + b_m jitter (SD 0.5)", NUL2)):
    v = dev(D)
    ax.hist(v, bins=bins, density=True, histtype="step", lw=2, color=col,
            label=f"{lab} (n = {len(v):,}; |dev| median {np.median(np.abs(v)):.1f} bp)")
ax.set_xlabel("footprint centre − site median (bp, motif strand)", color=INK2, fontsize=9); ax.set_ylabel("density", color=INK2)
ax.set_title("CTCF footprint centres scatter more than call sparsity explains", loc="left", fontsize=10, color=INK)
ax.legend(fontsize=7.5, frameon=False, loc="upper right")
ax = axes[1]; style(ax)
s2 = pd.read_csv(f"{HERE}/sites{NULL_TAG}.tsv.gz", sep="\t").set_index("site")
site = site.set_index("site")
ax.scatter(site.sd_null, site.sd_obs, s=18, color=OBS, edgecolor=SURF, linewidth=0.8, label="null: fixed footprint", zorder=3)
j = site.join(s2[["sd_null"]], rsuffix="_jit").dropna(subset=["sd_null_jit"])
ax.scatter(j.sd_null_jit, j.sd_obs, s=18, color=NUL2, edgecolor=SURF, linewidth=0.8, label="null + b_m jitter", zorder=3)
m = max(site.sd_obs.max(), site.sd_null.max()) * 1.05
ax.plot([0, m], [0, m], color=INK2, lw=1, ls=(0, (4, 3))); ax.set_xlim(0, m); ax.set_ylim(0, m)
ax.set_xlabel("per-site SD of footprint centre under the null (bp)", color=INK2, fontsize=9)
ax.set_ylabel("observed per-site SD (bp)", color=INK2, fontsize=9)
ax.set_title(f"Per site ({len(site)} sites with ≥ 5 footprint reads; dashed = equal)", loc="left", fontsize=10, color=INK)
ax.legend(fontsize=7.5, frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(f"{OUTD}/spread.png", dpi=110, facecolor=SURF); plt.close(fig)

# --- 1. examples: re-derive kernels for 3 sites
S = pd.read_csv(f"{V}/analysis_idr_r6/unique_motif_summaries.tsv.gz", sep="\t").set_index("site_key")
ex = site.sort_values("n_fp", ascending=False).index[:3]
rng = np.random.default_rng(1)
fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), facecolor=SURF, sharex=True)
with h5py.File(f"{V}/dataset.h5") as f:
    for col, key in enumerate(ex):
        s = S.loc[key]
        if isinstance(s, pd.DataFrame):
            s = s.iloc[0]
        w = f["windows"][s.window_id]; obs = w["observations"][:]; dna = w["dna"][:].astype(np.int8); L = len(dna)
        x = np.full((len(obs), L), -1, np.int8); x[:, w["positions"][:]] = obs
        b = read_efficiencies(EFF, w["molecule_id"].asstr()[:])
        k = kernel_pbound(call_llr(x, dna, b, BETA, F, 0.02), 7)
        m0, m1 = int(s.start - s.window_start), int(s.end - s.window_start); c = (m0 + m1 - 1) / 2
        rd = O[O.site == key].sort_values("centre")
        reads = rd.read.to_numpy()
        lo_c, hi_c = int(np.median(rd.left)), int(np.median(rd.right))
        # null: one simulation per footprint read with the consensus footprint
        codes = kmer_codes(dna); sim = np.full((len(reads), L), np.nan)
        cl, ch = (c + lo_c, c + hi_c) if s.strand == "+" else (c - hi_c, c - lo_c)
        for i, r in enumerate(reads):
            called = np.flatnonzero((x[r] >= 0) & (codes >= 0))
            near = np.abs(called - c) <= 100; called = called[near]
            p = np.where((called >= cl) & (called <= ch), 0.02, 1 / (1 + np.exp(-(BETA + F[codes[called]] + b[r]))))
            xs = x[r:r + 1].copy(); xs[0, called] = np.where(rng.random(len(called)) < p, 0, 1)
            sim[i] = kernel_pbound(call_llr(xs, dna, b[r:r + 1], BETA, F, 0.02), 7)[0]
        span = np.arange(int(c) - 80, int(c) + 81)
        for row, (mat, lab) in enumerate(((k[reads][:, span], "observed reads"), (sim[:, span], "null (fixed footprint)"))):
            if s.strand == "-":
                mat = mat[:, ::-1]
            if row == 1:  # sort null rows by their own detected centre for a fair visual comparison
                with np.errstate(all="ignore"):
                    cen = np.array([np.nanmean(np.flatnonzero(v > 0.5)) if (v > 0.5).any() else 0 for v in mat])
                mat = mat[np.argsort(cen)]
            ax = axes[row, col]
            ax.imshow(mat, aspect="auto", interpolation="nearest", cmap=KC, vmin=0, vmax=1, extent=[-80.5, 80.5, len(mat) - 0.5, -0.5])
            mo0, mo1 = (m0 - c, m1 - 1 - c) if s.strand == "+" else (c - (m1 - 1), c - m0)
            ax.axvline(mo0 - 0.5, color="#eda100", lw=1.2); ax.axvline(mo1 + 0.5, color="#eda100", lw=1.2)
            ax.set_title(f"{key.split(':', 2)[2]} — {lab} (n = {len(mat)})", loc="left", fontsize=8.5, color=INK)
            ax.set_ylabel("read (sorted by footprint centre)" if col == 0 else "", color=INK2, fontsize=8)
            ax.tick_params(colors=INK2, labelsize=8)
        axes[1, col].set_xlabel("bp from motif centre (motif strand; yellow = motif)", color=INK2, fontsize=8)
fig.suptitle("Debiased kernel P(bound) of CTCF footprint reads (dark = protected): observed vs simulated with one fixed footprint",
             x=0.01, ha="left", fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(f"{OUTD}/examples.png", dpi=100, facecolor=SURF); plt.close(fig)
print("wrote", OUTD)
