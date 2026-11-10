"""Fine-scale views of the sampled-read metrics, k=9 vs k=5 attention U-Net
(and the original SMFNet for scale), from each run's 32-step metrics output:
(1) mean pairwise correlation by separation (covariation_decay.npz), pooled
over regions, model vs observed; (2) run-length histograms by length bin,
pooled, model vs observed."""
import numpy as np
M = "/scratch/users/diamant/smf_models/"
runs = {
    "original": "SMFNet_full_run_2026-08-11/metrics_val_ac_epoch377_current_code_plots",
    "k5": "ONT_mESC_alphagenome_unet_attn_muon_ema/metrics_val_ac_epoch377_plots",
    "k9": "ONT_mESC_alphagenome_unet_attn_muon_ema_k9/metrics_val_ac_epoch398_plots",
}
cov = {k: np.load(M + v + "/covariation_decay.npz") for k, v in runs.items()}
rl = {k: np.load(M + v + "/run_length_histograms.npz") for k, v in runs.items()}
ids = cov["k5"]["region_ids"]
for k in runs:
    assert (cov[k]["region_ids"] == ids).all() and (rl[k]["region_ids"] == ids).all()
edges = cov["k5"]["bin_edges"]
mean = lambda d, w: d[f"covariation_decay_sum_{w}"].sum(0) / np.maximum(d[f"covariation_decay_count_{w}"].sum(0), 1)
obs = mean(cov["k5"], "observed")
print("Mean pairwise correlation by separation (pooled over 7,559 regions)")
print(f"{'sep bp':>9s} {'observed':>9s} " + " ".join(f"{k + ' err':>12s}" for k in runs))
for i in list(range(0, 10)) + list(range(10, 50, 5)):
    print(f"{edges[i]:4d}-{edges[i + 1]:<4d} {obs[i]:9.4f} " + " ".join(f"{mean(cov[k], 'model')[i] - obs[i]:+12.4f}" for k in runs))
# Per-region paired: |model - observed| mean correlation in separation groups
def region_err(d, lo, hi):
    sl = slice(np.searchsorted(edges, lo), np.searchsorted(edges, hi))
    m = d["covariation_decay_sum_model"][:, sl].sum(1) / np.maximum(d["covariation_decay_count_model"][:, sl].sum(1), 1)
    o = d["covariation_decay_sum_observed"][:, sl].sum(1) / np.maximum(d["covariation_decay_count_observed"][:, sl].sum(1), 1)
    ok = (d["covariation_decay_count_observed"][:, sl].sum(1) > 0) & (d["covariation_decay_count_model"][:, sl].sum(1) > 0)
    return np.abs(m - o), ok
print("\nPer-region |model - observed| mean correlation, k9 - k5 paired: mean (SE)")
for lo, hi in [(0, 20), (20, 50), (50, 150), (150, 500)]:
    a, oka = region_err(cov["k9"], lo, hi); b, okb = region_err(cov["k5"], lo, hi); ok = oka & okb
    d = a[ok] - b[ok]
    print(f"  {lo:3d}-{hi:<3d} bp: k5 {b[ok].mean():.4f}  k9 {a[ok].mean():.4f}  diff {d.mean():+.5f} ({d.std(ddof=1) / np.sqrt(ok.sum()):.5f})")
redges = rl["k5"]["bin_edges"]
print("\nRun-length histograms pooled over regions: share of runs per length bin, model - observed")
for kind in ("zeros", "ones"):
    o = rl["k5"][f"run_length_hist_{kind}_observed"].sum(0); o = o / o.sum()
    print(f"  {kind}: {'bin bp':>9s} {'observed':>9s} " + " ".join(f"{k:>9s}" for k in runs))
    for i in list(range(0, 8)) + list(range(8, 40, 4)):
        row = []
        for k in runs:
            m = rl[k][f"run_length_hist_{kind}_model"].sum(0); m = m / m.sum(); row.append(f"{m[i] - o[i]:+9.4f}")
        print(f"         {redges[i]:4d}-{redges[i + 1]:<4d} {o[i]:9.4f} " + " ".join(row))
    tv = {}
    for k in runs:
        m = rl[k][f"run_length_hist_{kind}_model"].sum(0); m = m / m.sum(); tv[k] = 0.5 * np.abs(m - o).sum()
    print("    total variation distance (pooled):", {k: round(v, 4) for k, v in tv.items()})
