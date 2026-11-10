"""Is the shared short-range error a sampler effect? Same model and regions
(the 1,000-region subset) at 32 vs 64 sampling steps: pooled mean pairwise
correlation by separation and accessible-run histograms, model - observed."""
import numpy as np
M = "/scratch/users/diamant/smf_models/"
for run, full in [("ONT_mESC_alphagenome_unet_attn_muon_ema", "metrics_val_ac_epoch377_plots"),
                  ("ONT_mESC_alphagenome_unet_attn_muon_ema_k9", "metrics_val_ac_epoch398_plots")]:
    a = np.load(M + run + "/" + full + "/covariation_decay.npz")
    b = np.load(M + run + "/metrics_val_ac_subset1000_steps64/covariation_decay.npz")
    ra = np.load(M + run + "/" + full + "/run_length_histograms.npz")
    rb = np.load(M + run + "/metrics_val_ac_subset1000_steps64/run_length_histograms.npz")
    idx = {r: i for i, r in enumerate(a["region_ids"])}
    sel = np.array([idx[r] for r in b["region_ids"]])
    mean = lambda d, w, s=slice(None): d[f"covariation_decay_sum_{w}"][s].sum(0) / np.maximum(d[f"covariation_decay_count_{w}"][s].sum(0), 1)
    o = mean(b, "observed"); e = a["bin_edges"]
    print(f"\n{run.split('muon_ema')[-1] or '_k5'}: correlation error by separation, same 1,000 regions")
    print(f"{'sep bp':>9s} {'observed':>9s} {'32 steps':>9s} {'64 steps':>9s}")
    for i in [0, 1, 2, 3, 4, 6, 8, 10, 20]:
        print(f"{e[i]:4d}-{e[i + 1]:<4d} {o[i]:9.4f} {mean(a, 'model', sel)[i] - o[i]:+9.4f} {mean(b, 'model')[i] - o[i]:+9.4f}")
    for kind in ("ones", "zeros"):
        ob = rb[f"run_length_hist_{kind}_observed"].sum(0); ob = ob / ob.sum()
        m32 = ra[f"run_length_hist_{kind}_model"][sel].sum(0); m32 = m32 / m32.sum()
        m64 = rb[f"run_length_hist_{kind}_model"].sum(0); m64 = m64 / m64.sum()
        print(f"  {kind} runs: share 0-10 bp observed {ob[0]:.4f}, 32 steps {m32[0] - ob[0]:+.4f}, 64 steps {m64[0] - ob[0]:+.4f};"
              f" total variation 32 {0.5 * abs(m32 - ob).sum():.4f}, 64 {0.5 * abs(m64 - ob).sum():.4f}")
