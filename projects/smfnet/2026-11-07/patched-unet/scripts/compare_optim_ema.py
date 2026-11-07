"""Val active-region metrics for all models: medians with the pseudo-replicate
ceiling, and paired per-region differences vs the resampled baseline."""
import numpy as np
import pandas as pd

M = "/scratch/users/diamant/smf_models"
runs = {
    "baseline": f"{M}/SMFNet_full_run_2026-08-11/metrics_val_ac_epoch377_current_code_plots",
    "ag_frozen": f"{M}/ONT_mESC_alphagenome_frozen/metrics_val_ac_epoch377_plots",
    "unet_attn": f"{M}/ONT_mESC_alphagenome_unet_attn/metrics_val_ac_epoch366_plots",
    "attn_muon_ema": f"{M}/ONT_mESC_alphagenome_unet_attn_muon_ema/metrics_val_ac_epoch377_plots",
    "attn_adamw_ema": f"{M}/ONT_mESC_alphagenome_unet_attn_adamw_ema/metrics_val_ac_epoch396_plots",
    "atac": f"{M}/ONT_mESC_ATAC_TKO_pooled/metrics_val_ac_epoch377_plots",
}
metrics = {  # name: higher is better
    "pearson_matched": True,
    "pearson_flank_residual_matched": True,
    "covariation_corr_of_corr_pearson": True,
    "covariation_residual_corr_of_corr_pearson": True,
    "run_length_zeros_wasserstein": False,
    "run_length_ones_wasserstein": False,
}
tabs = {}
for k, v in runs.items():
    try:
        tabs[k] = pd.read_csv(f"{v}/metrics.tsv", sep="\t").set_index("region")
    except FileNotFoundError:
        print("missing", k)
base = tabs["baseline"]
rows = []
for m in metrics:
    row = {"metric": m}
    for k, t in tabs.items():
        row[k] = t[m].median() if m in t else float("nan")  # residual metrics need the baselines
    rep = base[f"{m}_pseudo_replicate"].median()
    row["ceiling"] = rep
    for k in tabs:
        if k != "baseline":
            row[f"{k} %gap"] = 100 * (row[k] - row["baseline"]) / (rep - row["baseline"])
    rows.append(row)
pd.set_option("display.width", 250)
print("Medians (7,559 val active regions); %gap = share of baseline-to-ceiling gap closed")
print(pd.DataFrame(rows).set_index("metric").round(4).to_string())
print("\nPaired per-region differences vs baseline: mean (SE) [share better]")
rows = []
for m, higher in metrics.items():
    row = {"metric": m}
    for k, t in tabs.items():
        if k == "baseline":
            continue
        if m not in t:
            row[k] = "n/a (no baselines)"
            continue
        d = (t[m] - base[m]).dropna()
        better = ((d > 0) if higher else (d < 0)).mean()
        row[k] = f"{d.mean():+.4f} ({d.std(ddof=1) / np.sqrt(len(d)):.4f}) [{better:.0%}]"
    rows.append(row)
print(pd.DataFrame(rows).set_index("metric").to_string())

a, b = tabs["attn_adamw_ema"], tabs["attn_muon_ema"]
print("\nPaired per-region, AdamW+EMA minus Muon+EMA: mean (SE) [share where AdamW+EMA is better]")
for m, higher in metrics.items():
    d = (a[m] - b[m]).dropna()
    better = ((d > 0) if higher else (d < 0)).mean()
    print(f"  {m:45s} {d.mean():+.4f} ({d.std(ddof=1) / np.sqrt(len(d)):.4f}) [{better:.0%}]")
