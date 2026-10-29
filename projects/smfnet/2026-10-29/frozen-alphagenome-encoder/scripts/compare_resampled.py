import sys
import numpy as np
import pandas as pd

runs = {
    "baseline": "/scratch/users/diamant/smf_models/SMFNet_full_run_2026-08-11/metrics_val_ac_epoch377_current_code_plots",
    "baseline_old_samples": "/scratch/users/diamant/evaluate_SMFNet_full_run_2026-08-11/metrics_val_ac_current_code_plots",
    "alphagenome": "/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_frozen/metrics_val_ac_epoch377_plots",
    "atac": "/scratch/users/diamant/smf_models/ONT_mESC_ATAC_TKO_pooled/metrics_val_ac_epoch377_plots",
}
pd.set_option("display.width", 200)
summ = {k: pd.read_csv(f"{v}/metrics_summary.tsv", sep="\t").set_index("metric") for k, v in runs.items()}
print("Summary (model value; fraction = share of the null-to-floor gap closed)")
print(pd.concat({k: s[["model", "fraction"]] for k, s in summ.items()}, axis=1).round(4).to_string())
print("floor/null:", summ["baseline"][["floor", "null"]].round(4).to_dict("index"))

metrics = {  # name: higher is better
    "pearson_matched": True,
    "pearson_flank_residual_matched": True,
    "run_length_zeros_wasserstein": False,
    "run_length_ones_wasserstein": False,
    "covariation_corr_of_corr_pearson": True,
    "covariation_residual_corr_of_corr_pearson": True,
}
tabs = {k: pd.read_csv(f"{v}/metrics.tsv", sep="\t").set_index("region") for k, v in runs.items()}
base = tabs["baseline"]
print("\nPer-region paired differences vs baseline (mean, SE, median, share of regions better)")
rows = []
for other in ("baseline_old_samples", "alphagenome", "atac"):
    t = tabs[other]
    for m, higher in metrics.items():
        d = (t[m] - base[m]).dropna()
        better = (d > 0) if higher else (d < 0)
        rows.append((other, m, len(d), d.mean(), d.std(ddof=1) / np.sqrt(len(d)), d.median(), better.mean()))
print(pd.DataFrame(rows, columns=["model", "metric", "n", "mean_diff", "se", "median_diff", "frac_better"]).round(4).to_string(index=False))

sites = {k: pd.read_csv(f"{v}/site_correlations.tsv", sep="\t").set_index("comparison") for k, v in runs.items()}
print("\nPooled site correlations (model_all_reads):")
for k, s in sites.items():
    print(f"  {k}: pearson {s.loc['model_all_reads', 'pearson']:.4f}, spearman {s.loc['model_all_reads', 'spearman']:.4f}")
