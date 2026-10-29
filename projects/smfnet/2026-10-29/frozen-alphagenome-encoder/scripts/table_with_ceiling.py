import pandas as pd
runs = {
    "baseline": "/scratch/users/diamant/smf_models/SMFNet_full_run_2026-08-11/metrics_val_ac_epoch377_current_code_plots",
    "alphagenome": "/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_frozen/metrics_val_ac_epoch377_plots",
    "atac": "/scratch/users/diamant/smf_models/ONT_mESC_ATAC_TKO_pooled/metrics_val_ac_epoch377_plots",
}
tabs = {k: pd.read_csv(f"{v}/metrics.tsv", sep="\t") for k, v in runs.items()}
metrics = ["pearson_matched", "pearson_flank_residual_matched", "covariation_corr_of_corr_pearson",
           "covariation_residual_corr_of_corr_pearson", "run_length_zeros_wasserstein", "run_length_ones_wasserstein"]
b = tabs["baseline"]
print("check summary stat for pearson_matched: mean", b.pearson_matched.mean(), "median", b.pearson_matched.median())
rows = []
for m in metrics:
    row = {"metric": m}
    for k, t in tabs.items():
        row[k + "_mean"] = t[m].mean(); row[k + "_median"] = t[m].median(); row[k + "_n"] = t[m].notna().sum()
        rep = f"{m}_pseudo_replicate"
        if rep in t:
            row[k + "_rep_mean"] = t[rep].mean(); row[k + "_rep_median"] = t[rep].median()
    rows.append(row)
pd.set_option("display.width", 250)
print(pd.DataFrame(rows).set_index("metric").T.round(4).to_string())
