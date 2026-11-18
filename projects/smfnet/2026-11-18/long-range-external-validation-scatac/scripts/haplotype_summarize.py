import sys, numpy as np, pandas as pd
t = pd.read_csv(sys.argv[1], sep="\t")
t["n"] = t.n_hap1 + t.n_hap2
def block(d, label):
    out = {"group": label, "elements": len(d), "median n/hap": int(np.median(np.minimum(d.n_hap1, d.n_hap2)))}
    for m in ("b", "bl"):
        p, p0 = d[f"p_{m}"], d[f"pnull_{m}"]
        out[f"p<.001 {m}"] = f"{(p<1e-3).sum()} (null {(p0<1e-3).sum()})"
        out[f"q<.05 {m}"] = f"{(d[f'q_{m}']<.05).sum()} (null {(d[f'qnull_{m}']<.05).sum()})"
        lam = 0.5; out[f"1-pi0 {m}"] = round(1 - (p > lam).mean() / (1 - lam), 3)
    s = d[d.q_b < .05]
    out["|delta| sig median"] = round(s.delta_b.abs().median(), 3) if len(s) else np.nan
    out["|delta|>=0.1 & q<.05"] = int(((d.q_b < .05) & (d.delta_b.abs() >= .1)).sum())
    return out
rows = []
for chrom, d in t.groupby("chrom"): rows.append(block(d, chrom))
a = t[t.chrom != "chrX"].copy()
if "nearest_het" in a:
    a = a[a.confident_1kb >= 0.9]
    bins = [-1, 100, 500, 1000, 2000, 1e9]
    labels = ["het in element", "100-500", "500-1k", "1-2k", ">2k"]
    a["bin"] = pd.cut(a.nearest_het, bins, labels=labels)
    for b_, d in a.groupby("bin", observed=True): rows.append(block(d, f"auto, nearest het {b_}"))
print(pd.DataFrame(rows).to_string(index=False))
print("\nnull p uniformity (autosomes, model b): KS vs U(0,1):")
from scipy import stats
print(stats.kstest(t[t.chrom!='chrX'].pnull_b.dropna(), 'uniform'))
print("\nlevel_diff (domain shift) |median|: auto", t[t.chrom!='chrX'].level_diff.abs().median().round(4), " chrX", t[t.chrom=='chrX'].level_diff.abs().median().round(4))
g = t[(t.chrom=="chr20")&(t.centre.between(58_800_000, 58_950_000))]
print("\nGNAS region (chr20:58.80-58.95 Mb):"); print(g[["centre","n_hap1","n_hap2","acc_hap1","acc_hap2","delta_b","p_b","q_b","nearest_het"]].round(4).to_string(index=False))
