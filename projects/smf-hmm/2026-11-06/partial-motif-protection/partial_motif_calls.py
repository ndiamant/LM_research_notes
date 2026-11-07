"""Within-read partial protection of motifs with >=2 calls, on non-nucleosomal reads.

    PYTHONPATH=<smf_hmm> python partial_motif_calls.py OUT_DIR DATASET_DIR [DATASET_DIR ...]

Calls are CGC-merged (as in training). A read is "non-nucleosomal" at a motif
if both flanks (5-30 bp from the motif) have >=1 call and all are accessible.
Reads covering >=2 motif calls are classed all-accessible / mixed /
all-protected. The independence expectation uses each site's per-call
protection rates among its non-nucleosomal reads.
"""
import sys
import h5py, numpy as np, pandas as pd
from evaluation.tf_chip.analyze_covariation import read_table
from smf_hmm.data import merge_cgc_sparse

OUT, datasets = sys.argv[1], sys.argv[2:]
FLANK = (5, 30)
site_rows, read_rows = [], []
for D in datasets:
    u = pd.read_csv(f"{D}/analysis_clip_ema/unique_motif_summaries.tsv.gz", sep="\t")
    u["chip_pct"] = u.groupby("tf").chip_pooled_sum.rank(pct=True)
    with h5py.File(f"{D}/dataset.h5") as f:
        win = read_table(f["window_metadata"]).set_index("window_id")
        for wid, rows in u.groupby("window_id"):
            w = f["windows"][wid]
            obs, pos = merge_cgc_sparse(w["dna"][:].astype(np.int64), w["observations"][:], w["positions"][:])
            ws = int(win.loc[wid, "window_start"])
            for r in rows.itertuples():
                l, rr = r.start - ws, r.end - ws
                inm = (pos >= l) & (pos < rr)
                if inm.sum() < 2:
                    continue
                sides = [obs[:, (pos >= lo) & (pos < hi)] for lo, hi in ((l - FLANK[1], l - FLANK[0]), (rr + FLANK[0], rr + FLANK[1]))]
                nonnuc = np.all([(s >= 0).any(1) & ~(s == 1).any(1) for s in sides], axis=0)
                m = obs[:, inm]
                n_obs = (m >= 0).sum(1)
                keep = nonnuc & (n_obs >= 2)
                if keep.sum() == 0:
                    continue
                mk = m[keep].astype(float); mk[mk < 0] = np.nan
                p = np.nanmean(mk, 0)  # per-call protection rate at this site
                n_prot = np.nansum(mk, 1); n_o = (~np.isnan(mk)).sum(1)
                # independence expectation per read, over its observed calls
                pp = np.where(np.isnan(mk), 1.0, np.broadcast_to(np.nan_to_num(p), mk.shape))
                qq = np.where(np.isnan(mk), 1.0, 1 - np.broadcast_to(np.nan_to_num(p), mk.shape))
                e_all_p, e_all_a = pp.prod(1), qq.prod(1)
                site_rows.append(dict(tf=r.tf, site_key=r.site_key, chip=r.chip_pooled_sum, chip_pct=r.chip_pct,
                    n_motif_calls=int(inm.sum()), n_reads=int(keep.sum()),
                    all_acc=(n_prot == 0).mean(), mixed=((n_prot > 0) & (n_prot < n_o)).mean(), all_prot=(n_prot == n_o).mean(),
                    exp_mixed=(1 - e_all_p - e_all_a).mean(), exp_all_prot=e_all_p.mean(), call_prot_rate=np.nanmean(mk)))
s = pd.DataFrame(site_rows)
s.to_csv(f"{OUT}/sites.tsv", sep="\t", index=False)
s["chip_q"] = pd.cut(s.chip_pct, [0, .25, .5, .75, 1], labels=["Q1", "Q2", "Q3", "Q4"])
def summarize(g):
    w = g.n_reads
    avg = lambda c: np.average(g[c], weights=w)
    return pd.Series(dict(sites=len(g), reads=int(w.sum()), call_prot=avg("call_prot_rate"), all_acc=avg("all_acc"),
        mixed=avg("mixed"), all_prot=avg("all_prot"), exp_mixed=avg("exp_mixed"), exp_all_prot=avg("exp_all_prot")))
t = s.groupby(["tf", "chip_q"], observed=True).apply(summarize, include_groups=False)
t["mixed/exp"] = t.mixed / t.exp_mixed
t["all_prot/exp"] = t.all_prot / t.exp_all_prot
t.to_csv(f"{OUT}/summary_by_chip_quartile.tsv", sep="\t")
pd.set_option("display.width", 250)
print(t.round(3).to_string())
