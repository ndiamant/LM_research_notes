"""Prior P(open), P(TF), P(histone) of both models at the AlphaGenome model's MAP TF calls (2026-11-24)."""
import sys
import h5py, numpy as np, pandas as pd
sys.path.insert(0, "/home/users/diamant/repos/smf_hmm")
from evaluation.tf_chip.plot_model_comparison import group_marginals
V = "/scratch/users/diamant/gm12878_validation/covariation_windows_v2"
C = pd.read_csv("calls.tsv.gz", sep="\t")
ag = C[C.model == "AlphaGenome"].copy()
ag["grp"] = np.where(ag.on_strong, "on strong motif", np.where(ag.length <= 20, "off, <=20 bp", "off, >20 bp"))
rows = []
for ds, d in ag.groupby("chrom"):
    pf = {m: h5py.File(f"{V}/{ds}-512-valid90-reads10-flank128/posteriors/{m}_prior.h5") for m in ("r5_tf60_ol61", "r6_ol61_const")}
    for wid, g in d.groupby("window_id"):
        P = {m: group_marginals(f, wid)[0] for m, f in pf.items()}
        for r in g.itertuples():
            s, e = r.start, r.start + r.length
            rows.append(dict(grp=r.grp, **{f"{lab}_{grp}": P[m][s:e, i].mean() for m, lab in (("r5_tf60_ol61", "AG"), ("r6_ol61_const", "SF"))
                                         for i, grp in enumerate(["open", "TF", "other_long", "histone"]) if grp != "other_long"}))
R = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print(R.groupby("grp").median().round(4).to_string())
