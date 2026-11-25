"""What the sequence-free model calls at the AlphaGenome model's off-motif MAP TF calls, same reads (2026-11-24)."""
import numpy as np, pandas as pd
T = "/scratch/users/diamant/fiber_eval_work/tf_mechanism/viterbi"
C = pd.read_csv("calls.tsv.gz", sep="\t")
sf = pd.read_csv(f"{T}/viterbi_segments_r6_ol61_const_x1.tsv.gz", sep="\t")
idx = {k: g for k, g in sf.groupby(["chrom", "window_id", "read"])}
ag = C[(C.model == "AlphaGenome")].copy()
ag["grp"] = np.where(ag.on_strong, "on strong motif", np.where(ag.length <= 20, "off, <=20 bp", "off, >20 bp"))
rows = []
for r in ag.itertuples():
    g = idx[(r.chrom, r.window_id, r.read)]
    a, b = r.start, r.start + r.length
    ov = np.minimum(g.start + g.length, b) - np.maximum(g.start, a)
    cover = g.assign(ov=ov.clip(lower=0)).groupby("state").ov.sum()
    rows.append(dict(grp=r.grp, majority=cover.idxmax(), **{f"bp_{s}": cover.get(s, 0) / r.length for s in ("open", "TF", "other_long", "histone")}))
R = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("share of AG call covered by each sequence-free state:")
print(R.groupby("grp")[["bp_open", "bp_TF", "bp_other_long", "bp_histone"]].mean().round(3).to_string())
print("\nmajority sequence-free state:")
print(pd.crosstab(R.grp, R.majority, normalize="index").round(3).to_string())
