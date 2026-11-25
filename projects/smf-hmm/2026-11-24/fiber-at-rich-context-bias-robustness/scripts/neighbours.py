"""Context of MAP TF calls: nearest non-open segment on each side and the kernel evidence pattern (2026-11-24).
Reads calls.tsv.gz (off_motif_calls.py) and the Viterbi tables."""
import numpy as np, pandas as pd
T = "/scratch/users/diamant/fiber_eval_work/tf_mechanism/viterbi"
MODELS = {"r5_tf60_ol61": "AlphaGenome", "r6_ol61_const": "sequence-free"}
C = pd.read_csv("calls.tsv.gz", sep="\t")
out = []
for m, name in MODELS.items():
    s = pd.read_csv(f"{T}/viterbi_segments_{m}_x1.tsv.gz", sep="\t").sort_values(["chrom", "window_id", "read", "start"]).reset_index(drop=True)
    key = ["chrom", "window_id", "read"]
    prot = s[s.state != "open"].copy()
    prot["prev_state"] = prot.groupby(key).state.shift(1); prot["prev_end"] = prot.groupby(key).start.shift(1) + prot.groupby(key).length.shift(1)
    prot["next_state"] = prot.groupby(key).state.shift(-1); prot["next_start"] = prot.groupby(key).start.shift(-1)
    t = prot[prot.state == "TF"].copy()
    t["gap_prev"] = t.start - t.prev_end; t["gap_next"] = t.next_start - (t.start + t.length)
    t = t.merge(C[C.model == name], on=key + ["start", "length"])
    out.append(t)
T2 = pd.concat(out)
T2["between_histones"] = (T2.prev_state == "histone") & (T2.next_state == "histone") & (T2.gap_prev <= 30) & (T2.gap_next <= 30)
T2["next_to_histone"] = ((T2.prev_state == "histone") & (T2.gap_prev <= 30)) | ((T2.next_state == "histone") & (T2.gap_next <= 30))
T2["no_protection_inside"] = T2.kernel_in < 0.5
T2["protected_flank"] = (T2.kernel_left > 0.5) | (T2.kernel_right > 0.5)
T2["flank_unknown"] = T2.kernel_left.isna() | T2.kernel_right.isna()
T2["grp"] = np.where(T2.on_strong, "on strong motif", np.where(T2.length <= 20, "off, <=20 bp", "off, >20 bp"))
cols = ["between_histones", "next_to_histone", "no_protection_inside", "protected_flank", "flank_unknown", "kernel_fp", "access"]
pd.set_option("display.width", 220)
print(T2.groupby(["model", "grp"])[cols].agg(lambda v: v.median() if v.name == "access" else v.mean()).assign(n=T2.groupby(["model", "grp"]).size()).round(3).to_string())
