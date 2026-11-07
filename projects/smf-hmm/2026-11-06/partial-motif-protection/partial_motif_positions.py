"""Per-call detail for partial_motif_calls.py: protection by position within the motif.

    PYTHONPATH=<smf_hmm> python partial_motif_positions.py OUT_DIR DATASET_DIR [DATASET_DIR ...]

Writes OUT_DIR/calls.tsv.gz: one row per (non-nucleosomal read, motif call) for
motifs with >=2 CGC-merged calls, positions oriented 5'->3' on the motif strand.
"""
import sys
import h5py, numpy as np, pandas as pd
from evaluation.tf_chip.analyze_covariation import read_table
from smf_hmm.data import merge_cgc_sparse

OUT, datasets = sys.argv[1], sys.argv[2:]
FLANK = (5, 30)
out = []
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
                keep = np.flatnonzero(nonnuc & ((m >= 0).sum(1) >= 2))
                rel = pos[inm] - l
                if r.strand == "-":
                    rel = (rr - l - 1) - rel
                for i in keep:
                    o = m[i] >= 0
                    n_prot = int((m[i] == 1).sum()); n_o = int(o.sum())
                    for j in np.flatnonzero(o):
                        out.append((r.tf, r.site_key, r.chip_pct, rr - l, int(rel[j]), int(m[i, j]), int(i), n_o, n_prot))
c = pd.DataFrame(out, columns=["tf", "site_key", "chip_pct", "motif_width", "rel_pos", "protected", "read", "n_calls", "n_prot"])
c.to_csv(f"{OUT}/calls.tsv.gz", sep="\t", index=False)
print(len(c), "calls")
