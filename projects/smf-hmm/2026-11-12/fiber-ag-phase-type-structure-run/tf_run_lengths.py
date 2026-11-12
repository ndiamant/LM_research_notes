"""TF call lengths and TF-open-TF chaining on the val benchmark windows.

Usage: tf_run_lengths.py MODEL...   (posteriors/MODEL.h5 on chr8 + chr20)

Per read, positions with P(TF) > 0.5 form TF runs. A chain merges TF runs separated by gaps of
at most MAX_GAP bp in which no position has P(histone) > 0.5 (or P(other_long) > 0.5). The
capped TF state (40 bp) can only exceed 40 bp of continuous "TF-like" protection by chaining.
Reports, per model: TF-called bases per read-base, the run-length distribution, and the share of
TF-called bases in runs >= 35 bp and in chains > 40 / > 60 bp.
"""

import sys

import h5py
import numpy as np
import pandas as pd

ROOT = "/scratch/users/diamant/gm12878_validation/covariation_windows_v2"
MAX_GAP = 10
MAX_READS = 60000  # per chromosome, first reads in window order


def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    d = np.diff(np.concatenate(([0], mask.astype(np.int8), [0])))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


rows = []
for model in sys.argv[1:]:
    lengths, chains, tf_bases, total = [], [], 0, 0
    for chrom in ("chr20", "chr8"):
        with h5py.File(f"{ROOT}/{chrom}-512-valid90-reads10-flank128/posteriors/{model}.h5") as f:
            names = list(f["group_names"].asstr()[:])
            tf = [i for i, n in enumerate(names) if n == "TF" or n.startswith("TF_")]
            block = [i for i, n in enumerate(names) if n in ("histone", "other_long")]
            n = min(MAX_READS, f["posterior"].shape[0])
            for a in range(0, n, 4096):
                p = f["posterior"][a:min(a + 4096, n)]
                p_tf = p[..., tf].sum(-1)
                p_block = p[..., block].max(-1) if block else np.zeros_like(p_tf)
                for r in range(len(p)):
                    is_tf = p_tf[r] > 0.5
                    tf_bases += int(is_tf.sum())
                    total += is_tf.size
                    rr = runs(is_tf)
                    lengths += [e - s for s, e in rr]
                    if not rr:
                        continue
                    cs, ce = rr[0]
                    for s, e in rr[1:]:
                        if s - ce <= MAX_GAP and not (p_block[r, ce:s] > 0.5).any():
                            ce = e
                        else:
                            chains.append(ce - cs)
                            cs, ce = s, e
                    chains.append(ce - cs)
    L, C = np.array(lengths), np.array(chains)
    rows.append(dict(
        model=model, tf_base_share=tf_bases / total, n_runs=len(L),
        run_median=np.median(L) if len(L) else np.nan,
        run_p90=np.percentile(L, 90) if len(L) else np.nan,
        bases_in_runs_ge35=L[L >= 35].sum() / max(L.sum(), 1),
        chain_median=np.median(C) if len(C) else np.nan,
        chain_p90=np.percentile(C, 90) if len(C) else np.nan,
        bases_in_chains_gt40=C[C > 40].sum() / max(C.sum(), 1),
        bases_in_chains_gt60=C[C > 60].sum() / max(C.sum(), 1),
        bases_in_chains_gt100=C[C > 100].sum() / max(C.sum(), 1),
    ))
pd.set_option("display.width", 250)
print(f"reads: first {MAX_READS:,} per chromosome; chains merge TF runs across <= {MAX_GAP} bp non-histone gaps")
print(pd.DataFrame(rows).round(3).to_string(index=False))
