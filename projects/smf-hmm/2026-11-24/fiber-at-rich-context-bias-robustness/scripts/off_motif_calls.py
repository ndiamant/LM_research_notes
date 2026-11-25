"""Why does the AlphaGenome option-B model make more off-motif MAP TF calls than the sequence-free one? (2026-11-24)

For every MAP TF call (Viterbi, every 8th val window, chr8 + chr20) of r5_tf60_ol61 and r6_ol61_const:
  on_strong   overlaps the locked mask (20 strong-footprint PWMs, FPR 1e-6, +-10 bp)
  on_any      overlaps any of the 100 PWMs (91 screen + 9 benchmark, FPR 1e-6, +-10 bp)
  kernel_in   mean debiased kernel P(bound) (+-7 bp) over the call; kernel_flank: mean over 10-30 bp each side
  kernel_fp   footprint evidence: kernel_in > 0.8 and both flanks < 0.2 (as in specificity_metrics' read check)
  prior_tf    the model's mean prior P(TF) over the call
  access      read-pooled m6A fraction of calls within +-100 bp of the call centre (window accessibility)
Writes calls.tsv.gz and prints summaries.
"""
import sys
import h5py, numpy as np, pandas as pd
sys.path.insert(0, "/home/users/diamant/repos/smf_hmm")
from evaluation.tf_chip.specificity_metrics import BENCHMARK_PWMS, DEFAULT_MOTIFS, REFERENCES, meme_scanner, motif_mask
from evaluation.tf_chip.plot_model_comparison import motif_scanners, group_marginals, read_efficiencies, DEFAULT_KERNEL_REFERENCE
from evaluation.tf_chip.kernel_pbound import call_llr, kernel_pbound, load_reference
V = "/scratch/users/diamant/gm12878_validation/covariation_windows_v2"
T = "/scratch/users/diamant/fiber_eval_work/tf_mechanism/viterbi"
MODELS = {"r5_tf60_ol61": "AlphaGenome", "r6_ol61_const": "sequence-free"}
strong = motif_scanners(DEFAULT_MOTIFS, 1e-6, set((REFERENCES / "gm12878_strong_footprint_tfs.txt").read_text().split()))
anyp = motif_scanners(DEFAULT_MOTIFS, 1e-6) + [meme_scanner(REFERENCES / f"{m}.meme", 1e-6) for m in BENCHMARK_PWMS]
EFF = h5py.File("/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/active_efficiency.h5")
REF = load_reference(DEFAULT_KERNEL_REFERENCE)
segs = {m: pd.read_csv(f"{T}/viterbi_segments_{m}_x1.tsv.gz", sep="\t").query("state == 'TF'") for m in MODELS}
rows = []
for ds in ("chr20", "chr8"):
    D = f"{V}/{ds}-512-valid90-reads10-flank128"
    prio = {m: h5py.File(f"{D}/posteriors/{m}_prior.h5") for m in MODELS}
    by = {m: dict(tuple(s[s.chrom == ds].groupby("window_id"))) for m, s in segs.items()}
    wins = sorted(set().union(*[set(b) for b in by.values()]))
    with h5py.File(f"{D}/dataset.h5") as f:
        for wid in wins:
            w = f["windows"][wid]; dna = w["dna"][:].astype(np.int64); obs = w["observations"][:]; L = len(dna); R = len(obs)
            x = np.full((R, L), -1, np.int8); x[:, w["positions"][:]] = obs
            k = kernel_pbound(call_llr(x, dna, read_efficiencies(EFF, w["molecule_id"].asstr()[:]), *REF, 0.02), 7)
            ms, ma = motif_mask(dna, strong, 10), motif_mask(dna, anyp, 10)
            m6a = np.r_[0, np.cumsum((x == 0).sum(0))]; called = np.r_[0, np.cumsum((x >= 0).sum(0))]
            for m in MODELS:
                g = by[m].get(wid)
                if g is None:
                    continue
                ptf = group_marginals(prio[m], wid)[0][:, 1]
                for r, s, n in zip(g.read.to_numpy(), g.start.to_numpy(), g.length.to_numpy()):
                    e = s + n; c = (s + e) // 2
                    lo, hi = max(c - 100, 0), min(c + 100, L)
                    kin = k[r, s:e]; kl = k[r, max(s - 30, 0):max(s - 10, 0)]; kr = k[r, min(e + 10, L):min(e + 30, L)]
                    with np.errstate(all="ignore"):
                        rows.append(dict(model=MODELS[m], chrom=ds, window_id=wid, read=r, start=s, length=n,
                                         on_strong=ms[s:e].any(), on_any=ma[s:e].any(),
                                         kernel_in=np.nanmean(kin) if np.isfinite(kin).any() else np.nan,
                                         kernel_left=np.nanmean(kl) if np.isfinite(kl).any() else np.nan,
                                         kernel_right=np.nanmean(kr) if np.isfinite(kr).any() else np.nan,
                                         prior_tf=ptf[s:e].mean(),
                                         access=(m6a[hi] - m6a[lo]) / max(called[hi] - called[lo], 1)))
    print(ds, "done", flush=True)
C = pd.DataFrame(rows)
C["kernel_fp"] = (C.kernel_in > 0.8) & (C.kernel_left < 0.2) & (C.kernel_right < 0.2)
C["kernel_known"] = C[["kernel_in", "kernel_left", "kernel_right"]].notna().all(axis=1)
C["short"] = C.length <= 20
C["where"] = np.select([C.on_strong, C.on_any], ["strong motif", "other motif"], "no motif")
C.to_csv("calls.tsv.gz", sep="\t", index=False)
pd.set_option("display.width", 220)
n_reads_kb = {}
for m in MODELS:
    s = pd.read_csv(f"{T}/viterbi_segments_{m}_x1.tsv.gz", sep="\t", usecols=["chrom", "window_id", "read"]).drop_duplicates()
    n_reads_kb[MODELS[m]] = len(s) * 0.512
g = C.groupby(["model", "where"])
S = g.size().rename("calls").to_frame()
S["per_read_kb"] = [S.calls[i] / n_reads_kb[i[0]] for i in S.index]
S["share_short_le20"] = g.short.mean(); S["median_len"] = g.length.median()
S["kernel_fp_share"] = C[C.kernel_known].groupby(["model", "where"]).kernel_fp.mean()
S["kernel_in_median"] = g.kernel_in.median(); S["prior_tf_median"] = g.prior_tf.median(); S["access_median"] = g.access.median()
print(S.round(3).to_string())
print("\nby length, off strong motif:")
C["len_bin"] = pd.cut(C.length, [14, 20, 30, 45, 59, 60])
off = C[~C.on_strong]
t = off.groupby(["model", "len_bin"], observed=True).agg(calls=("length", "size"), kernel_fp=("kernel_fp", "mean"), on_other=("on_any", "mean"), prior_tf=("prior_tf", "median"))
t["per_read_kb"] = [t.calls[i] / n_reads_kb[i[0]] for i in t.index]
print(t.round(3).to_string())
