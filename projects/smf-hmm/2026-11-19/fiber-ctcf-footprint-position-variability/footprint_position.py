"""How much does CTCF's footprint position vary between reads, beyond what call sparsity explains? (2026-11-19)

Usage: footprint_position.py [CHROMS=chr20] [MAX_SITES=0 (all)] [N_SIM=20] [BM_JITTER=0] [CLEAN=0] [TAG=]
  BM_JITTER: SD (logit units) of a per-simulation perturbation of the read's b_m used to generate the null calls
             (decoding still uses the recorded b_m), standing in for per-read efficiency miscalibration.
  CLEAN=1: keep only footprints whose 10 bp beyond each edge all have kernel < 0.2 (no adjacent protection),
           in both the data and the null. TAG suffixes the output files.
Sites: IDR-bound CTCF motifs (analysis_idr_r6) in windows with >= 20 reads, one per locus (the strongest motif
match; matches within 20 bp of a stronger one are dropped).
Per read, the debiased kernel P(bound) (evaluation.tf_chip.kernel_pbound, +-7 bp, eps 0.02, read b_m) defines the
footprint: the run of kernel > 0.5 through the motif position with the highest kernel, bounded on both sides by
kernel < 0.5 (evidence of access; not NaN), 15-80 bp long. Its centre and edges are recorded in motif coordinates
(motif strand; 0 = motif centre, + = 3').
Null: per site, a consensus footprint (median edges over footprint reads); for each footprint read, N_SIM
simulations of its calls at the same A/T positions within +-100 bp of the motif (protected inside the consensus:
m6A prob eps; accessible outside: m6A prob e_t = sigmoid(beta + f + b_m)); the same detection is applied.
Writes reads.tsv.gz (observed footprints), null.tsv.gz (simulated), sites.tsv.gz (per-site SDs), and prints summaries.
"""
import os, sys
import h5py, numpy as np, pandas as pd
sys.path.insert(0, "/home/users/diamant/repos/smf_hmm")
from evaluation.tf_chip.kernel_pbound import call_llr, kernel_pbound, kmer_codes, load_reference
from evaluation.tf_chip.plot_model_comparison import DEFAULT_KERNEL_REFERENCE, read_efficiencies

OUT = os.path.dirname(os.path.abspath(__file__))
V = "/scratch/users/diamant/gm12878_validation/covariation_windows_v2"
CHROMS = os.environ.get("CHROMS", "chr20").split(",")
MAX_SITES = int(os.environ.get("MAX_SITES", "0"))
N_SIM = int(os.environ.get("N_SIM", "20"))
BM_JITTER = float(os.environ.get("BM_JITTER", "0"))
CLEAN = os.environ.get("CLEAN") == "1"
TAG = os.environ.get("TAG", "")
EPS, H, REGION = 0.02, 7, 100
BETA, F = load_reference(DEFAULT_KERNEL_REFERENCE)
EFF = h5py.File("/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/active_efficiency.h5")
rng = np.random.default_rng(0)


def footprint(k, m0, m1):
    """(left, right) of the motif-overlapping kernel run, or None."""
    seg = k[m0:m1]
    if not np.isfinite(seg).any() or np.nanmax(seg) <= 0.5:
        return None
    p = m0 + int(np.nanargmax(seg))
    lo, hi = p, p
    while lo - 1 >= 0 and np.isfinite(k[lo - 1]) and k[lo - 1] > 0.5:
        lo -= 1
    while hi + 1 < len(k) and np.isfinite(k[hi + 1]) and k[hi + 1] > 0.5:
        hi += 1
    if lo - 1 < 0 or hi + 1 >= len(k) or not (np.isfinite(k[lo - 1]) and np.isfinite(k[hi + 1])):
        return None  # an edge without evidence of access
    if not 15 <= hi - lo + 1 <= 80:
        return None
    if CLEAN:
        flank = np.r_[k[max(lo - 10, 0):lo], k[hi + 1:hi + 11]]
        if len(flank) < 20 or not (np.isfinite(flank).all() and (flank < 0.2).all()):
            return None
    return lo, hi


def to_motif(lo, hi, c, strand):
    """Edges and centre in motif coordinates (motif strand)."""
    if strand == "+":
        return lo - c, hi - c, (lo + hi) / 2 - c
    return c - hi, c - lo, c - (lo + hi) / 2


obs_rows, null_rows = [], []
n_sites = 0
for ds in CHROMS:
    D = f"{V}/{ds}-512-valid90-reads10-flank128"
    S = pd.read_csv(f"{D}/analysis_idr_r6/unique_motif_summaries.tsv.gz", sep="\t")
    S = S[(S.tf == "ctcf") & (S.idr_label == "bound")].sort_values("motif_score", ascending=False)
    # one motif per locus: drop matches whose centre is within 20 bp of a stronger kept match (overlapping or
    # opposite-strand duplicates of the same site would otherwise count the same reads twice)
    keep, kept = [], {}
    for s_ in S.itertuples():
        c_ = (s_.start + s_.end) / 2
        if any(abs(c_ - k_) <= 20 for k_ in kept.get(s_.chrom, [])):
            continue
        kept.setdefault(s_.chrom, []).append(c_); keep.append(s_.Index)
    S = S.loc[keep].sort_values("chip_pooled_sum", ascending=False)
    with h5py.File(f"{D}/dataset.h5") as f:
        for s in S.itertuples():
            if MAX_SITES and n_sites >= MAX_SITES:
                break
            w = f["windows"][s.window_id]; obs = w["observations"][:]
            if len(obs) < 20:
                continue
            dna = w["dna"][:].astype(np.int8); L = len(dna)
            x = np.full((len(obs), L), -1, np.int8); x[:, w["positions"][:]] = obs
            b = read_efficiencies(EFF, w["molecule_id"].asstr()[:])
            k = kernel_pbound(call_llr(x, dna, b, BETA, F, EPS), H)
            m0, m1 = int(s.start - s.window_start), int(s.end - s.window_start)
            c = (m0 + m1 - 1) / 2
            fps = {r: fp for r in range(len(obs)) if (fp := footprint(k[r], m0, m1))}
            n_sites += 1
            for r, (lo, hi) in fps.items():
                le, re, ce = to_motif(lo, hi, c, s.strand)
                obs_rows.append(dict(site=s.site_key, chip=s.chip_pooled_sum, read=r, b_m=b[r], left=le, right=re, centre=ce,
                                     length=hi - lo + 1, n_reads=len(obs)))
            if len(fps) < 5:
                continue
            los = np.array([v[0] for v in fps.values()]); his = np.array([v[1] for v in fps.values()])
            clo, chi = int(np.median(los)), int(np.median(his))
            codes = kmer_codes(dna)
            r0, r1 = max(int(c) - REGION, 0), min(int(c) + REGION + 1, L)
            for r in fps:
                called = np.flatnonzero((x[r, r0:r1] >= 0) & (codes[r0:r1] >= 0)) + r0
                inside = (called >= clo) & (called <= chi)
                for sim in range(N_SIM):
                    z = BETA + F[codes[called]] + b[r] + (rng.normal(0, BM_JITTER) if BM_JITTER else 0.0)
                    p_m6a = np.where(inside, EPS, 1 / (1 + np.exp(-z)))
                    xs = x[r:r + 1].copy()
                    xs[0, called] = np.where(rng.random(len(called)) < p_m6a, 0, 1)
                    ks = kernel_pbound(call_llr(xs, dna, b[r:r + 1], BETA, F, EPS), H)[0]
                    fp = footprint(ks, m0, m1)
                    if fp:
                        le, re, ce = to_motif(*fp, c, s.strand)
                        null_rows.append(dict(site=s.site_key, read=r, sim=sim, left=le, right=re, centre=ce, length=fp[1] - fp[0] + 1))
    print(ds, "sites so far", n_sites, flush=True)

O, N = pd.DataFrame(obs_rows), pd.DataFrame(null_rows)
O.to_csv(f"{OUT}/reads{TAG}.tsv.gz", sep="\t", index=False); N.to_csv(f"{OUT}/null{TAG}.tsv.gz", sep="\t", index=False)
site = O.groupby("site").agg(n_fp=("centre", "size"), sd_obs=("centre", "std"), sd_left_obs=("left", "std"),
                             sd_right_obs=("right", "std"), med_centre=("centre", "median"), med_len=("length", "median"))
site = site.join(N.groupby("site").agg(sd_null=("centre", "std"), sd_left_null=("left", "std"), sd_right_null=("right", "std"),
                                       null_detect=("centre", "size")))
mad = lambda v: float(np.median(np.abs(v - np.median(v))))
site = site.join(O.groupby("site").centre.agg(mad_obs=mad)).join(N.groupby("site").centre.agg(mad_null=mad))
site.to_csv(f"{OUT}/sites{TAG}.tsv.gz", sep="\t")
t = site.dropna(subset=["sd_null"])
pd.set_option("display.width", 200)
print(f"sites scanned {n_sites}; reads with a footprint {len(O)}; sites with >= 5 footprint reads {len(t)}")
print("median footprint length", O.length.median(), "| median centre (motif coords)", O.centre.median(),
      "| median left / right edge", O.left.median(), O.right.median())
print("per-site SD of centre: median observed", round(t.sd_obs.median(), 2), "median null", round(t.sd_null.median(), 2),
      "| median ratio", round((t.sd_obs / t.sd_null).median(), 2), "| sites with obs > null:", round((t.sd_obs > t.sd_null).mean(), 3))
print("per-site MAD of centre: median observed", round(t.mad_obs.median(), 2), "null", round(t.mad_null.median(), 2),
      "| median ratio", round((t.mad_obs / t.mad_null.replace(0, np.nan)).median(), 2))
print("left-edge SD: obs", round(t.sd_left_obs.median(), 2), "null", round(t.sd_left_null.median(), 2),
      "| right-edge SD: obs", round(t.sd_right_obs.median(), 2), "null", round(t.sd_right_null.median(), 2))
dev_obs = O.groupby("site").centre.transform(lambda v: v - v.median())
dev_null = N.groupby("site").centre.transform(lambda v: v - v.median())
print("pooled |centre - site median|: obs quantiles", np.quantile(np.abs(dev_obs[O.site.isin(t.index)]), [.5, .9, .99]).round(1),
      "null", np.quantile(np.abs(dev_null), [.5, .9, .99]).round(1))
