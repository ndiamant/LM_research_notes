"""Are low-efficiency TF calls real protection or residual context bias? Held-out calls inside MAP TF calls (2026-11-24).

Every 8th val window (chr20, chr8; sorted ids). Each read's observed A/T calls are masked independently at
RATE (seeded per window), the masked reads are Viterbi-decoded (prior x1) by each model, and every masked
call is recorded with the MAP group at its position. Masked calls did not inform the path, so their m6A rate
inside a TF call is an unbiased test of protection there.
Per masked call: truth (m6A), ctx = sigmoid(beta + f(context)) (position efficiency, no b_m),
e = sigmoid(beta + f(context) + b_m) (expected rate if accessible, kernel reference), group; for TF
positions also the call's start / length, whether it overlaps the strong-footprint mask (FPR 1e-6, +-10 bp),
and the call's mean ctx over its A/T sites.
Writes heldout_calls.npz; analysis in analyze.py.
"""
import os, sys
import h5py, numpy as np
sys.path.insert(0, "/home/users/diamant/repos/smf_hmm")
from evaluation.tf_chip.kernel_pbound import kmer_codes, load_reference
from evaluation.tf_chip.plot_model_comparison import DEFAULT_KERNEL_REFERENCE, GROUPS, Decoder, motif_scanners, read_efficiencies
from evaluation.tf_chip.specificity_metrics import DEFAULT_MOTIFS, REFERENCES, motif_mask
from pathlib import Path
RATE = float(os.environ.get("RATE", "0.2"))
V = "/scratch/users/diamant/gm12878_validation/covariation_windows_v2"
S = Path("/scratch/users/diamant")
MODELS = {"AlphaGenome": S / "2026-11-17_fiber_hist190_tf60_ol61_kl300_alphagenome_bf16",
          "sequence-free": S / "2026-11-18_fiber_hist190_tf60_ol61_kl300_constant_prior"}
beta, f = load_reference(DEFAULT_KERNEL_REFERENCE)
strong = motif_scanners(DEFAULT_MOTIFS, 1e-6, set((REFERENCES / "gm12878_strong_footprint_tfs.txt").read_text().split()))
EFF = h5py.File("/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/active_efficiency.h5")
TF = GROUPS.index("TF")
out = {k: [] for k in ("model", "window", "read", "pos", "truth", "ctx", "e", "group", "call_start", "call_len", "call_strong", "call_ctx")}
for mi, (name, run) in enumerate(MODELS.items()):
    dec = Decoder(run, 1.0)
    wi = 0
    for ds in ("chr20", "chr8"):
        with h5py.File(f"{V}/{ds}-512-valid90-reads10-flank128/dataset.h5") as fh:
            for wid in sorted(fh["windows"])[::8]:
                w = fh["windows"][wid]; dna = w["dna"][:].astype(np.int64); L = len(dna)
                calls = np.full((len(w["observations"]), L), -1, np.int8); calls[:, w["positions"][:]] = w["observations"][:]
                b = read_efficiencies(EFF, w["molecule_id"].asstr()[:])
                rng = np.random.default_rng([0, wi])
                masked = (calls >= 0) & (rng.random(calls.shape) < RATE)
                x = np.where(masked, -1, calls).astype(np.int32)
                z = dec(dna, x, b)  # (reads, L) GROUPS index
                codes = kmer_codes(dna); fz = beta + np.where(codes >= 0, f[np.maximum(codes, 0)], 0.0)
                ctx = 1 / (1 + np.exp(-fz)); has = codes >= 0
                smask = motif_mask(dna, strong, 10)
                r_idx, p_idx = np.nonzero(masked)
                g = z[r_idx, p_idx]
                cs = np.full(len(r_idx), -1); cl = np.zeros(len(r_idx), int); cst = np.zeros(len(r_idx), bool); cc = np.full(len(r_idx), np.nan)
                for j in np.flatnonzero(g == TF):
                    r, p = r_idx[j], p_idx[j]
                    a = p
                    while a > 0 and z[r, a - 1] == TF: a -= 1
                    e_ = p
                    while e_ < L - 1 and z[r, e_ + 1] == TF: e_ += 1
                    cs[j], cl[j], cst[j] = a, e_ - a + 1, smask[a:e_ + 1].any()
                    sites = has[a:e_ + 1]
                    cc[j] = ctx[a:e_ + 1][sites].mean() if sites.any() else np.nan
                for k, v in (("model", np.full(len(r_idx), mi)), ("window", np.full(len(r_idx), wi)), ("read", r_idx), ("pos", p_idx),
                             ("truth", calls[r_idx, p_idx] == 0), ("ctx", ctx[p_idx]), ("e", 1 / (1 + np.exp(-(fz[p_idx] + b[r_idx])))),
                             ("group", g), ("call_start", cs), ("call_len", cl), ("call_strong", cst), ("call_ctx", cc)):
                    out[k].append(v)
                wi += 1
                if wi % 100 == 0:
                    print(name, wi, flush=True)
    del dec
np.savez_compressed("heldout_calls.npz", models=np.array(list(MODELS)), **{k: np.concatenate(v) for k, v in out.items()})
print("done", sum(len(v) for v in out["truth"]))
