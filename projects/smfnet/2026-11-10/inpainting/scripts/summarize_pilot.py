"""Summarize the span-inpainting pilot: per model and span length, pooled
nats per hidden call (marginal, joint, reference), the within-span dependence
used (marginal - joint), and paired per-region differences vs the k=5
attention U-Net (Muon+EMA) with standard errors."""
import json, os
import numpy as np
S = "/scratch/users/diamant/smfnet_inpaint"
models = ["original", "ag_frozen", "atac", "k5", "k9", "k9_1down", "p4_muon_ema"]
data = {m: json.load(open(f"{S}/pilot_{m}.json")) for m in models if os.path.exists(f"{S}/pilot_{m}.json")}
ref = data["k5"]
for s in ref["span_lengths"]:
    print(f"\nspan {s} bp (1,000 val active regions; nats per hidden call)")
    print(f"{'model':>12s} {'marginal':>9s} {'joint':>9s} {'dep. used':>9s} | {'Δ joint vs k5 (SE)':>22s} {'Δ marginal vs k5 (SE)':>22s}")
    r = np.array(ref["per_region"][str(s)])
    for m, d in data.items():
        a = np.array(d["per_region"][str(s)])
        assert np.allclose(a[:, 3], r[:, 3]), f"{m}: different hidden calls"
        tot = a.sum(0)
        n = r[:, 3]; ok = n > 0
        pair = lambda c: ((a[ok, c] - r[ok, c]) / n[ok])
        dj, dm = pair(1), pair(0)
        print(f"{m:>12s} {tot[0] / tot[3]:9.4f} {tot[1] / tot[3]:9.4f} {(tot[0] - tot[1]) / tot[3]:9.4f} | "
              f"{dj.mean():+11.5f} ({dj.std(ddof=1) / np.sqrt(ok.sum()):.5f}) {dm.mean():+11.5f} ({dm.std(ddof=1) / np.sqrt(ok.sum()):.5f})")
    t = r.sum(0)
    print(f"{'reference':>12s} {t[2] / t[3]:9.4f}  (region mean of the other reads, no flanks)")
