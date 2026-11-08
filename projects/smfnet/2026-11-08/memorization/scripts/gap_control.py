"""Train - val gap with the region's own DNA vs with another region's DNA
(rolled), per mask rate. A gap that remains with the wrong DNA is not
DNA memorization."""
import json, numpy as np
O = "results"
def x(name, split):
    a = sum(np.array(json.load(open(f"{O}/{name}_{k}_{split}.json"))["sums"]) for k in ("ac", "bg"))
    return a[..., 0].sum(-1) / a[..., 1].sum(-1)
rates = json.load(open(f"{O}/baseline_ac_train.json"))["mask_rates"]
models = ["baseline", "ag_frozen", "ag_ftall", "attn_muon_ema", "p4_muon_ema", "atac"]
print("val - train, own DNA / rolled DNA")
print("mask  " + "  ".join(f"{m:>17s}" for m in models))
for i, r in enumerate(rates):
    row = []
    for m in models:
        t, v = x(m, "train"), x(m, "val")
        row.append(f"{v[0, i] - t[0, i]:+.4f}/{v[1, i] - t[1, i]:+.4f}")
    print(f"{r:4.2f}  " + "  ".join(f"{c:>17s}" for c in row))
