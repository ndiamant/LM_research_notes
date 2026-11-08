"""Summarize memorization_gap.py results: per model, cross-entropy per hidden
call (all site types pooled) by mask rate on train and val regions, the
val - train gap, the DNA benefit (rolled - own) per split, and the loss
increase from mutating 2% of bases, per split. Active and background regions
are pooled (summed over both h5s)."""
import json, os, sys
import numpy as np
O = "/scratch/users/diamant/smfnet_memorization/results"
models = sys.argv[1:] or ["baseline", "ag_frozen", "ag_ftall", "attn_muon_ema", "atac"]
def load(name, split, kinds=("ac", "bg")):
    tot = None
    for kind in kinds:
        f = f"{O}/{name}_{kind}_{split}.json"
        if not os.path.exists(f):
            return None, None
        d = json.load(open(f))
        a = np.array(d["sums" if "sums" in d else "ceiling"])
        tot = a if tot is None else tot + a
    return tot, d
ceil = {}
for split in ("train", "val"):
    c, _ = load("ceiling", split)
    if c is not None:
        ceil[split] = c[..., 0].sum() / c[..., 1].sum()
if ceil:
    print("region-mean ceiling (all valid calls, leave-one-out over 100 reads): "
          + ", ".join(f"{s} {v:.4f}" for s, v in ceil.items()))
for name in models:
    tr, d = load(name, "train")
    va, _ = load(name, "val")
    if tr is None or va is None:
        print(f"\n{name}: missing"); continue
    rates = d["mask_rates"]; conds = d["conditions"]
    xent = lambda a: a[..., 0].sum(-1) / a[..., 1].sum(-1)  # [condition, rate]
    t, v = xent(tr), xent(va)
    own, rolled, mut = (conds.index(c) for c in ("own", "rolled", "mutated"))
    print(f"\n{name}")
    print(f"{'mask':>5s} {'train':>7s} {'val':>7s} {'gap':>7s} | {'DNA benefit tr':>14s} {'val':>7s} | {'mutation cost tr':>16s} {'val':>7s}")
    for i, r in enumerate(rates):
        print(f"{r:5.2f} {t[own, i]:7.4f} {v[own, i]:7.4f} {v[own, i] - t[own, i]:+7.4f} | "
              f"{t[rolled, i] - t[own, i]:14.4f} {v[rolled, i] - v[own, i]:7.4f} | "
              f"{t[mut, i] - t[own, i]:16.4f} {v[mut, i] - v[own, i]:7.4f}")
