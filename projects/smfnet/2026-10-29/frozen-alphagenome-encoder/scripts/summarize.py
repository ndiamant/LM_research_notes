import json, sys
from pathlib import Path
import numpy as np

D = Path(sys.argv[1])
PREFIX = sys.argv[2] if len(sys.argv) > 2 else ""
models = ["baseline", "alphagenome", "atac"]
for kind in ["ac", "bg"]:
    data = {m: json.load(open(D / f"{PREFIX}{m}_{kind}.json")) for m in models if (D / f"{PREFIX}{m}_{kind}.json").exists()}
    if not data:
        continue
    any_ = next(iter(data.values()))
    rates, sites = any_["mask_rates"], any_["sites"]
    print(f"\n=== {kind} regions (n={any_['n_regions']}): cross-entropy per hidden call ===")
    for label, site_idx in [("all sites", slice(None)), ("CpG only", 0), ("GpC only", 1)]:
        print(f"-- {label}")
        print("mask   " + "  ".join(f"{m:>10}" for m in data) + "  " + "  ".join(f"{m+' -base':>14}" for m in data if m != "baseline") + "  " + "  ".join(f"{m+' rolled-own':>17}" for m in data))
        for ri, rate in enumerate(rates):
            own, rolled = {}, {}
            for m, d in data.items():
                s = np.array(d["sums"])  # cond, rate, site, 2
                sel = s[:, ri, site_idx]
                if sel.ndim == 3:
                    sel = sel.sum(1)
                own[m] = sel[0, 0] / sel[0, 1]
                rolled[m] = sel[1, 0] / sel[1, 1]
            line = f"{rate:>4} " + "  ".join(f"{own[m]:>10.4f}" for m in data)
            line += "  " + "  ".join(f"{own[m]-own['baseline']:>+14.4f}" for m in data if m != "baseline")
            line += "  " + "  ".join(f"{rolled[m]-own[m]:>+17.4f}" for m in data)
            print(line)
    s = np.array(any_["sums"])
    print("hidden-call share by site (rate 0.95):", dict(zip(sites, np.round(s[0, -1, :, 1] / s[0, -1, :, 1].sum(), 3))))
