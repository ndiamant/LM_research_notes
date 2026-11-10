"""50-epoch inpainting variants of the k=9 attention U-Net (Muon+EMA): the
epoch-49 checkpoint's fixed-mask val loss and span-inpainting scores, vs
the original k=9 run's epoch-49 checkpoint (k9_full_run) and the 50-epoch
rerun (control), with paired per-region differences vs the control."""
import json, os
import numpy as np
S = "/scratch/users/diamant/smfnet_inpaint"
names = ["k9_full_run", "control", "linear", "dilated", "local_attn", "bottleneck16", "read_masked", "window"]
val = {n: json.load(open(f"{S}/rescored_e49_{n}.json"))["val_loss_fixed_masks"]
       for n in names if os.path.exists(f"{S}/rescored_e49_{n}.json")}
pil = {n: json.load(open(f"{S}/pilot_e49_{n}.json")) for n in names if os.path.exists(f"{S}/pilot_e49_{n}.json")}
ref = pil.get("control")
# k9_full_run: the original k=9 run (same config as control), scored on an H100
print(f"{'variant':>13s} {'val loss':>9s} | " + " | ".join(f"{s} bp one-shot / joint (Δ joint vs control, SE)" for s in (50, 150)))
for n in names:
    if n not in pil:
        print(f"{n:>13s} {val.get(n, float('nan')):9.5f} | (pending)"); continue
    cells = []
    for s in (50, 150):
        a = np.array(pil[n]["per_region"][str(s)]); t = a.sum(0)
        cell = f"{t[0] / t[3]:.4f} / {t[1] / t[3]:.4f}"
        r = None if ref is None else np.array(ref["per_region"][str(s)])
        if r is not None and np.allclose(a[:, 3], r[:, 3]):  # same spans (the GPU type sets the CUDA draws)
            ok = r[:, 3] > 0
            d = (a[ok, 1] - r[ok, 1]) / r[ok, 3]
            cell += f" ({d.mean():+.5f}, {d.std(ddof=1) / np.sqrt(ok.sum()):.5f})"
        elif r is not None:
            cell += " (other spans: L40S)"
        cells.append(cell)
    print(f"{n:>13s} {val.get(n, float('nan')):9.5f} | " + " | ".join(cells))
