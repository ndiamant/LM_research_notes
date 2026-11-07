"""Val loss (EMA weights, deterministic val masks) by epoch: the attention
U-Net AdamW + EMA run vs the patched variants; plus mean train epoch time."""
import csv, os, sys
M = "/scratch/users/diamant/smf_models"
runs = {"base (1 bp)": f"{M}/ONT_mESC_alphagenome_unet_attn_adamw_ema"}
for v in ["p4_narrow", "p4_mid", "p4_wide", "p2", "p4_wide_h2"]:
    runs[v] = f"{M}/ONT_mESC_unet_attn_adamw_ema_{v}"
curves, times = {}, {}
for name, d in runs.items():
    f = f"{d}/csv/metrics.csv"
    if not os.path.exists(f):
        continue
    val, t = {}, []
    for r in csv.DictReader(open(f)):
        if r.get("val/loss"):
            val[int(r["epoch"])] = float(r["val/loss"])
        if r.get("train/epoch_time"):
            t.append(float(r["train/epoch_time"]))
    curves[name], times[name] = val, t
epochs = [0, 1, 2, 3, 4, 5, 7, 9, 12, 15, 19, 24, 29]
print(f"{'epoch':>6s} " + " ".join(f"{n:>12s}" for n in curves))
for e in epochs:
    print(f"{e:6d} " + " ".join(f"{c[e]:12.4f}" if e in c else f"{'':12s}" for c in curves.values()))
print(f"{'s/ep':>6s} " + " ".join(f"{(sum(t[1:]) / max(len(t) - 1, 1)):12.0f}" if len(t) > 1 else f"{'':12s}" for t in times.values()))
