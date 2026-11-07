"""Val cross-entropy by masking level at one epoch, baseline vs variants."""
import csv, os, sys
M = "/scratch/users/diamant/smf_models"
epoch = int(sys.argv[1])
runs = {"base (1 bp)": f"{M}/ONT_mESC_alphagenome_unet_attn_adamw_ema"}
for v in ["p4_narrow", "p4_mid", "p4_wide", "p2", "p4_wide_h2"]:
    runs[v] = f"{M}/ONT_mESC_unet_attn_adamw_ema_{v}"
bins = ["0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]
print(f"epoch {epoch}: val xent by mask rate\n{'run':>12s} " + " ".join(f"{b:>8s}" for b in bins))
for name, d in runs.items():
    f = f"{d}/csv/metrics.csv"
    if not os.path.exists(f):
        continue
    for r in csv.DictReader(open(f)):
        if r.get("val/loss") and int(r["epoch"]) == epoch:
            print(f"{name:>12s} " + " ".join(f"{float(r[f'val_mask/xent_m{b}']):8.4f}" for b in bins))
