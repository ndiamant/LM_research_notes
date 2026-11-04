"""Val/loss by epoch for the baseline (from its local .wandb file) and the
runs with CSV logs, to compare early training."""
import glob
import json
import sys

import pandas as pd
from wandb.proto import wandb_internal_pb2
from wandb.sdk.internal import datastore

M = "/scratch/users/diamant/smf_models"


def wandb_val_loss(path):
    ds = datastore.DataStore()
    ds.open_for_scan(path)
    rows = {}
    while True:
        data = ds.scan_data()
        if data is None:
            break
        rec = wandb_internal_pb2.Record()
        rec.ParseFromString(data)
        if rec.WhichOneof("record_type") != "history":
            continue
        items = {i.key or "/".join(i.nested_key): json.loads(i.value_json) for i in rec.history.item}
        if "val/loss" in items and "epoch" in items:
            rows[int(items["epoch"])] = items["val/loss"]
    return pd.Series(rows).sort_index()


def csv_val_loss(run):
    df = pd.read_csv(f"{M}/{run}/csv/metrics.csv")
    return df.dropna(subset=["val/loss"]).groupby("epoch")["val/loss"].first()


curves = {"baseline": wandb_val_loss(glob.glob("/scratch/users/diamant/wandb/wandb/run-20260811_175005-3ekr0isw/*.wandb")[0])}
for name, run in [("ag_frozen", "ONT_mESC_alphagenome_frozen"), ("ag_ftall", "ONT_mESC_alphagenome_ftall_run2"),
                  ("atac", "ONT_mESC_ATAC_TKO_pooled"), ("cheri_muon", "ONT_mESC_cheri_muon"), ("cheri_lr0.01", "ONT_mESC_cheri_muon_lr0.01"), ("ag_unet", "ONT_mESC_alphagenome_frozen_unet")]:
    try:
        curves[name] = csv_val_loss(run)
    except (FileNotFoundError, KeyError) as e:
        print("missing", name, e)
df = pd.DataFrame(curves)
epochs = [int(e) for e in sys.argv[1:]] or [0, 1, 2, 4, 9, 19, 49, 99, 199, 377, 399]
pd.set_option("display.width", 200)
print(df.loc[[e for e in epochs if e in df.index]].round(4).to_string())
print("\ncheri_muon epochs logged:", int(df["cheri_muon"].notna().sum()) if "cheri_muon" in df else 0)
