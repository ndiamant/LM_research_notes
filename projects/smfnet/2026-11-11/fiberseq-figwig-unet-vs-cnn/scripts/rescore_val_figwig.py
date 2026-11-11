"""Fixed-mask val loss of several checkpoints on the same figwig val batches.

The masks are seeded by batch index (SMFNet.validation_step), and the batches are the
figwig val split in its fixed order, so every checkpoint sees the same reads, crops and
masks. Prints and writes one JSON line per checkpoint.
"""

import argparse
import json
import os

import torch
from pytorch_lightning import Callback, Trainer


class BatchLosses(Callback):
    """Collects each val batch's loss, for paired comparisons between checkpoints."""

    def __init__(self):
        self.losses = []

    def on_validation_batch_end(self, trainer, pl_module, outputs, batch, batch_idx, *args):
        self.losses.append(float(outputs))

from smf_net.figwig_data import FigwigDataModuleConfig, FigwigSMFDataModule
from smf_net.model import load_model

ap = argparse.ArgumentParser()
ap.add_argument("--checkpoint", action="append", required=True, help="name=path, repeatable")
ap.add_argument("--out", required=True)
ap.add_argument("--num-workers", type=int, default=3)
args = ap.parse_args()

F = os.environ["SCRATCH"] + "/fiberseq/GM12878"
cfg = FigwigDataModuleConfig(
    fasta_path="/scratch/groups/btrippe/brian/projects/SMF/fiberseq/reference/hg38.fa",
    figwig_dir=F + "/figwig",
    regions_bed=F + "/kept_windows.bed",
    assay_type="fiberseq_m6a",
    region_length=2114,
    model_input_length=2048,
    min_overlap=500,
    min_spanning_molecules=32,
    val_chrm=["chr8", "chr20"],
    num_workers=args.num_workers,
)
data = FigwigSMFDataModule(cfg)
data.setup("validate")
print(f"{len(data.val_dataset)} val regions")

with open(args.out, "a") as out:
    for spec in args.checkpoint:
        name, path = spec.split("=", 1)
        model = load_model(path)
        if model.cfg.input_length != cfg.model_input_length or model.cfg.assay_type != cfg.assay_type:
            raise ValueError(f"{name}: model does not match the data")
        batches = BatchLosses()
        trainer = Trainer(
            accelerator="gpu", devices=1, precision="bf16-mixed", logger=False,
            enable_progress_bar=False, enable_model_summary=False, callbacks=[batches],
        )
        torch.backends.cudnn.benchmark = True
        (scores,) = trainer.validate(model, data.val_dataloader(), verbose=False)
        record = {"name": name, "checkpoint": path, "val_loss": scores["val/loss"],
                  "n_regions": len(data.val_dataset), "batch_losses": batches.losses}
        print(json.dumps(record))
        out.write(json.dumps(record) + "\n")
        del model, trainer
        torch.cuda.empty_cache()
