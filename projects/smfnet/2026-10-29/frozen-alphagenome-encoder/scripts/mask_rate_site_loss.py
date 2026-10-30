"""Val cross-entropy per hidden call by mask rate and site type (CpG/GpC).

Every model sees the same val batches (deterministic read choice per region)
and the same random masks (a uniform draw per call, seeded by batch index,
hidden where draw < mask rate). Each mask rate is evaluated with the region's
own conditioning and with the conditioning (DNA, and ATAC if used) rolled by
one region within the batch.

Works with both the nd_alpha_genome and nd_atac code (the latter via
PYTHONPATH pointing at its worktree).
"""

import argparse
import json

import numpy as np
import torch
import torch.nn.functional as F

from smf_net.data import SparseDataModuleConfig, SparseSMFDataModule, compute_masks_cg_gc
from smf_net.diffusion_schedule import MASK_TOKEN
from smf_net.model import load_model

MASK_RATES = [round(m, 2) for m in np.arange(0.05, 1.0, 0.1)]
SITES = ["cg_only", "gc_only", "gcg_both"]
CONDITIONS = ["own", "rolled"]

parser = argparse.ArgumentParser()
parser.add_argument("--checkpoint", required=True)
parser.add_argument("--h5", required=True)
parser.add_argument("--fasta", default="/scratch/groups/btrippe/ndiamant/mm10.fa")
parser.add_argument("--atac-bigwig", default=None)
parser.add_argument("--batch-size", type=int, default=128)
parser.add_argument("--out", required=True)
parser.add_argument("--mask-rates", type=float, nargs="+", default=MASK_RATES)
args = parser.parse_args()
MASK_RATES = [float(m) for m in args.mask_rates]

model = load_model(args.checkpoint, map_location="cuda").eval()
cfg = dict(
    fasta_path=args.fasta,
    h5_path=[args.h5],
    assay_type=model.cfg.assay_type,
    batch_size=args.batch_size,
    num_workers=4,
    reads_per_region=32,
    model_input_length=model.cfg.input_length,
)
if args.atac_bigwig:
    cfg.update(
        atac_bigwig_path=args.atac_bigwig,
        atac_target_library_size=model.cfg.atac_target_library_size,
    )
datamodule = SparseSMFDataModule(SparseDataModuleConfig(**cfg))
datamodule.setup("validate")

# [condition, mask rate, site, (summed cross-entropy, hidden calls)]
sums = torch.zeros(len(CONDITIONS), len(MASK_RATES), len(SITES), 2, dtype=torch.float64)
n_regions = 0
with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
    for batch_idx, batch in enumerate(datamodule.val_dataloader()):
        dna = batch.dna.cuda()
        atac = getattr(batch, "atac", None)
        atac = None if atac is None else atac.cuda()
        L = dna.shape[-1]
        x0 = batch.smf.cuda().long()[..., L // 4 : -L // 4]  # B x R x L/2
        valid = x0 != -1
        dna_int = torch.where(dna.sum(1) > 0, dna.argmax(1), -1)
        cg, gc = (m[:, None, L // 4 : -L // 4] for m in compute_masks_cg_gc(dna_int))
        site_masks = [cg & ~gc, gc & ~cg, cg & gc]
        generator = torch.Generator(device="cuda").manual_seed(batch_idx)
        draw = torch.rand(x0.shape, generator=generator, device="cuda")
        roll = torch.roll(torch.arange(dna.shape[0], device="cuda"), 1)
        n_regions += dna.shape[0]
        for ci, condition in enumerate(CONDITIONS):
            cond_dna = dna if condition == "own" else dna[roll]
            kwargs = {}
            if atac is not None:
                kwargs["atac"] = atac if condition == "own" else atac[roll]
            for ri, rate in enumerate(MASK_RATES):
                hidden = valid if rate >= 1.0 else valid & (draw < rate)
                x_t = torch.where(hidden | ~valid, MASK_TOKEN, x0)
                logits = model.model(model._encode_smf_tokens(x_t), cond_dna, **kwargs)
                bce = F.binary_cross_entropy_with_logits(
                    logits.squeeze(2).float(), x0.clamp(min=0).float(), reduction="none"
                )
                for si, site in enumerate(site_masks):
                    selected = hidden & site
                    sums[ci, ri, si, 0] += bce[selected].double().sum().cpu()
                    sums[ci, ri, si, 1] += selected.sum().cpu()

json.dump(
    {
        "checkpoint": args.checkpoint,
        "h5": args.h5,
        "atac_bigwig": args.atac_bigwig,
        "n_regions": n_regions,
        "mask_rates": MASK_RATES,
        "sites": SITES,
        "conditions": CONDITIONS,
        "sums": sums.tolist(),
    },
    open(args.out, "w"),
    indent=1,
)
print("wrote", args.out, "regions", n_regions)
