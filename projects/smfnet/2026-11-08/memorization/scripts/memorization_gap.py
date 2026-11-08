"""Train vs val cross-entropy per hidden call, by mask rate and site type, to
test whether models memorize their training regions.

For one h5 and one split (train or val), a fixed subset of regions is scored
the same way on both splits: deterministic read choice per region, the
middle of the region (no jitter), no reverse complement, and the same
seeded random masks (a uniform draw per call, seeded by batch index, hidden
where draw < mask rate; all valid calls at mask rate 1). Each mask rate is
scored with three DNA conditions:

- own: the region's DNA;
- rolled: the DNA of another region in the batch (how much the model uses
  its own region's DNA);
- mutated: the region's DNA with a fraction of bases substituted at random
  (how brittle its use of the exact sequence is).

Separately, the region-mean ceiling: each call predicted by the mean of the
same site over the region's other reads (100 reads per region,
leave-one-out, Jeffreys prior), the best a DNA-only model could do on that
region up to estimation noise.

Works with the nd_atac code too (PYTHONPATH pointing at its worktree).
"""

import argparse
import json

import numpy as np
import torch
import torch.nn.functional as F

from smf_net.data import SparseDataModuleConfig, SparseSMFDataModule, compute_masks_cg_gc
from smf_net.diffusion_schedule import MASK_TOKEN
from smf_net.model import load_model

MASK_RATES = [round(m, 2) for m in np.arange(0.05, 1.0, 0.1)] + [1.0]
SITES = ["cg_only", "gc_only", "gcg_both"]
CONDITIONS = ["own", "rolled", "mutated"]

parser = argparse.ArgumentParser()
parser.add_argument("--checkpoint", default=None, help="omit to compute only the region-mean ceiling")
parser.add_argument("--h5", required=True)
parser.add_argument("--split", choices=["train", "val"], required=True)
parser.add_argument("--max-regions", type=int, default=4000)
parser.add_argument("--fasta", default="/scratch/groups/btrippe/ndiamant/mm10.fa")
parser.add_argument("--atac-bigwig", default=None)
parser.add_argument("--assay-type", default="cg_gc_accessibility")
parser.add_argument("--input-length", type=int, default=2048)
parser.add_argument("--mutation-rate", type=float, default=0.02)
parser.add_argument("--batch-size", type=int, default=128)
parser.add_argument("--out", required=True)
args = parser.parse_args()

model = None if args.checkpoint is None else load_model(args.checkpoint, map_location="cuda").eval()


def datamodule(reads_per_region: int) -> SparseSMFDataModule:
    cfg = dict(
        fasta_path=args.fasta,
        h5_path=[args.h5],
        assay_type=model.cfg.assay_type if model else args.assay_type,
        batch_size=args.batch_size,
        num_workers=4,
        reads_per_region=reads_per_region,
        model_input_length=model.cfg.input_length if model else args.input_length,
        max_regions_per_h5_per_split=args.max_regions,
    )
    if args.atac_bigwig:
        cfg.update(
            atac_bigwig_path=args.atac_bigwig,
            atac_target_library_size=model.cfg.atac_target_library_size,
        )
    dm = SparseSMFDataModule(SparseDataModuleConfig(**cfg))
    # Score the split as validation is scored (train=False: no jitter or
    # reverse complement, reads chosen by region index).
    dm._setup_split(args.split, train=False)
    return dm


def loader(dm: SparseSMFDataModule):
    dm.val_dataset = getattr(dm, f"{args.split}_dataset")
    return dm.val_dataloader()


def middle(dna, smf):
    L = dna.shape[-1]
    x0 = smf.cuda().long()[..., L // 4 : -L // 4]  # B x R x L/2
    dna_int = torch.where(dna.sum(1) > 0, dna.argmax(1), -1)
    cg, gc = (m[:, None, L // 4 : -L // 4] for m in compute_masks_cg_gc(dna_int))
    return x0, [cg & ~gc, gc & ~cg, cg & gc]


def mutate(dna: torch.Tensor, rate: float, generator: torch.Generator) -> torch.Tensor:
    """Substitute a random different base at a fraction of the known positions."""
    known = dna.sum(1) > 0  # B x L
    hit = known & (torch.rand(known.shape, generator=generator, device=dna.device) < rate)
    base = dna.argmax(1)
    shift = torch.randint(1, 4, base.shape, generator=generator, device=dna.device)
    new = torch.where(hit, (base + shift) % 4, base)
    out = F.one_hot(new, 4).transpose(1, 2).to(dna.dtype)
    return torch.where(known[:, None], out, dna)


result = {
    "checkpoint": args.checkpoint,
    "h5": args.h5,
    "split": args.split,
    "max_regions": args.max_regions,
    "mutation_rate": args.mutation_rate,
    "mask_rates": MASK_RATES,
    "sites": SITES,
    "conditions": CONDITIONS,
}

with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16, enabled=model is not None):
    if model is not None:
        # [condition, mask rate, site, (summed cross-entropy, hidden calls)]
        sums = torch.zeros(len(CONDITIONS), len(MASK_RATES), len(SITES), 2, dtype=torch.float64)
        n_regions = 0
        for batch_idx, batch in enumerate(loader(datamodule(32))):
            dna = batch.dna.cuda()
            atac = getattr(batch, "atac", None)
            atac = None if atac is None else atac.cuda()
            x0, site_masks = middle(dna, batch.smf)
            valid = x0 != -1
            generator = torch.Generator(device="cuda").manual_seed(batch_idx)
            draw = torch.rand(x0.shape, generator=generator, device="cuda")
            mutated = mutate(dna, args.mutation_rate, generator)
            roll = torch.roll(torch.arange(dna.shape[0], device="cuda"), 1)
            n_regions += dna.shape[0]
            for ci, condition in enumerate(CONDITIONS):
                cond_dna = {"own": dna, "rolled": dna[roll], "mutated": mutated}[condition]
                kwargs = {}
                if atac is not None:
                    kwargs["atac"] = atac[roll] if condition == "rolled" else atac
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
        result.update(n_regions=n_regions, sums=sums.tolist())
    else:
        # Region-mean ceiling: [site, (summed cross-entropy, calls)]
        ceiling = torch.zeros(len(SITES), 2, dtype=torch.float64)
        n_regions = 0
        for batch in loader(datamodule(100)):
            dna = batch.dna.cuda()
            x0, site_masks = middle(dna, batch.smf)
            valid = x0 != -1
            ones = (x0 == 1).sum(1, keepdim=True)  # B x 1 x L/2, over reads
            count = valid.sum(1, keepdim=True)
            # Leave-one-out mean of the other reads, Jeffreys prior.
            others_ones = ones - (x0 == 1).long()
            others = count - valid.long()
            p = ((others_ones + 0.5) / (others + 1.0)).float()
            bce = F.binary_cross_entropy(p, x0.clamp(min=0).float(), reduction="none")
            n_regions += dna.shape[0]
            for si, site in enumerate(site_masks):
                selected = valid & site
                ceiling[si, 0] += bce[selected].double().sum().cpu()
                ceiling[si, 1] += selected.sum().cpu()
        result.update(n_regions=n_regions, ceiling=ceiling.tolist())

json.dump(result, open(args.out, "w"), indent=1)
print("wrote", args.out, "regions", result["n_regions"])
