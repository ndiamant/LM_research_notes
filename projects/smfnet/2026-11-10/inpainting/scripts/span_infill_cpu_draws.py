"""Span-inpainting likelihood: hide one contiguous span of each read and score
the model's predictions of the true calls in it, given the rest of the read
and the DNA. Sampling-free.

For each read and span length S, the span start is drawn uniformly (seeded by
batch index) so that the span lies within the read's covered extent; reads
too short for the span are skipped. Hidden calls are the read's valid calls
inside the span; all other calls of the read stay visible. Two scores per
span, in nats per hidden call:

- marginal: every hidden call predicted in one forward pass from the flanks
  alone (as if independent of the rest of the span);
- joint: the exact chain-rule likelihood of the span's true calls, revealing
  them one at a time in a random order (teacher forcing). joint <= marginal
  in expectation; marginal - joint is the within-span dependence the model
  uses.

Also a reference that needs no model: each hidden call predicted by the mean
of the same site over the region's other reads in the batch (leave-one-out,
Jeffreys prior), which ignores the read's own flanks.

Saves per-region sums so paired standard errors can be computed. Works with
the nd_atac code via PYTHONPATH (with --atac-bigwig).
"""

import argparse
import json

import numpy as np
import torch
import torch.nn.functional as F

from smf_net.data import SparseDataModuleConfig, SparseSMFDataModule
from smf_net.diffusion_schedule import MASK_TOKEN
from smf_net.model import load_model

parser = argparse.ArgumentParser()
parser.add_argument("--checkpoint", required=True)
parser.add_argument("--h5", default="/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions.h5")
parser.add_argument("--fasta", default="/scratch/groups/btrippe/ndiamant/mm10.fa")
parser.add_argument("--atac-bigwig", default=None)
parser.add_argument("--max-regions", type=int, default=1000)
parser.add_argument("--span-lengths", type=int, nargs="+", default=[20, 50, 150])
parser.add_argument("--batch-size", type=int, default=32)
parser.add_argument("--num-workers", type=int, default=3)
parser.add_argument("--out", required=True)
args = parser.parse_args()

model = load_model(args.checkpoint, map_location="cuda").eval()
cfg = dict(
    fasta_path=args.fasta,
    h5_path=[args.h5],
    assay_type=model.cfg.assay_type,
    batch_size=args.batch_size,
    num_workers=args.num_workers,
    reads_per_region=32,
    model_input_length=model.cfg.input_length,
    max_regions_per_h5_per_split=args.max_regions,
)
if args.atac_bigwig:
    cfg.update(
        atac_bigwig_path=args.atac_bigwig,
        atac_target_library_size=model.cfg.atac_target_library_size,
    )
datamodule = SparseSMFDataModule(SparseDataModuleConfig(**cfg))
datamodule.setup("validate")


def span_mask(valid: torch.Tensor, length: int, generator: torch.Generator) -> torch.Tensor:
    """(B, R, L) bool: one span of `length` bp per read within its covered extent."""
    L = valid.shape[-1]
    pos = torch.arange(L, device=valid.device)
    has = valid.any(-1)
    first = torch.where(valid, pos, L).amin(-1)
    last = torch.where(valid, pos, -1).amax(-1)
    room = last - first + 1 - length  # number of extra start positions
    ok = has & (room >= 0)
    u = torch.rand(valid.shape[:-1], generator=generator).to(valid.device)
    start = first + (u * (room.clamp(min=0) + 1)).floor().long()
    span = (pos >= start[..., None]) & (pos < start[..., None] + length)
    return span & ok[..., None] & valid


def logits_of(x_t, dna, atac):
    kwargs = {} if atac is None else {"atac": atac}
    return model.model(model._encode_smf_tokens(x_t), dna, **kwargs).squeeze(2).float()


region_ids = []
# per span length: per-region [marginal sum, joint sum, reference sum, hidden calls]
per_region = {s: [] for s in args.span_lengths}
with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
    for batch_idx, batch in enumerate(datamodule.val_dataloader()):
        dna = batch.dna.cuda()
        atac = getattr(batch, "atac", None)
        atac = None if atac is None else atac.cuda()
        L = dna.shape[-1]
        x0 = batch.smf.cuda().long()[..., L // 4 : -L // 4]  # B x R x L/2
        valid = x0 != -1
        target = x0.clamp(min=0).float()
        # Leave-one-out mean over the region's other reads, Jeffreys prior.
        ones = (x0 == 1).sum(1, keepdim=True)
        count = valid.sum(1, keepdim=True)
        p_ref = ((ones - (x0 == 1).long() + 0.5) / (count - valid.long() + 1.0)).float()
        ref_nll = -(target * p_ref.log() + (1 - target) * (1 - p_ref).log())
        # CPU draws, so the spans are the same on every GPU type.
        generator = torch.Generator().manual_seed(batch_idx)
        for s in args.span_lengths:
            hidden = span_mask(valid, s, generator)
            base = torch.where(hidden | ~valid, MASK_TOKEN, x0)
            nll = F.binary_cross_entropy_with_logits(logits_of(base, dna, atac), target, reduction="none")
            marginal = (nll * hidden).sum((1, 2))
            # Joint: reveal hidden calls one at a time in a random order.
            keys = torch.rand(hidden.shape, generator=generator).to(hidden.device)
            keys = torch.where(hidden, keys, 2.0)
            rank = keys.argsort(-1).argsort(-1)  # 0 = first revealed
            n_hidden = hidden.sum(-1)
            joint = torch.zeros_like(marginal)
            x_t = base.clone()
            for j in range(int(n_hidden.max())):
                current = hidden & (rank == j)
                if j > 0:
                    nll = F.binary_cross_entropy_with_logits(logits_of(x_t, dna, atac), target, reduction="none")
                joint += (nll * current).sum((1, 2))
                x_t = torch.where(current, x0, x_t)
            per_region[s].append(
                torch.stack(
                    [marginal, joint, (ref_nll * hidden).sum((1, 2)), hidden.sum((1, 2)).float()], -1
                ).cpu()
            )
        region_ids += [f"{batch_idx}:{i}" for i in range(dna.shape[0])]
        print(f"batch {batch_idx}", flush=True)

out = {
    "checkpoint": args.checkpoint,
    "h5": args.h5,
    "max_regions": args.max_regions,
    "span_lengths": args.span_lengths,
    "columns": ["marginal_nll_sum", "joint_nll_sum", "reference_nll_sum", "hidden_calls"],
    "region_ids": region_ids,
    "per_region": {str(s): torch.cat(v).tolist() for s, v in per_region.items()},
}
json.dump(out, open(args.out, "w"))
for s, v in per_region.items():
    t = torch.cat(v).sum(0)
    print(f"span {s:4d} bp: marginal {t[0] / t[3]:.4f}  joint {t[1] / t[3]:.4f}  reference {t[2] / t[3]:.4f}  nats/call over {int(t[3])} calls")
