"""Loss of the 16 kb and 2 kb models on the same sites: the central 1 kb of val regions.

Both models score the same reads (the 16 kb run's val items: reads covering >= 75% of the
16,384 bp SMF window), with the same noise level per read and the same masks at every
central site; only the context differs. The 16 kb model sees 18,432 bp of DNA and the
read's calls over 16 kb (masked at the same rate); the 2 kb model sees the central
2,048 bp of DNA and only the central 1 kb of calls (its SMF window). The loss is the
model's masked-diffusion loss (weighted cross-entropy on masked sites over observed sites),
restricted to the central 1 kb; the 16 kb model's whole-window loss is reported too.
Results per region go to --out (JSON lines) for stratifying by class and accessibility.
"""

import argparse
import json
import os

import numpy as np
import torch
import torch.nn.functional as F

from smf_net.diffusion_schedule import MASK_TOKEN
from smf_net.figwig_data import FigwigDataModuleConfig, FigwigSMFDataModule
from smf_net.model import load_model

ap = argparse.ArgumentParser()
ap.add_argument("--long", required=True, help="16 kb model checkpoint")
ap.add_argument("--short", required=True, help="2 kb model checkpoint")
ap.add_argument("--out", required=True)
ap.add_argument("--batch-size", type=int, default=16)
ap.add_argument("--max-regions", type=int, default=None)
ap.add_argument(
    "--mask-outside-centre", action="store_true",
    help="also score the 16 kb model with every call outside the central 1 kb masked, so it "
    "sees the same read calls as the 2 kb model and differs only in DNA context and training",
)
args = ap.parse_args()

F_DIR = os.environ["SCRATCH"] + "/fiberseq/GM12878"
cfg = FigwigDataModuleConfig(
    fasta_path="/scratch/groups/btrippe/ndiamant/reference/hg38.fa",
    figwig_dir=F_DIR + "/figwig", regions_bed=F_DIR + "/kept_windows.bed",
    assay_type="fiberseq_m6a", region_length=18498, model_input_length=18432, dna_flank=1024,
    min_read_fraction=0.75, min_overlap=500, min_spanning_molecules=32,
    val_chrm=["chr8", "chr20"], batch_size=args.batch_size, num_workers=4,
    max_regions_per_split=args.max_regions,
)
data = FigwigSMFDataModule(cfg)
data.setup("validate")
dataset = data.val_dataset
classes = {}
with open(F_DIR + "/kept_windows.bed") as handle:
    for line in handle:
        f = line.split("\t")
        classes[(f[0], int(f[1]))] = f[4]

device = "cuda"
long_model = load_model(args.long).to(device).eval()
short_model = load_model(args.short).to(device).eval()
assert long_model.cfg.smf_length == 16384 and short_model.cfg.input_length == 2048
schedule = long_model.schedule
L_IN, FLANK = 18432, 1024
MID = L_IN // 2  # 9216
centre = slice(MID - 512 - FLANK, MID + 512 - FLANK)  # central 1 kb in 16 kb SMF coordinates
short_dna = slice(MID - 1024, MID + 1024)  # the 2 kb model's DNA input


def weighted_terms(model, x0, dna, t, draw, keep_only=None):
    """Per-position weighted cross-entropy on masked sites, and the observed-site mask.
    ``keep_only`` (a position mask) masks every call outside it."""
    valid = x0 != -1
    keep = valid & (draw < schedule.alpha_t(t.float())[..., None])
    if keep_only is not None:
        keep &= keep_only
    x_t = torch.where(keep, x0, torch.full_like(x0, MASK_TOKEN))
    with torch.autocast("cuda", dtype=torch.bfloat16):
        logits = model.model(model._encode_smf_tokens(x_t), dna)
    bce = F.binary_cross_entropy_with_logits(
        logits.float(), x0.float().clamp_min(0).unsqueeze(2), reduction="none"
    ).squeeze(2)
    w = schedule.cross_entropy_loss_weight(t)
    return bce * w * (valid & (x_t == MASK_TOKEN)), valid


loader = data.val_dataloader()
n_done = 0
with open(args.out, "w") as out, torch.no_grad():
    for batch_idx, batch in enumerate(loader):
        smf = batch.smf.to(device).long()  # B x R x 18432
        dna = batch.dna.to(device).float()  # B x 4 x 18432
        x0 = smf[..., FLANK : L_IN - FLANK]  # 16 kb SMF window
        bsz, reads, _ = x0.shape
        g = torch.Generator(device=device).manual_seed(batch_idx)
        base = torch.linspace(0, 1 - 1 / reads, steps=reads, device=device)
        t = (base + torch.rand(bsz, 1, device=device, generator=g)) % 1
        draw = torch.rand(x0.shape, device=device, generator=g)

        long_terms, long_valid = weighted_terms(long_model, x0, dna, t, draw)
        short_terms, short_valid = weighted_terms(
            short_model, x0[..., centre], dna[..., short_dna], t, draw[..., centre]
        )
        assert torch.equal(short_valid, long_valid[..., centre])
        if args.mask_outside_centre:
            only = torch.zeros(x0.shape[-1], dtype=torch.bool, device=device)
            only[centre] = True
            dna_only_terms, _ = weighted_terms(long_model, x0, dna, t, draw, keep_only=only)
        for i in range(bsz):
            idx = n_done + i
            chrom, start = dataset.chroms[idx], int(dataset.starts[idx])
            kept_start = start + 18498 // 2 - 1057
            c = x0[i][..., centre]
            observed = c != -1
            out.write(json.dumps({
                "region": dataset.regions[idx],
                "class": classes.get((chrom, kept_start), "?"),
                "centre_calls": int(observed.sum()),
                "centre_accessible": float((c == 0).sum() / observed.sum().clamp_min(1)),
                "long_centre_sum": float(long_terms[i][..., centre].sum()),
                "long_dna_only_centre_sum": (
                    float(dna_only_terms[i][..., centre].sum()) if args.mask_outside_centre else None
                ),
                "short_centre_sum": float(short_terms[i].sum()),
                "long_window_sum": float(long_terms[i].sum()),
                "window_calls": int(long_valid[i].sum()),
            }) + "\n")
        n_done += bsz
        if batch_idx % 100 == 0:
            print(f"{n_done} regions", flush=True)

rows = [json.loads(line) for line in open(args.out)]


def loss(rows, key, calls):
    return sum(r[key] for r in rows) / sum(r[calls] for r in rows)


def summary(rows, label):
    if not rows:
        return
    lc = loss(rows, "long_centre_sum", "centre_calls")
    sc = loss(rows, "short_centre_sum", "centre_calls")
    line = (f"{label:32s} n={len(rows):6d}  centre 1 kb: 16 kb model {lc:.4f}, 2 kb model {sc:.4f} "
            f"({(lc - sc) / sc:+.1%})")
    if rows[0]["long_dna_only_centre_sum"] is not None:
        dc = loss(rows, "long_dna_only_centre_sum", "centre_calls")
        line += f"; 16 kb model, calls outside the centre masked {dc:.4f} ({(dc - sc) / sc:+.1%})"
    print(line)


print(f"\n{len(rows)} val regions; 16 kb model whole-window loss "
      f"{loss(rows, 'long_window_sum', 'window_calls'):.4f}")
summary(rows, "all")
for cls in ("active", "background"):
    summary([r for r in rows if r["class"] == cls], cls)
acc = np.array([r["centre_accessible"] for r in rows])
edges = np.quantile(acc, [0, 0.25, 0.5, 0.75, 0.9, 1.0])
for lo, hi in zip(edges[:-1], edges[1:]):
    summary([r for r in rows if lo <= r["centre_accessible"] <= hi],
            f"centre accessible {lo:.3f}-{hi:.3f}")
