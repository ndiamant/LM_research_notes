"""Same-sites loss of the spec run (emission head), the 16 kb pilot and the 2 kb model.

All three score the same reads: the pilot's val items (chr8, chr20; reads covering >= 75% of
the 16,384 bp SMF window), with no molecule-rate filter, so the low-labelling molecules the
spec run never trained on are included. Same noise level per read and the same masks for
every model. The spec model also gets each read's efficiency b_m (the others ignore it).
Scores: the central 1 kb (all three) and the whole 16 kb window (the two 16 kb models).
Per read, saved to --out (npz): each model's weighted masked cross-entropy sums, observed
calls, b_m and background-site m6A rate, and the region's class and centre accessibility.
"""

import argparse
import os

import numpy as np
import torch
import torch.nn.functional as F

from smf_net.diffusion_schedule import MASK_TOKEN
from smf_net.figwig_data import FigwigDataModuleConfig, FigwigSMFDataModule
from smf_net.model import load_model

ap = argparse.ArgumentParser()
ap.add_argument("--spec", required=True)
ap.add_argument("--pilot", required=True)
ap.add_argument("--short", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--batch-size", type=int, default=16)
ap.add_argument("--max-regions", type=int, default=None)
args = ap.parse_args()

F_DIR = os.environ["SCRATCH"] + "/fiberseq/GM12878"
cfg = FigwigDataModuleConfig(
    fasta_path="/scratch/groups/btrippe/ndiamant/reference/hg38.fa",
    figwig_dir=F_DIR + "/figwig", regions_bed=F_DIR + "/kept_windows.bed",
    assay_type="fiberseq_m6a", region_length=18498, model_input_length=18432, dna_flank=1024,
    min_read_fraction=0.75, min_overlap=500, min_spanning_molecules=32,
    val_chrm=["chr8", "chr20"], batch_size=args.batch_size, num_workers=6,
    max_regions_per_split=args.max_regions, efficiency_dir=F_DIR + "/efficiency",
)
data = FigwigSMFDataModule(cfg)
data.setup("validate")
dataset = data.val_dataset
classes = {}
for line in open(F_DIR + "/kept_windows.bed"):
    f = line.split("\t")
    classes[(f[0], int(f[1]))] = f[4]

device = "cuda"
models = {
    "spec": load_model(args.spec).to(device).eval(),
    "pilot": load_model(args.pilot).to(device).eval(),
    "short": load_model(args.short).to(device).eval(),
}
assert models["spec"].cfg.emission_head and models["short"].cfg.input_length == 2048
schedule = models["pilot"].schedule
L_IN, FLANK = 18432, 1024
MID = L_IN // 2
centre = slice(MID - 512 - FLANK, MID + 512 - FLANK)
short_dna = slice(MID - 1024, MID + 1024)


def terms(model, x0, dna, t, draw, efficiency):
    valid = x0 != -1
    keep = valid & (draw < schedule.alpha_t(t.float())[..., None])
    x_t = torch.where(keep, x0, torch.full_like(x0, MASK_TOKEN))
    with torch.autocast("cuda", dtype=torch.bfloat16):
        logits = model._predict(model._encode_smf_tokens(x_t), dna, efficiency, dna=dna)
    bce = F.binary_cross_entropy_with_logits(
        logits.float(), x0.float().clamp_min(0).unsqueeze(2), reduction="none"
    ).squeeze(2)
    return bce * schedule.cross_entropy_loss_weight(t) * (valid & (x_t == MASK_TOKEN)), valid


def molecules_of(idx):
    """The val item's molecules, drawn as FigwigSMFDataset.get_item draws them."""
    offset = dataset.centred_offset(L_IN)
    molecules, _ = dataset._window(idx, offset, L_IN)
    g = torch.Generator()
    g.manual_seed(idx)
    chosen = torch.multinomial(torch.ones(len(molecules)), dataset.reads_per_region,
                               replacement=len(molecules) < dataset.reads_per_region, generator=g)
    return molecules[chosen.numpy()]


rec = {k: [] for k in ("spec_c", "pilot_c", "short_c", "calls_c", "spec_w", "pilot_w", "calls_w",
                       "b", "rate", "region", "active", "centre_acc")}
n_done = 0
with torch.no_grad():
    for batch_idx, batch in enumerate(data.val_dataloader()):
        smf = batch.smf.to(device).long()
        dna = batch.dna.to(device).float()
        eff = batch.efficiency.to(device)
        x0 = smf[..., FLANK : L_IN - FLANK]
        bsz, reads, _ = x0.shape
        g = torch.Generator(device=device).manual_seed(batch_idx)
        base = torch.linspace(0, 1 - 1 / reads, steps=reads, device=device)
        t = (base + torch.rand(bsz, 1, device=device, generator=g)) % 1
        draw = torch.rand(x0.shape, device=device, generator=g)
        spec, valid = terms(models["spec"], x0, dna, t, draw, eff)
        pilot, _ = terms(models["pilot"], x0, dna, t, draw, None)
        short, _ = terms(models["short"], x0[..., centre], dna[..., short_dna], t, draw[..., centre], None)
        rec["spec_c"].append(spec[..., centre].sum(-1).cpu()); rec["pilot_c"].append(pilot[..., centre].sum(-1).cpu())
        rec["short_c"].append(short.sum(-1).cpu()); rec["calls_c"].append(valid[..., centre].sum(-1).cpu())
        rec["spec_w"].append(spec.sum(-1).cpu()); rec["pilot_w"].append(pilot.sum(-1).cpu())
        rec["calls_w"].append(valid.sum(-1).cpu()); rec["b"].append(eff.cpu())
        for i in range(bsz):
            idx = n_done + i
            lanes = dataset.lanes[dataset.chroms[idx]]
            rec["rate"].append(torch.from_numpy(lanes.background_rate[molecules_of(idx)])[None])
            kept_start = int(dataset.starts[idx]) + 18498 // 2 - 1057
            c = x0[i][..., centre]
            obs = c != -1
            rec["region"].append(torch.full((1, reads), idx))
            rec["active"].append(torch.full((1, reads), classes.get((dataset.chroms[idx], kept_start)) == "active"))
            rec["centre_acc"].append(torch.full((1, reads), float((c == 0).sum() / obs.sum().clamp_min(1))))
        n_done += bsz
        if batch_idx % 100 == 0:
            print(f"{n_done} regions", flush=True)

out = {k: torch.cat(v).flatten().numpy() for k, v in rec.items()}
np.savez(args.out, **out)


def L(sel, key, calls):
    return out[key][sel].sum() / out[calls][sel].sum()


def report(sel, label):
    s, p, q = L(sel, "spec_c", "calls_c"), L(sel, "pilot_c", "calls_c"), L(sel, "short_c", "calls_c")
    sw, pw = L(sel, "spec_w", "calls_w"), L(sel, "pilot_w", "calls_w")
    print(f"{label:34s} reads {sel.sum():7d} | centre 1 kb: spec {s:.4f} pilot {p:.4f} ({(s - p) / p:+.1%}) "
          f"2 kb {q:.4f} ({(s - q) / q:+.1%}) | 16 kb window: spec {sw:.4f} pilot {pw:.4f} ({(sw - pw) / pw:+.1%})")


every = np.ones(len(out["b"]), bool)
print(f"\n{n_done} val regions")
report(every, "all reads")
report(out["rate"] >= 0.02, "reads with rate >= 0.02 (trained on)")
report(out["rate"] < 0.02, "reads with rate < 0.02 (filtered)")
report(out["active"], "active regions")
report(~out["active"], "background regions")
edges = np.quantile(out["b"], [0, 0.2, 0.4, 0.6, 0.8, 1])
for k in range(5):
    report((out["b"] >= edges[k]) & (out["b"] <= edges[k + 1]), f"b quintile {k + 1} ({edges[k]:+.2f} to {edges[k + 1]:+.2f})")
