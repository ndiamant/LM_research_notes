"""Central-1 kb loss at CRE vs other sites: the 16 kb additive model against the 2 kb additive
model (whose SMF window is that 1 kb), with the 16 kb pilot and the original 2 kb model as
references, on the same reads and masks.

Reads: the 16 kb runs' val items (chr8, chr20; reads covering >= 75% of the 16,384 bp SMF
window; all molecules). Each read gets one noise level and one mask draw, shared by every
model; the additive models also get each read's b_m. Scored: weighted masked cross-entropy
(the training loss) over the central 1,024 bp, split by the IDR-peak mask (ENCFF614SMH) into
CRE and other sites, and by the region's class. "add16_centre_only" is the 16 kb additive
model with every call outside the central 1 kb hidden: the same read calls as the 2 kb
models, but 16 kb of DNA and its own training. Per-region sums go to --out (npz); 95%
intervals of the relative differences are from 2,000 region bootstraps.
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
ap.add_argument("--add16", required=True)
ap.add_argument("--pilot16", required=True)
ap.add_argument("--add2", required=True)
ap.add_argument("--orig2", required=True)
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
    cre_bed=F_DIR + "/inputs/ENCFF614SMH.bed.gz",
)
data = FigwigSMFDataModule(cfg)
data.setup("validate")
dataset = data.val_dataset
classes = {}
for line in open(F_DIR + "/kept_windows.bed"):
    f = line.split("\t")
    classes[(f[0], int(f[1]))] = f[4]

device = "cuda"
models = {k: load_model(getattr(args, k)).to(device).eval() for k in ("add16", "pilot16", "add2", "orig2")}
assert models["add16"].cfg.emission_form == "additive" and models["add2"].cfg.emission_form == "additive"
assert models["add2"].cfg.input_length == 2048 and models["add16"].cfg.smf_length == 16384
schedule = models["pilot16"].schedule
L_IN, FLANK = 18432, 1024
MID = L_IN // 2
centre = slice(MID - 512 - FLANK, MID + 512 - FLANK)  # central 1 kb, 16 kb SMF coordinates
short_dna = slice(MID - 1024, MID + 1024)


def terms(model, x0, dna, t, draw, efficiency, keep_only=None):
    valid = x0 != -1
    keep = valid & (draw < schedule.alpha_t(t.float())[..., None])
    if keep_only is not None:
        keep &= keep_only
    x_t = torch.where(keep, x0, torch.full_like(x0, MASK_TOKEN))
    with torch.autocast("cuda", dtype=torch.bfloat16):
        logits = model._predict(model._encode_smf_tokens(x_t), dna, efficiency, dna=dna)
    bce = F.binary_cross_entropy_with_logits(
        logits.float(), x0.float().clamp_min(0).unsqueeze(2), reduction="none"
    ).squeeze(2)
    return bce * schedule.cross_entropy_loss_weight(t) * (valid & (x_t == MASK_TOKEN)), valid


KEYS = ("add16", "add16_centre_only", "pilot16", "add2", "orig2")
rec = {f"{k}_{s}": [] for k in KEYS for s in ("cre", "bg")}
rec.update(calls_cre=[], calls_bg=[], active=[])
only = torch.zeros(L_IN - 2 * FLANK, dtype=torch.bool, device=device)
only[centre] = True
n_done = 0
with torch.no_grad():
    for batch_idx, batch in enumerate(data.val_dataloader()):
        smf = batch.smf.to(device).long()
        dna = batch.dna.to(device).float()
        eff = batch.efficiency.to(device)
        x0 = smf[..., FLANK : L_IN - FLANK]
        cre = batch.cre.to(device)[..., FLANK : L_IN - FLANK][..., centre][:, None]  # B x 1 x 1024
        bsz, reads, _ = x0.shape
        g = torch.Generator(device=device).manual_seed(batch_idx)
        base = torch.linspace(0, 1 - 1 / reads, steps=reads, device=device)
        t = (base + torch.rand(bsz, 1, device=device, generator=g)) % 1
        draw = torch.rand(x0.shape, device=device, generator=g)
        out = {
            "add16": terms(models["add16"], x0, dna, t, draw, eff),
            "add16_centre_only": terms(models["add16"], x0, dna, t, draw, eff, keep_only=only),
            "pilot16": terms(models["pilot16"], x0, dna, t, draw, None),
        }
        out = {k: (v[0][..., centre], v[1][..., centre]) for k, v in out.items()}
        xc, dc, drc = x0[..., centre], dna[..., short_dna], draw[..., centre]
        out["add2"] = terms(models["add2"], xc, dc, t, drc, eff)
        out["orig2"] = terms(models["orig2"], xc, dc, t, drc, None)
        valid = out["pilot16"][1]
        for k in KEYS:
            loss = out[k][0]
            assert torch.equal(out[k][1], valid)
            rec[f"{k}_cre"].append((loss * cre).sum((1, 2)).cpu())
            rec[f"{k}_bg"].append((loss * ~cre).sum((1, 2)).cpu())
        rec["calls_cre"].append((valid & cre).sum((1, 2)).cpu())
        rec["calls_bg"].append((valid & ~cre).sum((1, 2)).cpu())
        for i in range(bsz):
            idx = n_done + i
            kept_start = int(dataset.starts[idx]) + 18498 // 2 - 1057
            rec["active"].append(torch.tensor([classes.get((dataset.chroms[idx], kept_start)) == "active"]))
        n_done += bsz
        if batch_idx % 100 == 0:
            print(f"{n_done} regions", flush=True)

r = {k: torch.cat(v).double().numpy() for k, v in rec.items()}
r["active"] = r["active"].astype(bool)
np.savez(args.out, **r)
rng = np.random.default_rng(0)
boot = rng.integers(0, len(r["active"]), (2000, len(r["active"])))


def loss(key, part, sel):
    return r[f"{key}_{part}"][sel].sum() / r[f"calls_{part}"][sel].sum()


def rel_ci(a, b, part, sel):
    idx = np.flatnonzero(sel)
    vals = []
    for row in boot[:, : len(idx)]:
        s = idx[row % len(idx)]
        vals.append(r[f"{a}_{part}"][s].sum() / r[f"{b}_{part}"][s].sum() - 1)
    return np.percentile(vals, [2.5, 97.5]) * 100


print(f"\n{n_done} val regions; central 1 kb; loss per observed call")
for label, sel in (("all regions", np.ones_like(r["active"])), ("active regions", r["active"]),
                   ("background regions", ~r["active"])):
    for part, name in (("cre", "CRE (IDR peak) sites"), ("bg", "other sites")):
        if r[f"calls_{part}"][sel].sum() == 0:
            continue
        line = f"{label:19s} {name:21s} calls {int(r[f'calls_{part}'][sel].sum()):>10,d} |"
        for k in KEYS:
            line += f" {k} {loss(k, part, sel):.4f}"
        ref = loss("add2", part, sel)
        for k in ("add16", "add16_centre_only"):
            lo, hi = rel_ci(k, "add2", part, sel)
            line += f" | {k} vs add2 {100 * (loss(k, part, sel) / ref - 1):+.1f}% [{lo:+.1f}, {hi:+.1f}]"
        print(line)
