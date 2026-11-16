"""Does the additive emission head restore between-molecule (labelling) variation in generated
reads? 2 kb models, the first --regions val regions (chr8, chr20; the metrics' 1,000).

Per region, R observed reads (the val item drawn with --reads, seed = region index) and, row
for row, generated reads from: the additive model with each observed read's b_m ("matched"),
the additive model at b_m = 0 ("fixed0"), and the original model ("original"). Each generated
read takes its observed read's no-call positions. Per read: mean accessibility (m6A share of
called sites in the SMF window). Per region and set: sd of read means, its ratio to observed,
and Pearson(read mean, b_m). Co-accessibility: 100 bp bins, Pearson across reads between bins
500-1,000 bp apart, raw and adjusted for the read's mean (per-bin regression, as the pilot
note). Writes per-region rows to --out (npz) and prints medians/means.
"""
import argparse
import os
import time

import numpy as np
import torch

from smf_net.figwig_data import FigwigDataModuleConfig, FigwigSMFDataModule
from smf_net.model import load_model

ap = argparse.ArgumentParser()
ap.add_argument("--additive", required=True)
ap.add_argument("--original", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--regions", type=int, default=1000)
ap.add_argument("--reads", type=int, default=128)
ap.add_argument("--batch", type=int, default=16)
args = ap.parse_args()

F_DIR = os.environ["SCRATCH"] + "/fiberseq/GM12878"
cfg = FigwigDataModuleConfig(
    fasta_path="/scratch/groups/btrippe/ndiamant/reference/hg38.fa",
    figwig_dir=F_DIR + "/figwig", regions_bed=F_DIR + "/kept_windows.bed",
    assay_type="fiberseq_m6a", region_length=2114, model_input_length=2048,
    min_overlap=500, min_spanning_molecules=32, val_chrm=["chr8", "chr20"],
    reads_per_region=args.reads, max_regions_per_split=args.regions,
    efficiency_dir=F_DIR + "/efficiency", num_workers=0)
data = FigwigSMFDataModule(cfg)
data.setup("validate")
ds = data.val_dataset
FLANK, L = 512, 1024
additive = load_model(args.additive).cuda().eval()
original = load_model(args.original).cuda().eval()
assert additive.cfg.emission_head and additive.cfg.emission_form == "additive"
BIN = 100
NB = L // BIN


def read_means(r):
    called = r != -1
    return np.where(called.sum(1) > 0, (r == 0).sum(1) / np.maximum(called.sum(1), 1), np.nan)


def binned(r):
    c = r[:, : NB * BIN].reshape(len(r), NB, BIN)
    n = (c != -1).sum(-1)
    return np.where(n > 0, (c == 0).sum(-1) / np.maximum(n, 1), np.nan)


def coacc(r, adjust):
    x = binned(r)
    m = read_means(r)[:, None]
    if adjust:
        ok = ~np.isnan(x)
        xm = np.nanmean(x, 0); mm = np.nanmean(np.where(ok, m, np.nan), 0)
        cov = np.nanmean((x - xm) * (m - mm), 0); var = np.nanmean(np.where(ok, (m - mm) ** 2, np.nan), 0)
        x = x - xm - np.where(var > 0, cov / np.where(var > 0, var, 1), 0) * (m - mm)
    x = np.where(np.isnan(x), np.nanmean(x, 0), x)
    x = x - x.mean(0); sd = np.sqrt((x * x).mean(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (x.T @ x / len(x)) / np.outer(sd, sd)
    i, j = np.triu_indices(NB, 1)
    far = (j - i >= 5) & (j - i <= 10)
    return np.nanmean(c[i[far], j[far]])


SETS = ("observed", "matched", "fixed0", "original")
rows = {f"{s}_{k}": [] for s in SETS for k in ("sd", "r_b", "raw", "adj")}
rows["mean_b"] = []
t0 = time.time()
with torch.no_grad():
    for start in range(0, len(ds), args.batch):
        idx = range(start, min(start + args.batch, len(ds)))
        items = [ds[i] for i in idx]
        dna = torch.stack([it.dna for it in items]).cuda()
        obs = torch.stack([it.smf for it in items])[..., FLANK:-FLANK].numpy()
        b = torch.stack([it.efficiency for it in items]).cuda()
        g = torch.Generator(device="cuda").manual_seed(start)
        gen = {
            "matched": additive.sample(dna, args.reads, generator=g, efficiency=b),
            "fixed0": additive.sample(dna, args.reads, generator=g, efficiency=0.0),
            "original": original.sample(dna, args.reads, generator=g),
        }
        gen = {k: v.cpu().numpy() for k, v in gen.items()}
        bn = b.cpu().numpy()
        for n in range(len(items)):
            nocall = obs[n] == -1
            reads = {"observed": obs[n]}
            for k, v in gen.items():
                reads[k] = np.where(nocall, -1, v[n])
            for s, r in reads.items():
                m = read_means(r)
                ok = ~np.isnan(m)
                rows[f"{s}_sd"].append(np.std(m[ok]))
                rows[f"{s}_r_b"].append(np.corrcoef(m[ok], bn[n][ok])[0, 1])
                rows[f"{s}_raw"].append(coacc(r, False))
                rows[f"{s}_adj"].append(coacc(r, True))
            rows["mean_b"].append(bn[n].mean())
        print(f"{idx[-1] + 1} regions, {time.time() - t0:.0f} s", flush=True)

out = {k: np.array(v) for k, v in rows.items()}
np.savez(args.out, **out)
print(f"\n{len(out['mean_b'])} regions x {args.reads} reads")
print(f"{'set':9s} {'sd ratio (median)':>18s} {'r(read mean, b_m)':>18s} {'raw coacc 0.5-1kb':>18s} {'adjusted':>9s}")
for s in SETS:
    ratio = np.nanmedian(out[f"{s}_sd"] / out["observed_sd"])
    print(f"{s:9s} {ratio:18.3f} {np.nanmean(out[f'{s}_r_b']):18.3f} "
          f"{np.nanmean(out[f'{s}_raw']):18.3f} {np.nanmean(out[f'{s}_adj']):9.3f}")
