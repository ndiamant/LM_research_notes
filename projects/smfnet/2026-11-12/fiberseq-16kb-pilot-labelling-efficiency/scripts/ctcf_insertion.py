"""In-silico insertion of two CTCF elements, convergent and controls, into background windows.

Till's suggestion: does the 16 kb model learn interactions at > 5 kb? Each element is three
copies of the CTCF consensus (JASPAR MA0139.1, TGGCCACCAGGGGGCGCTA) at offsets 0, 50 and
100 bp, as in Mola (Wang et al. 2026). The two elements are centred in the 16,384 bp SMF
window, d = 2, 5, 8, 12 kb apart. Orientations: convergent (forward motif on the left
element, reverse complement on the right, both pointing inward), divergent, tandem (both
forward), each element alone, and no insertion. Backgrounds: val windows (chr8, chr20) of
the H5's GC-matched background class.

Per sequence and orientation: 256 generated reads; accessibility at each element (+/-100 bp
around its centre, excluding nothing; generated reads are called at every A/T); the raw
correlation of the two elements' accessibility across reads, and the partial correlation
given each read's labelling level (its accessibility over the window excluding +/-300 bp
around both element positions); and the read-level-adjusted co-accessibility map at 100 bp
bins, summed over backgrounds. Writes results.jsonl, maps.npz and the reads of the first
backgrounds for plotting.
"""

import argparse
import json
import os
import time

import numpy as np
import torch

from smf_net.cli import load_model_for_inference
from smf_net.figwig_data import load_region_dna, read_bed_regions

ap = argparse.ArgumentParser()
ap.add_argument("--checkpoint", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--n-backgrounds", type=int, default=50)
ap.add_argument("--reads", type=int, default=256)
ap.add_argument("--keep-reads-for", type=int, default=3, help="backgrounds whose reads are saved")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
FASTA = "/scratch/groups/btrippe/ndiamant/reference/hg38.fa"
MOTIF = "TGGCCACCAGGGGGCGCTA"
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
RC = str.maketrans("ACGT", "TGCA")
L_IN, FLANK = 18432, 1024
L_SMF = L_IN - 2 * FLANK
OFFSETS = (0, 50, 100)
ELEMENT = OFFSETS[-1] + len(MOTIF)  # 119 bp
DISTANCES = (2000, 5000, 8000, 12000)
BIN, HALF, EXCL = 100, 100, 300

model = load_model_for_inference(args.checkpoint, "cuda").eval()
assert model.cfg.input_length == L_IN and model.cfg.smf_length == L_SMF

classes = {}
for line in open(S + "/kept_windows.bed"):
    f = line.split("\t")
    classes[(f[0], int(f[1]))] = f[4]
rng = np.random.default_rng(0)
background = [(c, s + 1057 - L_IN // 2) for c, s in read_bed_regions(S + "/kept_windows.bed", 2114)
              if c in ("chr8", "chr20") and classes[(c, s)] == "background"]
background = [background[i] for i in rng.choice(len(background), args.n_backgrounds, replace=False)]
dna_bg = load_region_dna(FASTA, background, L_IN)  # [n, L_IN] int8


def element(orientation):
    motif = MOTIF if orientation == "F" else MOTIF.translate(RC)[::-1]
    seq = ["N"] * ELEMENT
    for o in OFFSETS:
        seq[o : o + len(MOTIF)] = motif
    return seq


def insert(dna, centre_smf, orientation):
    """Insert an element centred at SMF-window coordinate centre_smf; N positions keep the background."""
    start = FLANK + centre_smf - ELEMENT // 2
    for k, base in enumerate(element(orientation)):
        if base != "N":
            dna[start + k] = CODE[base]


conditions = [("none", None, None, None)]
for d in DISTANCES:
    left, right = L_SMF // 2 - d // 2, L_SMF // 2 + d // 2
    conditions += [
        (f"convergent_{d}", d, ("F", left), ("R", right)),
        (f"divergent_{d}", d, ("R", left), ("F", right)),
        (f"tandem_{d}", d, ("F", left), ("F", right)),
        (f"left_only_{d}", d, ("F", left), None),
        (f"right_only_{d}", d, None, ("R", right)),
    ]


def accessible_fraction(reads, lo, hi):
    c = reads[:, lo:hi]
    sites = c != 2
    with np.errstate(invalid="ignore", divide="ignore"):
        return (sites & (c == 0)).sum(1) / sites.sum(1)


def binned(reads):
    n = L_SMF // BIN
    c = reads[:, : n * BIN].reshape(len(reads), n, BIN)
    sites = (c != 2).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(sites > 0, ((c == 0) & (c != 2)).sum(-1) / sites, np.nan)


def adjust(x):
    m = np.nanmean(x, axis=1, keepdims=True)
    ok = ~np.isnan(x)
    n = ok.sum(0).clip(1)
    xbar = np.where(ok, x, 0).sum(0) / n
    mbar = np.where(ok, m, 0).sum(0) / n
    cov = np.where(ok, (x - xbar) * (m - mbar), 0).sum(0) / n
    var = np.where(ok, (m - mbar) ** 2, 0).sum(0) / n
    slope = np.where(var > 0, cov / np.where(var > 0, var, 1), 0)
    return np.where(ok, x - xbar - slope * (m - mbar), np.nan)


def corr(x):
    x = np.where(np.isnan(x), np.nanmean(x, axis=0, keepdims=True), x)
    x = x - x.mean(0)
    sd = np.sqrt((x * x).mean(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        return (x.T @ x / len(x)) / np.outer(sd, sd)


maps = {name: np.zeros((L_SMF // BIN, L_SMF // BIN)) for name, *_ in conditions}
kept = {}
t0 = time.time()
with open(f"{args.out}/results.jsonl", "w") as out:
    for b in range(len(background)):
        seqs = []
        for name, d, left, right in conditions:
            dna = dna_bg[b].clone()
            for part in (left, right):
                if part is not None:
                    insert(dna, part[1], part[0])
            seqs.append(dna)
        dna = torch.stack(seqs)
        reads = []
        for chunk in range(0, len(dna), 7):  # 7 sequences x 128 reads per forward pass
            parts = []
            for rep in range(args.reads // 128):
                g = torch.Generator(device="cuda").manual_seed(1000 * b + 10 * chunk + rep)
                parts.append(model.sample(dna[chunk : chunk + 7], reads_per_region=128, generator=g))
            reads.append(torch.cat(parts, dim=1).to(torch.int8).cpu().numpy())
        reads = np.concatenate(reads)  # [conditions, reads, L_SMF]
        for (name, d, left, right), r in zip(conditions, reads, strict=True):
            row = {"background": f"{background[b][0]}:{background[b][1]}", "condition": name, "distance": d}
            if d is not None:
                pl, pr = L_SMF // 2 - d // 2, L_SMF // 2 + d // 2
                x = accessible_fraction(r, pl - HALF, pl + HALF)
                y = accessible_fraction(r, pr - HALF, pr + HALF)
                rest = r.copy()
                for c in (pl, pr):
                    rest[:, c - EXCL : c + EXCL] = 2
                level = accessible_fraction(rest, 0, L_SMF)
                design = np.column_stack([np.ones(len(level)), level])
                rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
                ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
                row.update(left_access=float(x.mean()), right_access=float(y.mean()),
                           raw_r=float(np.corrcoef(x, y)[0, 1]), partial_r=float(np.corrcoef(rx, ry)[0, 1]))
            maps[name] += np.nan_to_num(corr(adjust(binned(r))))
            out.write(json.dumps(row) + "\n")
        if b < args.keep_reads_for:
            kept[f"bg{b}"] = reads
        print(f"background {b + 1}/{len(background)} done, {time.time() - t0:.0f} s", flush=True)

np.savez_compressed(f"{args.out}/maps.npz", **{k: v / len(background) for k, v in maps.items()})
np.savez_compressed(f"{args.out}/reads_first_backgrounds.npz",
                    conditions=np.array([c[0] for c in conditions]), **kept)
print("done", time.time() - t0)
