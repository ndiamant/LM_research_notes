"""Per-molecule labelling efficiency b_m from logit P(m6A at site s on molecule m) = a_s + b_m.

a_s, the per-site intercept, comes from the pileup of all molecules at the site (smoothed
logit of its m6A fraction): it absorbs the region's activity and the enzyme/caller's
sequence-context bias. b_m is then each molecule's offset from the other molecules at the
same sites, fitted on its background sites only: sites outside +/-500 bp of IDR ATAC peaks
(ENCFF614SMH) and of active cCREs, and with >= 10 molecules calling them. So b_m carries
labelling efficiency (and any molecule-wide state) but not which elements the molecule has
open, nor how active its region is.

The fit uses a per-molecule histogram of a_s over its background sites (BINS logit bins):
b_m solves sum_k n_k sigma(a_k + b) = m6A calls. Two passes over each chromosome's lanes,
in CHUNK-base pieces. Writes <out>/<chrom>.efficiency.npz aligned to the molecule table.

Usage: molecule_efficiency.py <figwig_dir> <out_dir> chrom [chrom ...]
"""

import gzip
import os
import sys
import time
from collections import defaultdict

import numpy as np

from smf_net.figwig_data import MoleculeLanes

figwig_dir, out_dir, chroms = sys.argv[1], sys.argv[2], sys.argv[3:]
os.makedirs(out_dir, exist_ok=True)
S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
CHUNK, MARGIN, MIN_DEPTH = 1 << 18, 500, 10
LOGIT_LO, LOGIT_HI, BINS = -8.0, 4.0, 96
centres = LOGIT_LO + (np.arange(BINS) + 0.5) * (LOGIT_HI - LOGIT_LO) / BINS

elements = defaultdict(list)
for line in gzip.open(S + "/inputs/ENCFF614SMH.bed.gz", "rt"):
    f = line.split("\t")
    elements[f[0]].append((int(f[1]), int(f[2])))
for line in open("/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/candidates.bed"):
    f = line.split("\t")
    if f[4] == "active":
        c = (int(f[1]) + int(f[2])) // 2
        elements[f[0]].append((c, c + 1))


def solve(hist, m6a, iters=50):
    """b with sum_k hist[:, k] sigma(centres[k] + b) = m6a, by bisection on [-12, 12]."""
    lo = np.full(len(hist), -12.0)
    hi = np.full(len(hist), 12.0)
    for _ in range(iters):
        mid = (lo + hi) / 2
        expected = (hist / (1 + np.exp(-(centres[None, :] + mid[:, None])))).sum(1)
        too_low = expected < m6a
        lo = np.where(too_low, mid, lo)
        hi = np.where(too_low, hi, mid)
    return (lo + hi) / 2


for chrom in chroms:
    t0 = time.time()
    lanes = MoleculeLanes(f"{figwig_dir}/{chrom}.fw")
    L = lanes.chrom_length
    n_lanes = int(lanes.lane.max()) + 1
    background = np.ones(L, bool)
    for s, e in elements[chrom]:
        background[max(s - MARGIN, 0) : e + MARGIN] = False

    # pass 1: per-site pileup -> binned intercept (255 = too shallow or no call)
    a_bin = np.full(L, 255, np.uint8)
    for c0 in range(0, L, CHUNK):
        w = min(CHUNK, L - c0)
        v = lanes.reader.read(chrom, [c0], w, tracks=list(range(n_lanes)), n_jobs=4)[0]
        acc = (v == 1).sum(0)
        called = acc + (v == 2).sum(0)
        ok = called >= MIN_DEPTH
        logit = np.log((acc + 0.5) / (called - acc + 0.5))
        k = np.clip(((logit - LOGIT_LO) / (LOGIT_HI - LOGIT_LO) * BINS).astype(int), 0, BINS - 1)
        a_bin[c0 : c0 + w] = np.where(ok, k, 255)
    t1 = time.time()

    # pass 2: per-molecule histogram of a_s over background sites, and m6A calls there
    n = len(lanes)
    hist = np.zeros((n, BINS), np.int32)
    m6a = np.zeros(n, np.int64)
    for c0 in range(0, L, CHUNK):
        w = min(CHUNK, L - c0)
        lo = np.searchsorted(lanes.start, c0 - lanes.max_length)
        hi = np.searchsorted(lanes.start, c0 + w)
        idx = np.arange(lo, hi)
        idx = idx[lanes.end[idx] > c0]
        if len(idx) == 0:
            continue
        v = lanes.reader.read(chrom, [c0], w, tracks=list(range(n_lanes)), n_jobs=4)[0]
        bg = background[c0 : c0 + w]
        ab = a_bin[c0 : c0 + w]
        for i in idx:
            a, z = max(lanes.start[i], c0) - c0, min(lanes.end[i], c0 + w) - c0
            seg = v[lanes.lane[i], a:z]
            use = (seg > 0) & bg[a:z] & (ab[a:z] != 255)
            if use.any():
                hist[i] += np.bincount(ab[a:z][use], minlength=BINS)
                m6a[i] += int((seg[use] == 1).sum())
    t2 = time.time()

    n_sites = hist.sum(1)
    b = np.full(n, np.nan, np.float32)
    fit = n_sites >= 500
    b[fit] = solve(hist[fit].astype(float), m6a[fit])
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = (m6a / n_sites).astype(np.float32)
    np.savez(f"{out_dir}/{chrom}.efficiency.npz", efficiency=b, background_sites=n_sites.astype(np.int32),
             background_m6a=m6a.astype(np.int32), background_rate=rate)
    print(f"{chrom}: {n:,} molecules, {fit.sum():,} fitted ({fit.mean():.1%}); "
          f"b median {np.nanmedian(b):+.2f}, IQR {np.nanpercentile(b, 25):+.2f} to {np.nanpercentile(b, 75):+.2f}; "
          f"{L / 1e6:.0f} Mb in {t1 - t0:.0f} + {t2 - t1:.0f} s", flush=True)
