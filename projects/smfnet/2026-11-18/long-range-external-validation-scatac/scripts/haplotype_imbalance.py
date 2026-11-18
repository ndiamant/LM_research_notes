"""Haplotype imbalance of element accessibility in phased GM12878 Fiber-seq reads.

For each active cCRE, the molecules covering its centre +/- half_width are split by their
haplotype tag (1/2; untagged dropped). The element's accessible fraction is regressed on a
haplotype indicator plus covariates:

- model "b": b_m (per-molecule labelling efficiency) only;
- model "bl": b_m and the molecule's labelling level (accessible fraction over +/- pad,
  excluding +/- exclude around the element). This also absorbs haplotype differences in
  the surrounding domain (e.g. the inactive X), so it isolates element-specific imbalance.

The haplotype coefficient (hap1 - hap2, in accessible-fraction units), its OLS t and
two-sided p are reported, plus the same p with the haplotype labels permuted once within
the element (an empirical null). Each element is annotated with the distance to the
nearest GIAB HG001 heterozygous variant and the coverage of +/- 1 kb by GIAB confident
regions (hets outside them are unknown). HP labels are only consistent within a phase
block, so the sign of the coefficient is arbitrary across elements; only magnitudes and
p-values are comparable.

Usage: python haplotype_imbalance.py <out.tsv> chr20 chrX ...
"""

import gzip
import math
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

from smf_net.evaluation.cre_pairs import accessible_fraction, benjamini_hochberg, read_elements
from smf_net.figwig_data import MoleculeLanes

F = "/scratch/users/diamant/fiberseq/GM12878"
CANDIDATES = "/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/candidates.bed"
VCF = f"{F}/variants/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
CONFIDENT = f"{F}/variants/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
HALF_WIDTH, PAD, EXCLUDE = 100, 4000, 300
MIN_CALLS, MIN_LEVEL_CALLS, MIN_PER_HAP = 10, 200, 10


def het_positions(chroms):
    out = defaultdict(list)
    with gzip.open(VCF, "rt") as handle:
        for line in handle:
            if line[0] == "#":
                continue
            f = line.split("\t", 10)
            if f[0] not in chroms:
                continue
            gt = f[9].split(":", 1)[0].replace("|", "/").split("/")
            if len(gt) == 2 and gt[0] != gt[1]:
                out[f[0]].append(int(f[1]) - 1)
    return {c: np.array(sorted(v)) for c, v in out.items()}


def confident_intervals(chroms):
    out = defaultdict(list)
    with open(CONFIDENT) as handle:
        for line in handle:
            f = line.split("\t")
            if f[0] in chroms:
                out[f[0]].append((int(f[1]), int(f[2])))
    return {c: np.array(sorted(v)).reshape(-1, 2) for c, v in out.items()}


def covered(intervals, lo, hi):
    if intervals is None or not len(intervals):
        return math.nan
    i = np.searchsorted(intervals[:, 1], lo, side="right")
    total = 0
    while i < len(intervals) and intervals[i, 0] < hi:
        total += min(intervals[i, 1], hi) - max(intervals[i, 0], lo)
        i += 1
    return total / (hi - lo)


def hap_test(x, hap1, covariates):
    design = np.column_stack([np.ones(len(x)), hap1, *covariates])
    beta, *_ = np.linalg.lstsq(design, x, rcond=None)
    resid = x - design @ beta
    dof = len(x) - design.shape[1]
    sigma2 = resid @ resid / dof
    try:
        cov = sigma2 * np.linalg.inv(design.T @ design)
    except np.linalg.LinAlgError:
        return math.nan, math.nan, math.nan
    se = math.sqrt(cov[1, 1]) if cov[1, 1] > 0 else math.nan
    t = beta[1] / se if se else math.nan
    return float(beta[1]), se, float(2 * stats.t.sf(abs(t), dof)) if not math.isnan(t) else math.nan


def element_rows(chrom, centres, ctcf, hets, conf, rng):
    lanes = MoleculeLanes(f"{F}/figwig/{chrom}.fw")
    lanes.attach_efficiency(f"{F}/efficiency")
    with np.load(f"{F}/figwig/{chrom}.molecules.npz") as table:
        haplotype = table["haplotype"]
    rows = []
    h = HALF_WIDTH
    for c in centres:
        start, end = c - PAD, c + PAD
        if start < 0 or end > lanes.chrom_length:
            continue
        mol = lanes.overlapping(c - h, c + h, 2 * h)
        if len(mol) < 2 * MIN_PER_HAP:
            continue
        calls = lanes.calls(mol, start, end - start).numpy()
        x, nx = accessible_fraction(calls, PAD - h, PAD + h)
        rest = calls.copy()
        rest[:, PAD - EXCLUDE : PAD + EXCLUDE] = -1
        level, nl = accessible_fraction(rest, 0, end - start)
        hp, b = haplotype[mol], lanes.efficiency[mol]
        ok = (nx >= MIN_CALLS) & (nl >= MIN_LEVEL_CALLS) & (hp > 0) & ~np.isnan(b)
        n1, n2 = int((ok & (hp == 1)).sum()), int((ok & (hp == 2)).sum())
        if min(n1, n2) < MIN_PER_HAP:
            continue
        x, level, b, hap1 = x[ok], level[ok], b[ok], (hp[ok] == 1).astype(float)
        row = {"chrom": chrom, "centre": c, "ctcf": ctcf[c], "n_hap1": n1, "n_hap2": n2,
               "untagged": int(((nx >= MIN_CALLS) & (hp == 0)).sum()),
               "mean_acc": float(x.mean()), "acc_hap1": float(x[hap1 == 1].mean()),
               "acc_hap2": float(x[hap1 == 0].mean()),
               "level_diff": float(level[hap1 == 1].mean() - level[hap1 == 0].mean())}
        for name, cov in (("b", [b]), ("bl", [b, level])):
            d, se, p = hap_test(x, hap1, cov)
            row[f"delta_{name}"], row[f"se_{name}"], row[f"p_{name}"] = d, se, p
            row[f"pnull_{name}"] = hap_test(x, rng.permutation(hap1), cov)[2]
        if chrom in hets:
            pos = hets[chrom]
            i = np.searchsorted(pos, c)
            near = [abs(pos[j] - c) for j in (i - 1, i) if 0 <= j < len(pos)]
            row["nearest_het"] = min(near) if near else math.nan
            row["het_in_element"] = int(np.searchsorted(pos, c + h) - np.searchsorted(pos, c - h))
            row["het_1kb"] = int(np.searchsorted(pos, c + 1000) - np.searchsorted(pos, c - 1000))
            row["het_2_15kb"] = int(
                np.searchsorted(pos, c + 15000) - np.searchsorted(pos, c + 2000)
                + np.searchsorted(pos, c - 2000) - np.searchsorted(pos, c - 15000)
            )
            row["confident_1kb"] = covered(conf.get(chrom), c - 1000, c + 1000)
        rows.append(row)
    return rows


def main():
    out, chroms = sys.argv[1], sys.argv[2:]
    elements = [e for e in read_elements(CANDIDATES, "active") if e.chrom in chroms]
    by_chrom = defaultdict(dict)
    for e in elements:
        by_chrom[e.chrom][e.centre] = by_chrom[e.chrom].get(e.centre, False) or e.ctcf_bound
    hets = het_positions(set(chroms))
    conf = confident_intervals(set(chroms))
    rng = np.random.default_rng(0)
    rows = []
    for chrom in chroms:
        rows += element_rows(chrom, sorted(by_chrom[chrom]), by_chrom[chrom], hets, conf, rng)
        print(chrom, len(by_chrom[chrom]), "elements,", sum(r["chrom"] == chrom for r in rows),
              "tested", flush=True)
    table = pd.DataFrame(rows)
    for name in ("b", "bl"):
        for k in ("p", "pnull"):
            table[f"{k.replace('p', 'q', 1)}_{name}"] = np.nan
            for auto in (True, False):
                sel = (table.chrom != "chrX") == auto
                table.loc[sel, f"{k.replace('p', 'q', 1)}_{name}"] = benjamini_hochberg(
                    table.loc[sel, f"{k}_{name}"].to_numpy())
    table.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
