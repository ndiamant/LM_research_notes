"""Empirical distribution of noise-corrected cell-type differences on the test split.

For each shared region, each cell type's observed molecules are split at random
into a selection set (default 30%) and an evaluation set. For every cell-type pair,
the noise-corrected RMS difference is estimated separately in each set:

    effect = sqrt(max(0, mean(d**2) - mean(noise variance of d)))

where d is the per-site difference in mean accessibility (profile) or the
per-pair difference in within-molecule correlation (covariation). Profile noise
is binomial, p(1-p)/(n-1) per site and cell, checked against a whole-molecule
bootstrap on a subset of regions. Covariation noise comes from a whole-molecule
bootstrap. The two sets share no molecules, so their agreement shows how well
a selection set of that size places a region in the right bin. Floors are the
differential defaults (30 calls per site, 50 joint reads, 5 minor calls)
scaled by each set's share of molecules.

Run on a compute node with load_smf active:

    python analyses/differential_bins_20261006/effect_sizes.py \
      --output "$SCRATCH/differential_bins_20261006"
"""

from argparse import ArgumentParser
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
import json
import os

import h5py
import numpy as np
import pandas as pd
import torch

from smf_net.evaluation.covariation import PairIndex, pair_correlations
from smf_net.evaluation.differential import _has_minor_calls, validate_sample_files
from smf_net.evaluation.sample_io import _read_sample_region, _read_source_registry


ROOT = Path(os.path.expandvars("$SCRATCH/evaluate_SMFNet_cell_lines"))
CELLS = ("MEL", "C2C12", "NP", "ES")
PATHS = {cell: ROOT / cell / "test_sample" / "samples.h5" for cell in CELLS}
SELECTION_FRACTION = 0.3
SITE_FLOOR = {"selection": 9, "evaluation": 21}
JOINT_FLOOR = {"selection": 15, "evaluation": 35}
MINOR_FLOOR = {"selection": 2, "evaluation": 4}
MAX_SEPARATION = 500
MIN_FEATURES = 10

_handles = _indexes = _registries = None


def _set_selection_fraction(fraction):
    """Rescale the floors to each set's share of molecules (workers fork after this)."""
    global SELECTION_FRACTION
    SELECTION_FRACTION = fraction
    for floors, full in ((SITE_FLOOR, 30), (JOINT_FLOOR, 50), (MINOR_FLOOR, 5)):
        floors["selection"] = max(2, round(full * fraction))
        floors["evaluation"] = max(2, round(full * (1 - fraction)))


def _init():
    global _handles, _indexes, _registries
    torch.set_num_threads(1)
    _handles = {c: h5py.File(p, "r") for c, p in PATHS.items()}
    _, _indexes = validate_sample_files(_handles)
    _registries = {c: _read_source_registry(_handles[c]) for c in CELLS}


def _means(reads):
    valid = reads != -1
    n = valid.sum(0)
    return reads.clamp_min(0).sum(0) / n.clamp_min(1), n


def _binomial_noise(reads):
    mean, n = _means(reads.double())
    return mean * (1 - mean) / (n - 1).clamp_min(1)


def _bootstrap(reads, statistic, resamples, generator):
    draws = [statistic(reads[torch.randint(len(reads), (len(reads),), generator=generator)])
             for _ in range(resamples)]
    return torch.stack(draws).var(dim=0)


def _effect(d, noise):
    energy = float(d.square().mean())
    noise = float(noise.mean())
    return np.sqrt(max(0.0, energy - noise)), energy, noise


def score(task):
    region_id, seed, bootstrap_profile, cov_resamples = task
    generator = torch.Generator().manual_seed(seed)
    regions = {c: _read_sample_region(_handles[c]["regions"][f"{_indexes[c][region_id]:08d}"],
                                      _registries[c]) for c in CELLS}
    positions = regions[CELLS[0]].observed_positions.numpy()
    for r in regions.values():
        positions = np.intersect1d(positions, r.observed_positions.numpy())
        positions = np.intersect1d(positions, r.predicted_positions.numpy())
    sets = {"selection": {}, "evaluation": {}}
    for c, r in regions.items():
        reads = r.observed_smf[:, np.searchsorted(r.observed_positions.numpy(), positions)].long()
        order = torch.randperm(len(reads), generator=generator)
        cut = int(round(SELECTION_FRACTION * len(reads)))
        sets["selection"][c] = reads[order[:cut]]
        sets["evaluation"][c] = reads[order[cut:]]
    keep = np.ones(len(positions), dtype=bool)
    for name, by_cell in sets.items():
        for reads in by_cell.values():
            keep &= ((reads != -1).sum(0) >= SITE_FLOOR[name]).numpy()
    sites = torch.as_tensor(np.nonzero(keep)[0])
    positions = positions[keep]
    rows = []
    base = dict(region_id=region_id, n_sites=len(positions),
                n_reads=int(min(len(r.observed_smf) for r in regions.values())))
    if len(positions) < MIN_FEATURES:
        return [dict(base, cell_a=a, cell_b=b) for a, b in combinations(CELLS, 2)]
    sets = {name: {c: x[:, sites] for c, x in by_cell.items()} for name, by_cell in sets.items()}

    # Profile.
    profile = {}
    for name, by_cell in sets.items():
        means = {c: _means(x.double())[0] for c, x in by_cell.items()}
        noise = {c: _binomial_noise(x) for c, x in by_cell.items()}
        boot = ({c: _bootstrap(x.double(), lambda y: _means(y)[0], 100, generator)
                 for c, x in by_cell.items()}
                if bootstrap_profile and name == "selection" else None)
        profile[name] = (means, noise, boot)

    # Covariation: pairs eligible in both sets and every cell type.
    index = PairIndex.from_positions(torch.from_numpy(positions))
    keep_pairs = index.separation <= MAX_SEPARATION
    correlations = {"selection": {}, "evaluation": {}}
    for name, by_cell in sets.items():
        for c, x in by_cell.items():
            pairs = pair_correlations(x, index)
            keep_pairs &= pairs.usable(JOINT_FLOOR[name])
            keep_pairs &= _has_minor_calls(x, index, MINOR_FLOOR[name])
            correlations[name][c] = pairs.correlation
    n_pairs = int(keep_pairs.sum())
    cov_noise = {"selection": {}, "evaluation": {}}
    if n_pairs >= MIN_FEATURES and cov_resamples:
        sub = PairIndex(index.rows[keep_pairs], index.cols[keep_pairs],
                        index.separation[keep_pairs])
        for name, by_cell in sets.items():
            for c, x in by_cell.items():
                # nan_to_num: a resample can make a site constant; that draw
                # contributes zero rather than poisoning the variance.
                cov_noise[name][c] = _bootstrap(
                    x, lambda y: torch.nan_to_num(pair_correlations(y, sub).correlation),
                    cov_resamples, generator)

    for a, b in combinations(CELLS, 2):
        row = dict(base, cell_a=a, cell_b=b, n_pairs=n_pairs)
        for name in sets:
            means, noise, boot = profile[name]
            d = means[a] - means[b]
            effect, energy, noise_energy = _effect(d, noise[a] + noise[b])
            row[f"profile_effect_{name}"] = effect
            row[f"profile_energy_{name}"] = energy
            row[f"profile_noise_{name}"] = noise_energy
            row[f"profile_offset_fraction_{name}"] = float(d.mean() ** 2 / d.square().mean()) if energy > 0 else np.nan
            if boot is not None:
                row["profile_noise_selection_bootstrap"] = float((boot[a] + boot[b]).mean())
            if cov_noise[name]:
                d = correlations[name][a][keep_pairs] - correlations[name][b][keep_pairs]
                effect, energy, noise_energy = _effect(d, cov_noise[name][a] + cov_noise[name][b])
                row[f"covariation_effect_{name}"] = effect
                row[f"covariation_energy_{name}"] = energy
                row[f"covariation_noise_{name}"] = noise_energy
        rows.append(row)
    return rows


def main():
    parser = ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--regions", type=int, default=3000)
    parser.add_argument("--bootstrap-check", type=int, default=300)
    parser.add_argument("--cov-resamples", type=int, default=50)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--selection-fraction", type=float, default=SELECTION_FRACTION)
    args = parser.parse_args()
    _set_selection_fraction(args.selection_fraction)
    args.output.mkdir(parents=True, exist_ok=True)
    _init()
    ids = sorted(set.intersection(*(set(i) for i in _indexes.values())))
    rng = np.random.default_rng(0)
    chosen = sorted(rng.choice(ids, size=min(args.regions, len(ids)), replace=False).tolist())
    tasks = [(r, i, i < args.bootstrap_check, args.cov_resamples) for i, r in enumerate(chosen)]
    with Pool(args.workers, initializer=_init) as pool:
        rows = [row for result in pool.imap(score, tasks, chunksize=8) for row in result]
    table = pd.DataFrame(rows)
    table.to_csv(args.output / "effects.tsv", sep="\t", index=False)
    (args.output / "metadata.json").write_text(json.dumps(dict(
        paths={c: str(p) for c, p in PATHS.items()}, shared_regions=len(ids),
        analyzed_regions=len(chosen), selection_fraction=SELECTION_FRACTION,
        site_floor=SITE_FLOOR, joint_floor=JOINT_FLOOR, minor_floor=MINOR_FLOOR,
        max_separation=MAX_SEPARATION, min_features=MIN_FEATURES,
        cov_resamples=args.cov_resamples, bootstrap_check=args.bootstrap_check,
    ), indent=2) + "\n")
    print(f"Wrote {len(table)} rows for {len(chosen)} of {len(ids)} shared regions")


if __name__ == "__main__":
    main()
