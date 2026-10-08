"""Same-cell-type control: noise-corrected effect between two disjoint halves of one cell's molecules.

The true difference is zero, so the corrected effects show whether the noise
correction in effect_sizes.py removes sampling noise. Halves hold 35% of
molecules each, matching the depth behind each side of an evaluation-set
contrast at half the molecules.

    python analyses/differential_bins_20261006/control.py \
      --output "$SCRATCH/differential_bins_20261006/control"
"""

from argparse import ArgumentParser
from multiprocessing import Pool
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).parent))
import effect_sizes as es  # noqa: E402
from smf_net.evaluation.covariation import PairIndex, pair_correlations  # noqa: E402
from smf_net.evaluation.differential import _has_minor_calls  # noqa: E402


def control(task):
    region_id, seed = task
    generator = torch.Generator().manual_seed(seed)
    rows = []
    for c in es.CELLS:
        r = es._read_sample_region(es._handles[c]["regions"][f"{es._indexes[c][region_id]:08d}"],
                                   es._registries[c])
        reads = r.observed_smf.long()
        order = torch.randperm(len(reads), generator=generator)
        half = int(round(0.35 * len(reads)))
        a, b = reads[order[:half]], reads[order[half:2 * half]]
        keep = ((a != -1).sum(0) >= es.SITE_FLOOR["evaluation"] // 2) & (
            (b != -1).sum(0) >= es.SITE_FLOOR["evaluation"] // 2)
        if int(keep.sum()) < es.MIN_FEATURES:
            continue
        a, b = a[:, keep], b[:, keep]
        positions = r.observed_positions[keep]
        d = es._means(a.double())[0] - es._means(b.double())[0]
        effect, energy, noise = es._effect(d, es._binomial_noise(a) + es._binomial_noise(b))
        row = dict(region_id=region_id, cell=c, n_sites=int(keep.sum()),
                   profile_effect=effect, profile_energy=energy, profile_noise=noise)
        index = PairIndex.from_positions(positions)
        keep_pairs = index.separation <= es.MAX_SEPARATION
        corr = []
        for x in (a, b):
            pairs = pair_correlations(x, index)
            keep_pairs &= pairs.usable(es.JOINT_FLOOR["evaluation"] // 2)
            keep_pairs &= _has_minor_calls(x, index, es.MINOR_FLOOR["evaluation"])
            corr.append(pairs.correlation)
        if int(keep_pairs.sum()) >= es.MIN_FEATURES:
            sub = PairIndex(index.rows[keep_pairs], index.cols[keep_pairs],
                            index.separation[keep_pairs])
            stat = lambda y: torch.nan_to_num(pair_correlations(y, sub).correlation)  # noqa: E731
            noise = es._bootstrap(a, stat, 50, generator) + es._bootstrap(b, stat, 50, generator)
            d = corr[0][keep_pairs] - corr[1][keep_pairs]
            effect, energy, noise = es._effect(d, noise)
            row.update(n_pairs=int(keep_pairs.sum()), covariation_effect=effect,
                       covariation_energy=energy, covariation_noise=noise)
        rows.append(row)
    return rows


def main():
    parser = ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--regions", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    es._init()
    ids = sorted(set.intersection(*(set(i) for i in es._indexes.values())))
    chosen = sorted(np.random.default_rng(1).choice(ids, size=args.regions, replace=False).tolist())
    with Pool(args.workers, initializer=es._init) as pool:
        rows = [row for result in pool.imap(control, [(r, i) for i, r in enumerate(chosen)],
                                            chunksize=8) for row in result]
    table = pd.DataFrame(rows)
    table.to_csv(args.output / "control.tsv", sep="\t", index=False)
    for metric in ("profile", "covariation"):
        t = table.dropna(subset=[f"{metric}_effect"])
        effect = 100 * t[f"{metric}_effect"]
        ratio = t[f"{metric}_energy"] / t[f"{metric}_noise"]
        print(f"{metric}: n={len(t)}; corrected effect x100 quantiles "
              f"50% {effect.median():.2f}, 75% {effect.quantile(.75):.2f}, "
              f"90% {effect.quantile(.9):.2f}, 95% {effect.quantile(.95):.2f}; "
              f"zero in {(effect == 0).mean():.0%}; energy/noise median {ratio.median():.3f}")


if __name__ == "__main__":
    main()
