"""Does sampling only a short window fix out-of-distribution dense reads?

GCH training reads span about 165 bp of the 1,024-bp output window, but
default sampling generates every site of the window in each read. This
compares two sample files from the same model, regions and settings: a
default ("dense") run and one sampled with ``sampling.sampling_window``
("windowed"). Both are scored on the same window sites and pairs against the
same observed reads.

Statistics pool sites (profile: per-site mean accessibility) or site pairs
(covariation: within-molecule pair correlation) across regions:

- ``pooled``: all features, so between-region differences count;
- ``within``: features relative to their region's mean.

Each reports Pearson r, its noise ceiling (the r a perfect model reaches given
sampling noise in observed and predicted reads), r / ceiling, a
noise-corrected slope (profile), and amplitude, the predicted / observed
noise-corrected standard deviation (1 = right spread). Noise is binomial for
site means and a whole-molecule bootstrap for pair correlations and region
means. Differences (windowed - dense) get a paired region-bootstrap 95% CI.

    python analyses/windowed_sampling_20261006/compare_windows.py \\
      --dense "$OUT/dense/samples.h5" --windowed "$OUT/windowed/samples.h5" \\
      --window 384 640 --output "$OUT/comparison"
"""

from argparse import ArgumentParser
from pathlib import Path
import json
import zlib

import h5py
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

from smf_net.evaluation.covariation import PairIndex, pair_correlations
from smf_net.evaluation.differential import (
    _binomial_variance,
    _bootstrap_noise,
    _bootstrap_weights,
    _has_minor_calls,
    _means,
    _pooled,
    _weighted_site_means,
    weighted_pair_correlations,
)
from smf_net.evaluation.sample_io import _read_sample_region, _read_source_registry

RUNS = ("dense", "windowed")
MATCHED_ATTRS = ("checkpoint", "split", "reads_per_region", "sample_steps", "seed",
                 "model_input_length", "assay_type", "data_config")


def _check_files(handles: dict[str, h5py.File], window: tuple[int, int]) -> list[str]:
    for key in MATCHED_ATTRS:
        values = {run: handles[run].attrs.get(key) for run in RUNS}
        if values["dense"] != values["windowed"]:
            raise ValueError(f"{key} differs between runs: {values}")
    if "sampling_window_start" in handles["dense"].attrs:
        raise ValueError("The dense file was sampled with a window")
    stored = tuple(int(handles["windowed"].attrs.get(f"sampling_window_{end}", -1))
                   for end in ("start", "end"))
    if stored != window:
        raise ValueError(f"Windowed file was sampled in {stored}, not {window}")
    ids = {run: handles[run]["region_ids"].asstr()[:].tolist() for run in RUNS}
    if ids["dense"] != ids["windowed"]:
        raise ValueError("The runs hold different regions or region order")
    return ids["dense"]


def _pool_row(x, y, vx, vy) -> np.ndarray:
    return np.array([len(x), x.sum(), y.sum(), (x * x).sum(), (y * y).sum(), (x * y).sum(),
                     vx.sum(), vy.sum()], dtype=float)


def _within_row(x, y, vx, vy, mean_vx, mean_vy) -> np.ndarray:
    n = len(x)
    return np.array([n, 0.0, 0.0, (x * x).sum() - x.sum() ** 2 / n,
                     (y * y).sum() - y.sum() ** 2 / n, (x * y).sum() - x.sum() * y.sum() / n,
                     vx.sum() - n * mean_vx, vy.sum() - n * mean_vy], dtype=float)


def _summary(stats: np.ndarray) -> dict[str, float]:
    point = {k: float(v) for k, v in _pooled(stats).items()}
    n, sx, sy, sxx, syy, _, vx, vy = stats
    with np.errstate(invalid="ignore", divide="ignore"):
        var_x, var_y = sxx / n - (sx / n) ** 2, syy / n - (sy / n) ** 2
        point["amplitude"] = float(np.sqrt((var_x - vx / n) / (var_y - vy / n)))
    point["mean_predicted"], point["mean_observed"] = float(sx / n), float(sy / n)
    point.pop("evaluation_effect")
    return point


def main():
    parser = ArgumentParser()
    parser.add_argument("--dense", type=Path, required=True)
    parser.add_argument("--windowed", type=Path, required=True)
    parser.add_argument("--window", type=int, nargs=2, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-site-coverage", type=int, default=30)
    parser.add_argument("--min-predicted-site-coverage", type=int, default=10)
    parser.add_argument("--min-joint-reads", type=int, default=50)
    parser.add_argument("--min-minor-calls", type=int, default=5)
    parser.add_argument("--min-sites", type=int, default=5)
    parser.add_argument("--noise-resamples", type=int, default=50)
    parser.add_argument("--ci-resamples", type=int, default=1000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    window = tuple(args.window)
    torch.set_num_threads(1)

    handles = {"dense": h5py.File(args.dense, "r"), "windowed": h5py.File(args.windowed, "r")}
    ids = _check_files(handles, window)
    registries = {run: _read_source_registry(h) for run, h in handles.items()}
    rows = {(metric, scope, run): [] for metric in ("profile", "covariation")
            for scope in ("pooled", "within") for run in RUNS}
    counts = []
    for index, region_id in enumerate(tqdm(ids, desc="Scoring regions", unit="region")):
        regions = {run: _read_sample_region(handles[run]["regions"][f"{index:08d}"], registries[run])
                   for run in RUNS}
        observed = regions["dense"].observed_smf
        if not (torch.equal(observed, regions["windowed"].observed_smf) and torch.equal(
                regions["dense"].observed_positions, regions["windowed"].observed_positions)):
            raise ValueError(f"{region_id}: observed reads differ between runs")
        generator = torch.Generator().manual_seed(zlib.crc32(region_id.encode()))
        positions = regions["dense"].observed_positions.numpy()
        positions = positions[(positions >= window[0]) & (positions < window[1])]
        for run in RUNS:
            positions = np.intersect1d(positions, regions[run].predicted_positions.numpy())
        reads = {"observed": observed[:, np.searchsorted(
            regions["dense"].observed_positions.numpy(), positions)].long()}
        for run in RUNS:
            reads[run] = regions[run].predicted_smf[:, np.searchsorted(
                regions[run].predicted_positions.numpy(), positions)].long()
        keep = ((reads["observed"] != -1).sum(0) >= args.min_site_coverage)
        for run in RUNS:
            keep &= (reads[run] != -1).sum(0) >= args.min_predicted_site_coverage
        reads = {k: v[:, keep] for k, v in reads.items()}
        positions = positions[keep.numpy()]
        if len(positions) < args.min_sites:
            continue

        weights = {k: _bootstrap_weights(len(v), args.noise_resamples, generator)
                   for k, v in reads.items()}
        means = {k: _means(v).double().numpy() for k, v in reads.items()}
        site_noise = {k: _binomial_variance(v).numpy() for k, v in reads.items()}
        mean_noise = {k: float(_bootstrap_noise(_weighted_site_means(v, weights[k]))[1])
                      for k, v in reads.items()}
        for run in RUNS:
            x, y, vx, vy = means[run], means["observed"], site_noise[run], site_noise["observed"]
            rows[("profile", "pooled", run)].append((region_id, _pool_row(x, y, vx, vy)))
            rows[("profile", "within", run)].append((region_id, _within_row(
                x, y, vx, vy, mean_noise[run], mean_noise["observed"])))

        index_pairs = PairIndex.from_positions(torch.from_numpy(positions))
        usable = torch.ones(index_pairs.num_pairs, dtype=torch.bool)
        correlations = {}
        for k, v in reads.items():
            pairs = pair_correlations(v, index_pairs)
            usable &= pairs.usable(args.min_joint_reads)
            usable &= _has_minor_calls(v, index_pairs, args.min_minor_calls)
            correlations[k] = pairs.correlation
        n_pairs = int(usable.sum())
        counts.append(dict(region_id=region_id, sites=len(positions), pairs=n_pairs))
        if n_pairs < args.min_sites:
            continue
        rows_i, cols_i = index_pairs.rows[usable], index_pairs.cols[usable]
        pair_noise, pair_mean_noise = {}, {}
        for k, v in reads.items():
            feature, mean = _bootstrap_noise(weighted_pair_correlations(v, rows_i, cols_i, weights[k]))
            pair_noise[k], pair_mean_noise[k] = feature.numpy(), float(mean)
        y = correlations["observed"][usable].double().numpy()
        for run in RUNS:
            x = correlations[run][usable].double().numpy()
            rows[("covariation", "pooled", run)].append((region_id, _pool_row(
                x, y, pair_noise[run], pair_noise["observed"])))
            rows[("covariation", "within", run)].append((region_id, _within_row(
                x, y, pair_noise[run], pair_noise["observed"],
                pair_mean_noise[run], pair_mean_noise["observed"])))

    rng = np.random.default_rng(0)
    results = []
    for metric in ("profile", "covariation"):
        for scope in ("pooled", "within"):
            stats = {run: np.array([r for _, r in rows[(metric, scope, run)]]) for run in RUNS}
            n_regions = len(stats["dense"])
            points = {run: _summary(stats[run].sum(0)) for run in RUNS}
            for run in RUNS:
                results.append(dict(metric=metric, scope=scope, run=run, regions=n_regions,
                                    features=int(stats[run][:, 0].sum()), **points[run]))
            draws = rng.integers(n_regions, size=(args.ci_resamples, n_regions))
            difference = {"pearson": [], "fraction_of_ceiling": [], "amplitude": []}
            for draw in draws:
                boot = {run: _summary(stats[run][draw].sum(0)) for run in RUNS}
                for key in difference:
                    difference[key].append(boot["windowed"][key] - boot["dense"][key])
            row = dict(metric=metric, scope=scope, run="windowed - dense", regions=n_regions,
                       features=int(stats["dense"][:, 0].sum()))
            for key, values in difference.items():
                values = np.array(values)
                values = values[np.isfinite(values)]
                row[key] = points["windowed"][key] - points["dense"][key]
                row[f"{key}_low"], row[f"{key}_high"] = np.quantile(values, [.025, .975])
            results.append(row)
    table = pd.DataFrame(results)
    table.to_csv(args.output / "comparison.tsv", sep="\t", index=False)
    pd.DataFrame(counts).to_csv(args.output / "region_features.tsv", sep="\t", index=False)
    (args.output / "metadata.json").write_text(json.dumps(dict(
        dense=str(args.dense.resolve()), windowed=str(args.windowed.resolve()), window=window,
        regions_in_files=len(ids), checkpoint=handles["dense"].attrs["checkpoint"],
        floors=dict(site=args.min_site_coverage, predicted_site=args.min_predicted_site_coverage,
                    joint=args.min_joint_reads, minor=args.min_minor_calls,
                    min_features=args.min_sites),
        noise_resamples=args.noise_resamples, ci_resamples=args.ci_resamples,
    ), indent=2) + "\n")
    pd.set_option("display.width", 220)
    columns = ["metric", "scope", "run", "regions", "features", "pearson", "ceiling",
               "fraction_of_ceiling", "amplitude", "slope", "mean_predicted", "mean_observed"]
    print(table[[c for c in columns if c in table]].round(3).to_string(index=False))
    diffs = table[table.run == "windowed - dense"]
    print("\nwindowed - dense, 95% paired region-bootstrap CIs:")
    for row in diffs.itertuples():
        print(f"  {row.metric:11s} {row.scope:6s} r {row.pearson:+.3f} [{row.pearson_low:+.3f}, "
              f"{row.pearson_high:+.3f}]  fraction {row.fraction_of_ceiling:+.3f} "
              f"[{row.fraction_of_ceiling_low:+.3f}, {row.fraction_of_ceiling_high:+.3f}]  "
              f"amplitude {row.amplitude:+.3f} [{row.amplitude_low:+.3f}, {row.amplitude_high:+.3f}]")


if __name__ == "__main__":
    main()
