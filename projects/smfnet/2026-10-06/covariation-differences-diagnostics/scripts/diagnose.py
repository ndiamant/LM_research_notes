"""Why are predicted cell-type differences in pair correlations weak?

On random shared test regions, using the same molecule split and eligible
sites and pairs as scripts/differential_metrics.py, this compares predicted
and observed (evaluation-molecule) values three ways, for site means
(profile) and pair correlations (covariation):

1. absolute agreement per cell type, and the spread of predicted vs observed
   values;
2. the amplitude of cell-type differences: predicted / observed standard
   deviation of cell_a - cell_b, after removing sampling noise (1 = right size);
3. how similar two cell types are, in the predictions and in the data.

Every variance is corrected for sampling noise estimated by a whole-molecule
bootstrap (binomial for site means), and correlations are disattenuated, so
statistics measure the underlying values rather than read depth.

    python analyses/covariation_diagnostics_20261006/diagnose.py \
      --output "$SCRATCH/covariation_diagnostics_20261006"
"""

from argparse import ArgumentParser
from itertools import combinations
from multiprocessing import get_context
from pathlib import Path
import json
import os

import numpy as np
import pandas as pd
import torch

from smf_net.evaluation.differential import (
    DifferentialConfig,
    _bootstrap_noise,
    _bootstrap_weights,
    load_region_features,
    validate_sample_files,
    weighted_pair_correlations,
)

ROOT = Path(os.path.expandvars("$SCRATCH/evaluate_SMFNet_cell_lines"))
CELLS = ("MEL", "C2C12", "NP", "ES")
PATHS = {cell: str(ROOT / cell / "test_sample" / "samples.h5") for cell in CELLS}
SOURCES = {"observed": "evaluation", "predicted": "predicted"}


def region_values(region_id: str) -> dict | None:
    """Per-cell values and noise at one region's eligible sites and pairs."""
    torch.set_num_threads(1)
    config = DifferentialConfig()
    features = load_region_features(PATHS, [region_id], config)[region_id]
    keep = features.usable_pairs
    if len(features.positions) < config.min_correlation_sites or int(keep.sum()) < config.min_correlation_pairs:
        return None
    rows, cols = features.index.rows[keep], features.index.cols[keep]
    out = {"region_id": region_id}
    for name, source in SOURCES.items():
        out[f"profile_{name}"] = features.means[source].numpy()
        out[f"profile_{name}_noise"] = features.site_noise[source].numpy()
        out[f"covariation_{name}"] = features.pair_correlations[source][:, keep].numpy()
        noise = [_bootstrap_noise(weighted_pair_correlations(
            x, rows, cols, _bootstrap_weights(len(x), config.noise_resamples, features.generator)))[0]
            for x in features.reads[source]]
        out[f"covariation_{name}_noise"] = torch.stack(noise).numpy()
    return out


def disattenuated(x, y, vx, vy) -> dict:
    """Noise-corrected correlation and standard deviations of paired values."""
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(vx) & np.isfinite(vy)
    x, y, vx, vy = x[keep], y[keep], vx[keep], vy[keep]
    var_x = max(np.var(x) - vx.mean(), 1e-12)
    var_y = max(np.var(y) - vy.mean(), 1e-12)
    covariance = np.mean((x - x.mean()) * (y - y.mean()))
    return dict(n=int(keep.sum()), raw_r=float(np.corrcoef(x, y)[0, 1]),
                true_r=float(covariance / np.sqrt(var_x * var_y)),
                sd_x=float(np.sqrt(var_x)), sd_y=float(np.sqrt(var_y)),
                mean_x=float(x.mean()), mean_y=float(y.mean()))


def main():
    parser = ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--regions", type=int, default=1500)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    import h5py
    handles = {c: h5py.File(p, "r") for c, p in PATHS.items()}
    _, indexes = validate_sample_files(handles)
    ids = sorted(set.intersection(*(set(i) for i in indexes.values())))
    for handle in handles.values():
        handle.close()
    chosen = sorted(np.random.default_rng(7).choice(ids, args.regions, replace=False).tolist())
    with get_context("spawn").Pool(args.workers) as pool:
        regions = [r for r in pool.map(region_values, chosen, chunksize=8) if r is not None]

    def stack(key):
        return np.concatenate([r[key] for r in regions], axis=1)

    rows = []
    for metric in ("profile", "covariation"):
        obs, pred = stack(f"{metric}_observed"), stack(f"{metric}_predicted")
        v_obs, v_pred = stack(f"{metric}_observed_noise"), stack(f"{metric}_predicted_noise")
        for i, cell in enumerate(CELLS):
            stats = disattenuated(pred[i], obs[i], v_pred[i], v_obs[i])
            rows.append(dict(metric=metric, analysis="absolute", comparison=cell, **stats,
                             amplitude=stats["sd_x"] / stats["sd_y"]))
        for i, j in combinations(range(len(CELLS)), 2):
            name = f"{CELLS[i]}-{CELLS[j]}"
            stats = disattenuated(pred[i] - pred[j], obs[i] - obs[j],
                                  v_pred[i] + v_pred[j], v_obs[i] + v_obs[j])
            rows.append(dict(metric=metric, analysis="difference", comparison=name, **stats,
                             amplitude=stats["sd_x"] / stats["sd_y"]))
            for source, values, noise in (("observed", obs, v_obs), ("predicted", pred, v_pred)):
                stats = disattenuated(values[i], values[j], noise[i], noise[j])
                rows.append(dict(metric=metric, analysis=f"similarity_{source}", comparison=name,
                                 **stats, amplitude=float("nan")))
    table = pd.DataFrame(rows)
    table.to_csv(args.output / "diagnostics.tsv", sep="\t", index=False)
    (args.output / "metadata.json").write_text(json.dumps(dict(
        paths=PATHS, shared_regions=len(ids), sampled=len(chosen), analyzed=len(regions),
        seed=7, config="DifferentialConfig() defaults",
    ), indent=2) + "\n")
    pd.set_option("display.width", 200)
    print(f"{len(regions)} of {len(chosen)} sampled regions have enough sites and pairs")
    for metric in ("profile", "covariation"):
        frame = table[table.metric == metric]
        print(f"\n=== {metric}")
        print("Absolute values per cell type (x = predicted, y = observed):")
        print(frame[frame.analysis == "absolute"][
            ["comparison", "n", "raw_r", "true_r", "mean_x", "mean_y", "sd_x", "sd_y", "amplitude"]
        ].round(3).to_string(index=False))
        print("cell_a - cell_b differences:")
        print(frame[frame.analysis == "difference"][
            ["comparison", "raw_r", "true_r", "sd_x", "sd_y", "amplitude"]].round(3).to_string(index=False))
        similarity = frame[frame.analysis.str.startswith("similarity")].pivot(
            index="comparison", columns="analysis", values="true_r")
        print("Similarity between cell types (noise-corrected r):")
        print(similarity.round(3).to_string())


if __name__ == "__main__":
    main()
