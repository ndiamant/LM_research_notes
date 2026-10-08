"""Is one cell-type model much worse than the others?

For each cell type c, compares the predicted and observed deviation of c from
the mean of the other cell types, on the same random test regions, split and
eligible features as diagnose.py. Reports noise-corrected (disattenuated)
correlation and amplitude, pooled and within regions, plus each cell type's
observed read depth.

Within-region values are centered on their region's mean; their noise is
taken as the feature noise times (1 - 1/n), which treats a region's features
as independent and so slightly overstates within-region noise for
covariation.

    python analyses/covariation_diagnostics_20261006/one_vs_rest.py \
      --output "$SCRATCH/covariation_diagnostics_20261006"
"""

from argparse import ArgumentParser
from multiprocessing import get_context
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from diagnose import CELLS, PATHS, disattenuated, region_values  # noqa: E402
from smf_net.evaluation.differential import validate_sample_files  # noqa: E402


def one_vs_rest(values: np.ndarray, noise: np.ndarray, cell: int):
    others = [i for i in range(len(values)) if i != cell]
    deviation = values[cell] - values[others].mean(0)
    variance = noise[cell] + noise[others].sum(0) / len(others) ** 2
    return deviation, variance


def main():
    parser = ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--regions", type=int, default=1500)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    handles = {c: h5py.File(p, "r") for c, p in PATHS.items()}
    _, indexes = validate_sample_files(handles)
    ids = sorted(set.intersection(*(set(i) for i in indexes.values())))
    chosen = sorted(np.random.default_rng(7).choice(ids, args.regions, replace=False).tolist())
    depth = []
    for cell, handle in handles.items():
        for region_id in chosen[:500]:
            smf = handle["regions"][f"{indexes[cell][region_id]:08d}"]["observed/smf"][:]
            depth.append(dict(cell=cell, molecules=len(smf),
                              calls_per_site=float(np.median((smf != -1).sum(0)))))
        handle.close()
    with get_context("spawn").Pool(args.workers) as pool:
        regions = [r for r in pool.map(region_values, chosen, chunksize=8) if r is not None]

    rows = []
    for metric in ("profile", "covariation"):
        for i, cell in enumerate(CELLS):
            parts = {"pooled": [], "within": []}
            for region in regions:
                d_obs, v_obs = one_vs_rest(region[f"{metric}_observed"],
                                           region[f"{metric}_observed_noise"], i)
                d_pred, v_pred = one_vs_rest(region[f"{metric}_predicted"],
                                             region[f"{metric}_predicted_noise"], i)
                parts["pooled"].append((d_pred, d_obs, v_pred, v_obs))
                shrink = 1 - 1 / len(d_obs)
                parts["within"].append((d_pred - d_pred.mean(), d_obs - d_obs.mean(),
                                        v_pred * shrink, v_obs * shrink))
            for scope, chunks in parts.items():
                x, y, vx, vy = (np.concatenate(c) for c in zip(*chunks))
                if scope == "pooled":
                    x, y = x - np.nanmean(x), y - np.nanmean(y)
                stats = disattenuated(x, y, vx, vy)
                rows.append(dict(metric=metric, cell=cell, scope=scope, true_r=stats["true_r"],
                                 amplitude=stats["sd_x"] / stats["sd_y"],
                                 observed_sd=stats["sd_y"]))
    table = pd.DataFrame(rows)
    table.to_csv(args.output / "one_vs_rest.tsv", sep="\t", index=False)
    depth = pd.DataFrame(depth).groupby("cell", sort=False).median()
    depth.to_csv(args.output / "observed_depth.tsv", sep="\t")
    pd.set_option("display.width", 200)
    print(f"{len(regions)} regions")
    print("Observed depth per test region (median over 500 regions):")
    print(depth.round(1).to_string())
    print("\nEach cell type vs the mean of the others (noise-corrected):")
    print(table.pivot_table(index=["metric", "cell"], columns="scope",
                            values=["true_r", "amplitude"], sort=False).round(2).to_string())


if __name__ == "__main__":
    main()
