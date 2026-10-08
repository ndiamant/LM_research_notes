"""Summarize effects.tsv: effect-size distribution and candidate bin edges."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


CANDIDATES = {
    "profile": {  # percentage points
        "proposed 2/5/10": [2, 5, 10],
        "5/10/20": [5, 10, 20],
        "8/12/18": [8, 12, 18],
        "10/20": [10, 20],
        "10/15": [10, 15],
    },
    "covariation": {  # correlation units x 100
        "proposed 2/5/10": [2, 5, 10],
        "10/20/30": [10, 20, 30],
        "12/18/25": [12, 18, 25],
        "15/25": [15, 25],
        "20": [20],
    },
}


def _bins(values, edges):
    return np.digitize(values, edges)


def _labels(edges):
    bounds = [None, *edges, None]
    return [f"<{bounds[1]}" if lo is None else f">={lo}" if hi is None else f"{lo}-{hi}"
            for lo, hi in zip(bounds[:-1], bounds[1:])]


def main():
    parser = ArgumentParser()
    parser.add_argument("effects", type=Path)
    args = parser.parse_args()
    table = pd.read_csv(args.effects, sep="\t")
    regions = table.region_id.nunique()
    print(f"{regions} regions, {len(table)} region-pair rows")
    profiled = table.dropna(subset=["profile_effect_selection"])
    print(f"Regions with >=10 eligible sites: {profiled.region_id.nunique()} "
          f"({profiled.region_id.nunique() / regions:.0%})")
    cov = table.dropna(subset=["covariation_effect_selection"])
    print(f"Regions with >=10 eligible pairs: {cov.region_id.nunique()} "
          f"({cov.region_id.nunique() / regions:.0%})")
    check = profiled.dropna(subset=["profile_noise_selection_bootstrap"])
    ratio = check.profile_noise_selection / check.profile_noise_selection_bootstrap
    print(f"Profile noise, binomial/bootstrap: median {ratio.median():.3f} "
          f"(IQR {ratio.quantile(.25):.3f}-{ratio.quantile(.75):.3f}, n={len(check)})")

    for metric, candidates in CANDIDATES.items():
        frame = table.dropna(subset=[f"{metric}_effect_selection", f"{metric}_effect_evaluation"])
        sel = 100 * frame[f"{metric}_effect_selection"].to_numpy()
        ev = 100 * frame[f"{metric}_effect_evaluation"].to_numpy()
        unit = "pp" if metric == "profile" else "x100 correlation"
        print(f"\n=== {metric} ({unit}), n={len(frame)} region-pairs ===")
        quantiles = [.05, .1, .25, .5, .75, .9, .95]
        print("Quantiles " + "  ".join(f"{q:.0%}" for q in quantiles))
        print("selection  " + "  ".join(f"{v:5.1f}" for v in np.quantile(sel, quantiles)))
        print("evaluation " + "  ".join(f"{v:5.1f}" for v in np.quantile(ev, quantiles)))
        print(f"Spearman(selection, evaluation) = {spearmanr(sel, ev).statistic:.3f}")
        snr = frame[f"{metric}_energy_evaluation"] / frame[f"{metric}_noise_evaluation"]
        print(f"Evaluation energy/noise: median {snr.median():.2f}")
        if metric == "profile":
            off = frame.profile_offset_fraction_evaluation
            print(f"Offset fraction of evaluation energy: median {off.median():.2f}, "
                  f">0.5 in {(off > .5).mean():.0%}")
        per_pair = frame.assign(sel=sel).groupby(["cell_a", "cell_b"]).sel.median()
        print("Median selection effect by pair: "
              + ", ".join(f"{a}-{b} {v:.1f}" for (a, b), v in per_pair.items()))
        for name, edges in candidates.items():
            bs, be = _bins(sel, edges), _bins(ev, edges)
            print(f"\n  {name}: same bin in both sets {np.mean(bs == be):.0%}, "
                  f"within one bin {np.mean(abs(bs - be) <= 1):.0%}")
            for k, label in enumerate(_labels(edges)):
                mask = bs == k
                if not mask.any():
                    print(f"    {label:>8}: 0")
                    continue
                print(f"    {label:>8}: {mask.mean():5.1%} of rows; evaluation effect "
                      f"median {np.median(ev[mask]):5.1f} "
                      f"(IQR {np.quantile(ev[mask], .25):5.1f}-{np.quantile(ev[mask], .75):5.1f}); "
                      f"evaluation in same bin {np.mean(be[mask] == k):.0%}")


if __name__ == "__main__":
    main()
