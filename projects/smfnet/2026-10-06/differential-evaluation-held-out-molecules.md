---
title: Differential cell-type evaluation with held-out molecules, effect-size strata, and between/within-region scores
date: 2026-10-06
project: smfnet
agent: Claude Code
status: draft
sources:
  - SMFNet repo, scripts/differential_metrics.py and scripts/differential_metrics.md (branch nd/more_assays at 1151cbb, uncommitted at time of writing)
  - SMFNet repo, src/smf_net/evaluation/differential.py and differential_plotting.py (same branch, uncommitted)
  - SMFNet repo, tests/test_evaluation_differential.py (same branch, uncommitted)
  - SMFNet repo, analyses/differential_bins_20261006/ (effect_sizes.py, control.py, summarize.py)
  - /scratch/users/diamant/evaluate_SMFNet_cell_lines/differential_test_cell_lines_full/ (full test-split run)
  - /scratch/users/diamant/differential_bins_20261006/ (bin calibration runs)
  - differential-evaluation-test-split/ in this folder (copied tables and figures)
  - User discussion, 2026-10-06
tags: [smfnet, differential-effects, cell-types, gch, evaluation, noise-ceiling, held-out-molecules, bootstrap, strata, batch-effects]
---

# Summary

`scripts/differential_metrics.py` now tests whether the NP, MEL, ES and C2C12
models predict **how cell types differ**, scored against observed molecules
that were not used to pick the regions. Each region's observed molecules are
split 30/70 into a selection set, which only places the region in an
effect-size bin, and an evaluation set, which the model is scored against.
Results are Pearson r between predicted and observed `cell_a − cell_b`
differences beside a noise ceiling (the r a perfect model reaches given
sampling noise), split into a **between-region** part (each region's mean
difference) and a **within-region** part (which sites inside a region
differ). This follows the recommendation in
[the chr16 feasibility note](../2026-10-05/chr16-differential-effects-feasibility.md)
to select and score on independent molecules.

On the full test split (25,860 regions, 4 models), predictions improve
steadily with the size of the observed difference but stay well below what
the data can resolve. In the strongest profile bin (≥20 pp), the between-region
r is 0.52 (ceiling 0.96) and the within-region r is 0.35 (ceiling 0.89).
Covariation differences are weaker, explained in
[the companion diagnosis note](covariation-differences-diagnostics.md).

# Key Points

- **Design:** select on 30% of molecules, score on the other 70%, so strong
  regions are not chosen for noise that the score then rewards (winner's curse).
  The split is seeded by `seed` and region ID, so it doesn't depend on which
  regions are analyzed.
- **Bins (empirically chosen):** profile <10, 10–20, ≥20 pp; covariation <0.2,
  ≥0.2. The first-proposed 2/5/10 pp edges could not be told apart at 30% of
  the molecules, and almost no region differs by under 5 pp.
- **Ceiling:** √(reliability_observed × reliability_predicted), from noise
  estimates that a simulated perfect model validates (fraction of ceiling within
  0.01–0.02 of 1 for all six profile/covariation × pooled/between/within scores).
- **Between vs within:** pooled r mixes "which regions shift" with "which sites
  differ". The two behave differently, so report both.
- **Contrast centering:** `contrast=all` rows subtract each contrast's own mean
  difference before pooling. Without this, pooling rewards six cell-type-wide
  offsets (e.g. C2C12 is ~12 pp more accessible than NP overall), which looked
  like regional agreement and are also what batch effects would produce. The
  offsets are reported separately (`contrast_offset_*`).
- **Main result:** the model reaches 17–44% of the pooled ceiling (profile, by
  bin) and 23–35% (covariation). The within-region scores, the most specific test
  of cell-type-specific prediction, reach at most 39% of their ceiling.

# Details

## Method

All of this is in `src/smf_net/evaluation/differential.py` (scoring, strata)
and `differential_plotting.py` (figures), documented in
`scripts/differential_metrics.md`.

**Eligibility.** Observed floors are set for all of a region's molecules and
scaled to each set's share, rounding up. Predicted floors are not scaled. A site
or pair must pass in every cell type and in the selection, evaluation and
predicted reads.

| Floor | Set as | Selection | Evaluation | Predicted |
|---|---|---|---|---|
| Calls per site | 30 | 9 | 21 | 10 |
| Molecules covering both sites of a pair | 50 | 15 | 35 | 50 |
| Calls of each state at both sites, within those | 5 | 2 | 4 | 5 |
| Max pair separation | 500 bp | | | |
| Sites / pairs per region | 10 / 10 | | | |

These floors are stricter than `metrics.py`'s (10) because the between-cell
signal left after removing the shared mean is small next to sampling noise, and
noise adds directly to RMSE and explained variance. The minor-call floor is now
applied to predicted reads as well as observed, after we found the old filter
was asymmetric. 84% of test regions have ≥10 eligible sites and 81% have ≥10
eligible pairs.

**Effect size (bins).** For a region and contrast, the selection-set effect is
the noise-corrected RMS of the per-site difference,
`sqrt(max(0, mean(d²) − mean(noise)))`. Profile effects are in accessibility
units and reported as **percentage points (pp)**: 40% vs 60% accessible is
20 pp. Covariation effects are in units of within-molecule correlation.

**Noise.** Per-site noise for profiles is binomial, p(1 − p)/(n − 1). It
matched a whole-molecule bootstrap (median ratio 1.004). Pair-correlation noise,
and the noise of each region's mean difference for both metrics, come from a
whole-molecule bootstrap (50 resamples, vectorized as multinomial weights
times matrix products). Resampling whole molecules keeps the correlation
between sites on the same molecule. For that reason a region mean's noise is
not the sum of its sites' noise.

**Strata statistics.** Within each contrast × bin, features are pooled
across regions using per-row sums stored in `pairwise.tsv` (`pool_*`,
`mean_noise_*`), so `strata.tsv` can be recomputed without rescoring:

- `pearson`, `ceiling`, `fraction_of_ceiling`: all pooled features.
- `between_*`: one point per region and contrast (its mean difference).
- `within_*`: features after subtracting their region's mean. The expected
  noise of a centered value is the summed feature noise minus n times the
  noise of the mean.
- `slope` (profile only): observed on predicted, fit with an intercept and
  corrected for the predictions' noise. Covariation reports no slope, as in
  `metrics.py`'s correlation of correlations
  ([background](../2026-09-29/positional-covariation-baselines-and-noise-matched-corr-of-corr.md)).
- 95% CIs resample whole regions (1,000 draws).

**Contrast centering.** In `contrast=all` rows, each contrast is centered on
its own mean within each bin before pooling. The centering means come from the
whole stratum and stay fixed inside the bootstrap, which is negligible with
thousands of regions per contrast. A test confirms that adding a constant to one
contrast leaves every `contrast=all` statistic and that contrast's
correlations and slopes unchanged, and moves only its offset columns.

## Bin calibration

`analyses/differential_bins_20261006/` estimated selection- and
evaluation-set effects independently on 3,000 random test regions (18,000
region × contrast rows). Summaries are in
`differential-evaluation-test-split/bin_calibration_selection30.txt` and
`_selection50.txt`.

- **Same-cell-type control** (one cell type's molecules split in two): median
  corrected effect 0, energy/noise 0.92 (profile) and 0.99 (covariation). The
  correction is calibrated. But a truly zero region still scores > 8.5 pp in
  10% of cases at this depth, so single-region estimates are imprecise.
- **Profile:** the median evaluation-set effect is 12.8 pp, and only about 5%
  of rows are below 5 pp in the evaluation set.
  - With the 2/5/10 pp edges, regions binned <2, 2–5 and 5–10 had evaluation
    medians of 8.5, 8.4 and 9.7 pp, so those bins can't be told apart.
  - With 10/20 pp edges they have 9.3, 13.5 and 21.8 pp, with the same bin in
    both sets 66% of the time and within one bin 99%.
- **Covariation:** agreement between the sets is weaker (Spearman 0.45,
  profile 0.68). A single split at 0.2 is the most a 30% selection set can
  support: 71% land in the same bin.
- A 50% selection split helped little (Spearman 0.68 → 0.72) and lowered
  evaluation signal-to-noise from 4.6 to 3.5. We kept 30%.
- About 31% of profile difference energy is a region-wide offset.

## Results on the full test split

Contrast-centered, all contrasts (from `differential-evaluation-test-split/strata.tsv`):

| Metric | Bin | Regions | Pooled r / ceiling | Between r / ceiling | Within r / ceiling | Between slope |
|---|---|---|---|---|---|---|
| Profile | <10 pp | 22,397 | 0.13 / 0.77 | 0.22 / 0.84 | 0.11 / 0.76 | 0.26 |
| | 10–20 pp | 24,782 | 0.22 / 0.84 | 0.34 / 0.91 | 0.19 / 0.82 | 0.47 |
| | ≥20 pp | 12,195 | 0.40 / 0.91 | 0.52 / 0.96 | 0.35 / 0.89 | 0.95 |
| Covariation | <0.2 | 24,462 | 0.16 / 0.69 | 0.21 / 0.77 | 0.15 / 0.68 | – |
| | ≥0.2 | 20,610 | 0.27 / 0.76 | 0.38 / 0.87 | 0.24 / 0.74 | – |

Regions are counted once per bin they appear in for any contrast, so the bin
counts sum to more than the 25,860 analyzed regions.

**Per contrast, top profile bin:**
- **Between-region r:** 0.40 (C2C12 − NP) to 0.62 (MEL − NP).
- **Within-region r:** 0.22 (MEL − C2C12) to 0.45 (MEL − ES).
- **Best and worst contrasts:** MEL − ES and NP − ES are generally best;
  MEL − C2C12 and C2C12 − NP weakest.

**Cell-type-wide offsets (profile, observed vs predicted):** the model gets
these largely right.

| Contrast | Observed | Predicted |
|---|---|---|
| MEL − C2C12 | −5.8 pp | −3.0 pp |
| MEL − NP | +5.7 | +4.9 |
| MEL − ES | +7.0 | +8.1 |
| C2C12 − NP | +11.5 | +7.9 |
| C2C12 − ES | +12.8 | +11.1 |
| NP − ES | +1.3 | +3.2 |

Before centering, the pooled between-region r in the top bin was 0.71, above
every single contrast. After centering it is 0.52, inside the per-contrast
range.

**Examples:** `examples.pdf` (4 random regions per contrast from the top
profile bin) repeatedly shows nearly identical predicted profiles and
correlation matrices for both cell types, where the observed data differ by
20–30 pp.

## Implementation notes

- **Workers:** `differential.num_workers=8` scores regions in single-threaded
  spawned processes. Output is identical to a single-process run (checked on
  2,000 real regions), and the full test split takes about 6 minutes. A region
  costs about 66 ms single-threaded (12 ms of that is reading).
- **Seeding bug, fixed:** PyTorch's CPU generator keeps only a seed's low
  32 bits. The first version seeded with `seed << 32 + crc32(region_id)`, so
  changing `seed` did nothing. Seed and region ID are now hashed together with
  `crc32(f"{seed}:{region_id}")`.
- **Tests:** `tests/test_evaluation_differential.py` has 28 tests. They include a
  perfect, independent model that reaches every ceiling, offset invariance,
  bootstrap correctness against explicit resampling, and worker equivalence.
  The full suite passes (387 passed, 10 skipped).

## Commands

The environment is the `load_smf` alias, or
`$GROUP_HOME/uv-envs/SMFNet-cu124/bin/python`. The runs used Sherlock compute
node sh04-03n07 in job 46777507 (8 CPUs). Paths are relative to the SMFNet repo.

```bash
OUT_ROOT="$SCRATCH/evaluate_SMFNet_cell_lines"
python scripts/differential_metrics.py \
  "samples=[\"MEL=$OUT_ROOT/MEL/test_sample/samples.h5\",\"C2C12=$OUT_ROOT/C2C12/test_sample/samples.h5\",\"NP=$OUT_ROOT/NP/test_sample/samples.h5\",\"ES=$OUT_ROOT/ES/test_sample/samples.h5\"]" \
  differential.num_workers=8 hydra.run.dir="$OUT_ROOT/differential_test_cell_lines_full"

# Bin calibration (each under 2 minutes)
python analyses/differential_bins_20261006/effect_sizes.py --output "$SCRATCH/differential_bins_20261006/sel0.3" --selection-fraction 0.3
python analyses/differential_bins_20261006/effect_sizes.py --output "$SCRATCH/differential_bins_20261006/sel0.5" --selection-fraction 0.5
python analyses/differential_bins_20261006/control.py --output "$SCRATCH/differential_bins_20261006/control"
python analyses/differential_bins_20261006/summarize.py "$SCRATCH/differential_bins_20261006/sel0.3/effects.tsv"

# Tests
CUDA_VISIBLE_DEVICES="" MPLBACKEND=Agg python -m pytest -q tests/test_evaluation_differential.py
```

The full-run folder's `strata.tsv` and PDFs predate contrast centering. The
copied `strata.tsv` here (`strata_centered.tsv` in that folder) was
recomputed from its `pairwise.tsv` with
`differential_strata(pairs, DifferentialConfig())`. Rerunning the command above
regenerates everything with centering.

# Related Notes

- [Chr16 feasibility of differential effects](../2026-10-05/chr16-differential-effects-feasibility.md): Showed that 80–90% of regions differ significantly at full depth, so significance can't serve as a filter. It recommended selecting and scoring on independent molecules, which this design implements.
- [Covariation differences diagnosis](covariation-differences-diagnostics.md): Explains why covariation differences are predicted weakly.
- [Positional baselines and noise-matched correlation of correlations](../2026-09-29/positional-covariation-baselines-and-noise-matched-corr-of-corr.md): Background on noise matching and on why correlation-of-correlations reports no slope; the same reasoning is applied here.

# Open Questions

- How much of the between-region and offset agreement is biological rather
  than dataset-level (conversion efficiency, depth, handling)? The tables can't
  separate the two, and within-region scores are the most defensible.
- Would a validation-split run, or models trained jointly across cell types,
  change the within-region gap?
- Whether to pick the main-figure metric (within-region r with ceiling) for
  the paper, with between-region and offsets in the supplement.

# Sources

- SMFNet `scripts/differential_metrics.py`, `scripts/differential_metrics.md`
- SMFNet `src/smf_net/evaluation/differential.py`, `differential_plotting.py`
- SMFNet `tests/test_evaluation_differential.py`
- SMFNet `analyses/differential_bins_20261006/` (not tracked in SMFNet; copies: [effect_sizes.py](differential-evaluation-test-split/scripts/effect_sizes.py), [control.py](differential-evaluation-test-split/scripts/control.py), [summarize.py](differential-evaluation-test-split/scripts/summarize.py))
- Copied artifacts: [strata.tsv](differential-evaluation-test-split/strata.tsv), [examples.tsv](differential-evaluation-test-split/examples.tsv), [summary.tsv](differential-evaluation-test-split/summary.tsv), [analysis_metadata.json](differential-evaluation-test-split/analysis_metadata.json), [bin calibration (30%)](differential-evaluation-test-split/bin_calibration_selection30.txt), [bin calibration (50%)](differential-evaluation-test-split/bin_calibration_selection50.txt), [same-cell-type control](differential-evaluation-test-split/same_cell_control.txt), [strata figure](differential-evaluation-test-split/strata_page.png), [region mean differences, profile](differential-evaluation-test-split/region_mean_differences_profile.png)
- Full outputs (purged after 90 days unmodified): `/scratch/users/diamant/evaluate_SMFNet_cell_lines/differential_test_cell_lines_full/`
