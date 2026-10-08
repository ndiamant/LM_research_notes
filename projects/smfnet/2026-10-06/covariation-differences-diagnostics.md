---
title: Why SMFNet predicts cell-type differences in pair correlations weakly
date: 2026-10-06
project: smfnet
agent: Claude Code
status: draft
sources:
  - SMFNet repo, analyses/covariation_diagnostics_20261006/diagnose.py (branch nd/more_assays at 1151cbb, uncommitted)
  - /scratch/users/diamant/covariation_diagnostics_20261006/ (diagnostics.tsv, metadata.json)
  - covariation-differences-diagnostics/ in this folder (copied tables)
  - User discussion, 2026-10-06
tags: [smfnet, covariation, pair-correlations, cell-types, differential-effects, diagnostics, noise-correction]
---

# Summary

Within-molecule pair correlations are predicted very well in absolute terms:
noise-corrected r is 0.92–0.93 per cell type, with predicted spread 87–93% of
observed. But that structure is mostly shared across cell types, and the four
models make it even more alike than the data. Observed cell types correlate
0.85–0.89 with one another, predicted cell types 0.96–0.98. So the models keep
roughly a quarter of the cell-type-specific variation in pair correlations,
and their predicted covariation differences are 38–48% of the observed size
and only weakly aligned (true r 0.21–0.41). This explains the weak covariation
scores in
[the differential evaluation note](differential-evaluation-held-out-molecules.md)
and the nearly blank predicted-difference matrices in its examples.

# Key Points

- **Ruled out:** the models do not just shrink dependence overall. Absolute
  covariation is predicted better than absolute profiles.
- **Shared structure dominates.** Observed covariation is about 87% shared
  across cell types, compared with 57–74% for profiles.
- **Predictions are too alike.** In the predictions the cell-type-specific
  room shrinks from about 0.13 (1 − 0.87) to about 0.03 (1 − 0.97). Profiles
  lose much less (similarity 0.57–0.74 observed vs 0.69–0.77 predicted).
- **Too small:** predicted covariation differences are 0.38–0.48 of their
  observed size. Profile differences are 0.59–0.78.
- **Real signal is missed.** The observed differences are resolvable (noise
  ceilings about 0.7–0.8 in the strata), so the gap is in the models, not the
  data.

# Details

## Method

`analyses/covariation_diagnostics_20261006/diagnose.py` drew 1,500 random
shared test regions (seed 7); 1,205 had at least 10 eligible sites and pairs.
It used `load_region_features` with `DifferentialConfig()` defaults, so it
sees the same molecule split, eligible sites and pairs as
`scripts/differential_metrics.py`. Observed values come from evaluation
molecules.

Every variance is corrected by subtracting sampling noise: binomial for site
means, a 50-draw whole-molecule bootstrap for pair correlations. Correlations
are disattenuated. Three comparisons, for site means (profile) and pair
correlations (covariation):

1. **Absolute:** predicted vs observed values per cell type. Amplitude is the
   predicted / observed noise-corrected standard deviation.
2. **Differences:** predicted vs observed `cell_a − cell_b`, with the same
   amplitude ratio; 1 means the right size.
3. **Similarity:** the noise-corrected correlation between two cell types'
   values, computed separately on observed and on predicted.

## Results

| | Profile | Covariation |
|---|---|---|
| Absolute, true r | 0.64–0.77 | 0.92–0.93 |
| Absolute, amplitude | 0.66–0.86 | 0.87–0.93 |
| Differences, true r | 0.19–0.41 | 0.21–0.41 |
| Differences, amplitude | 0.59–0.78 | 0.38–0.48 |
| Similarity between cell types, observed | 0.57–0.74 | 0.85–0.89 |
| Similarity between cell types, predicted | 0.69–0.77 | 0.96–0.98 |

Per-contrast values are in
[diagnostics.tsv](covariation-differences-diagnostics/diagnostics.tsv). For
covariation, MEL − ES is the best contrast (true r 0.41, amplitude 0.47) and
C2C12 − NP the worst (0.21, 0.44).

## Possible causes (speculative, not tested)

- **Shared structure dominates training.** Most within-molecule dependence is
  probably footprint geometry shared by all cell types, such as
  nucleosome-sized protection. It dominates what a generative model fits, while
  the cell-specific part is a small residual.
- **Sequence-only inputs.** The models are trained separately per cell type,
  each with DNA sequence as its only input. Predicting a cell-specific
  co-occupancy pattern means inferring that cell's factor configuration from
  sequence.
- **Weak training signal.** Pair statistics need many reads per region to
  estimate, so cell-specific dependence gives a weak signal.
- **Uneven training data.** Data sizes differ between cell types, e.g. the
  C2C12 validation metrics scored 3,832 regions vs MEL's 17,640.

## Commands

```bash
# From the SMFNet repo with load_smf active, on a compute node (~1 minute on 8 CPUs)
python analyses/covariation_diagnostics_20261006/diagnose.py \
  --output "$SCRATCH/covariation_diagnostics_20261006" --regions 1500
```

# Related Notes

- [Differential evaluation with held-out molecules](differential-evaluation-held-out-molecules.md): Defines the split, eligibility and noise estimates reused here, and reports the weak covariation strata this note explains.
- [Chr16 feasibility of differential effects](../2026-10-05/chr16-differential-effects-feasibility.md): Showed that observed correlation patterns differ in most regions, so the differences the models miss are real and testable.

# Open Questions

Suggested checks, not yet run:

- **Profile vs covariation:** are covariation differences predicted better
  where profile differences are? That would suggest the correlation changes
  mostly follow accessibility changes.
- **Training data size:** does training-data size per cell type predict a
  model's difference amplitude?
- **Distance:** split by pair separation, under vs over 150 bp, to see whether
  the missing part is local footprint structure or longer-range co-occupancy.

# Sources

- SMFNet `analyses/covariation_diagnostics_20261006/diagnose.py` (not tracked in SMFNet; copies: [diagnose.py](covariation-differences-diagnostics/scripts/diagnose.py), [one_vs_rest.py](covariation-differences-diagnostics/scripts/one_vs_rest.py))
- Copied artifacts: [diagnostics.tsv](covariation-differences-diagnostics/diagnostics.tsv), [metadata.json](covariation-differences-diagnostics/metadata.json)
- Full outputs (purged after 90 days unmodified): `/scratch/users/diamant/covariation_diagnostics_20261006/`
