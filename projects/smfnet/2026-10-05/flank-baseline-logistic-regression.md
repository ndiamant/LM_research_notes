---
title: Flank baseline for profile metrics, logistic regression instead of a k-mer table
date: 2026-10-05
project: smfnet
agent: Claude Code
status: draft
sources:
  - SMFNet repo, src/smf_net/evaluation/flank_baselines.py (branch nd/more_assays, uncommitted at time of writing)
  - SMFNet repo, src/smf_net/evaluation/{metrics,reporting,metrics_plotting}.py
  - /scratch/users/diamant/gch_rc_analysis/ (strand_symmetry_control.py, plot_context_bias.py, kmer_baseline_check.py, context_bias_MEL_flank2.{png,tsv})
  - User discussion, 2026-10-04 to 2026-10-05
tags: [smfnet, evaluation, baselines, sequence-context, strand-bias, gch, logistic-regression]
---

# Summary

SMFNet's metrics now include a **flank baseline**. It predicts each site's
mean methylation call from its flanking bases alone, using a logistic
regression on one-hot flanks fit on the training chromosomes. It gives the
reference a model must beat to show it learned more than local sequence
context. We chose it over a k-mer lookup table because an additive model has
no rare-context problem, needs no backoff or shrinkage, and allows wider
context cheaply. The available evidence suggests the effects are close to
additive. It sits alongside the sequence-blind positional baselines from
[the 2026-09-29 note](../2026-09-29/positional-covariation-baselines-and-noise-matched-corr-of-corr.md).

# Key Points

- **Motivation.** In the GCH bait-capture SMF data, methylation rate depends
  on flanking context in a strand-specific way.
  - A context and its reverse complement differ by up to about 11 percentage
    points (MEL, ±2 bp flanks).
  - Across all contexts, the correlation between a context's rate and its
    reverse complement's is −0.50. The same context in two halves of the
    regions correlates at 0.997.
  - The likely cause is strand-specific capture before bisulfite conversion.
    This is inferred, not confirmed.
  - A model can therefore earn profile correlation from local context alone.
- **Model.** logit(site mean) = bias + Σ over flank positions of
  w[position, base].
  - Defaults: ±5 bp flanks and a small L2 penalty, fit with LBFGS.
  - Labels are soft: each training region's site mean counts once.
  - `cg_gc` fits separate models for CpG, GpC and GCG cytosines.
- **Strand.** Flanks are read on the plus strand and never tied to their
  reverse complement, so the strand bias is part of the baseline instead of
  credit for the model.
- **Missing flanks.** An unknown or out-of-window flank contributes its
  position's frequency-weighted mean effect, which is what a centred one-hot
  encoding gives a missing base. Window edges need no special handling.
- **Metrics.** Profile correlations only:
  - `pearson/spearman_flank_baseline`: the prediction against the observed
    site means.
  - `*_flank_residual`: model and observed means after subtracting the
    prediction.
  - `*_flank_residual_matched{,_pseudo_replicate,_shuffle}`: the same on
    matched reads, which gives a replicate ceiling.
- **No covariation version.** Sampling reads from per-site probabilities makes
  sites independent, so the expected pair correlations are all zero.
- **Validation check.** Fit on 3,000 MEL training regions, scored on 1,000
  validation regions, at sites with ≥10 reads. Median per-region Pearson r was
  0.17 at ±2 bp, 0.22 at ±5 bp and 0.22 at ±8 bp, so gains level off around
  ±5 bp.

# Details

**Why not a k-mer table.** A table of P(methylated | k-mer) is the obvious
first design. However:

- **It needs extra machinery.** With ±3 bp flanks there are 4,096 contexts per
  site type, and many are rare. The table then needs backoff or shrinkage
  toward shorter contexts, with a prior strength to tune. Window-edge sites
  also need fallback to shorter contexts.
- **Its size limits the context.** Every extra flank base multiplies the table
  by 4, so in practice the table stays narrow.

The additive model avoids both problems:

- **Few parameters.** It has 4 × 2F weights plus a bias: 40 parameters at
  F = 5.
- **Rare contexts are no problem.** Each weight is learned from every site
  carrying that base at that position.
- **Missing bases are built in.** The centred encoding gives unknown bases a
  principled default.

**Cost.** The model is blind to interactions between flank positions. The ±1
bp GCH rates (fraction of 1s, MEL) suggest this costs little:

| Context | Rate | Change |
|---|---|---|
| AGCA | 0.605 | baseline |
| TGCA | 0.636 | +0.031 for T on the left |
| AGCT | 0.622 | +0.017 for T on the right |
| TGCT (observed) | 0.660 | |
| TGCT (additive prediction) | 0.653 | 0.605 + 0.031 + 0.017 |

The prediction misses by under a point. Only this one ±1 bp quartet was
checked. If residuals later show interactions matter, one-hots for adjacent
base pairs can be added without changing anything else.

**Implementation.**

- `train_baselines_for` loads the training split once, through
  `load_training_split`, and fits both the positional and flank baselines.
- The flank fit is subsampled to at most 2M sites per site type and takes
  seconds.
- The fit is saved as `flank_baselines.npz` next to `positional_baselines.npz`.
- Settings are `MetricsConfig.include_flank_baselines`, `flank_baseline_bp`,
  `flank_baseline_l2` and `flank_baseline_max_sites`.
- It was first named "k-mer baseline" and renamed to "flank", because it isn't
  a k-mer table.

Reproduce the validation check, from the SMFNet repo root on a compute node:

```bash
CUDA_VISIBLE_DEVICES="" $GROUP_HOME/uv-envs/SMFNet-cu124/bin/python \
  /scratch/users/diamant/gch_rc_analysis/kmer_baseline_check.py
```

# Related Notes

- [Positional baselines and noise-matched correlation-of-correlations](../2026-09-29/positional-covariation-baselines-and-noise-matched-corr-of-corr.md):
  The sequence-blind baselines this one complements. It shares their training
  data loading and their soft-label, one-count-per-region-site convention.

# Open Questions

- How much of a trained model's headline Pearson does the flank baseline
  explain? The residual columns have not yet been computed on real samples.
- Do adjacent-base interactions matter at wider flanks? Only the ±1 bp
  additivity check has been done.
- Is ±5 bp the best width? It was chosen from one coarse MEL check, not tuned
  on validation data.
- Does the strand bias come from capture? The proposed check is to count
  plus- and minus-strand reads at the `CollapseStrandsSM` step of the GCH h5
  pipeline.

# Sources

- SMFNet `src/smf_net/evaluation/flank_baselines.py` and its tests in
  `tests/test_evaluation_flank_baselines.py`.
- Strand-asymmetry analyses: `/scratch/users/diamant/gch_rc_analysis/strand_symmetry_control.py`,
  `plot_context_bias.py` and `context_bias_MEL_flank2.png`.
- Validation check: `/scratch/users/diamant/gch_rc_analysis/kmer_baseline_check.py`.
