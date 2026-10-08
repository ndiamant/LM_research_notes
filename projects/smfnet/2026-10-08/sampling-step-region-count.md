---
title: Region count for the sampling-step appendix figure
date: 2026-10-08
project: smfnet
agent: Codex
status: complete
sources:
  - /scratch/users/diamant/evaluate_SMFNet_test_set/old_full_model/metrics_v3/metrics.tsv
  - /scratch/users/diamant/SMFNet_sample_step_eval_1024_samples/steps-{32,64,96}/metrics/metrics.tsv
tags: [sampling, paired-comparisons, standard-error, appendix]
---

# Summary

Choose **2,048 identical test regions at every sampling-step setting**, with **1,024 generated reads per region**, for an appendix figure showing the sampling-step trend. Region count and reads per region are different quantities. This is a compromise for showing the trend, not a design guaranteeing significance for every adjacent step comparison.

# Key Points

- A single-condition test table has 13,744 finite regions for all three metrics. Reducing to 1,024 regions increases ordinary mean SE by 3.66× while reducing sampling work about 13.4×.
- Paired differences across the same regions are more precise than comparing separate mean error bars. For 2,048 regions, the validation-based paired Wasserstein SE is about **0.057 bp**.
- 1,024 regions reliably show the 32→64 Wasserstein gain; 2,048 better capture the corr-of-corr gain. Resolving the smaller 64→96 Wasserstein gain favors 4,096 regions, but that was not the main appendix goal.

# Details

The test table's ordinary mean SEs (full / estimated at 1,024 regions) were: profile Pearson **0.00127 / 0.00467**, corr-of-corr Pearson **0.00069 / 0.00252**, ones-run-length Wasserstein **0.0334 / 0.1222 bp**.

For paired analysis, align the 32-, 64-, and 96-step validation TSVs by unique `region` IDs (7,559 matching regions). Compute higher-step minus lower-step values. Full-data differences were:

| Steps | Profile Pearson | Corr-of-corr Pearson | Ones Wasserstein (bp) |
|---|---:|---:|---:|
| 32→64 | −0.000153 | +0.002936 | −0.373879 |
| 64→96 | +0.000446 | +0.000695 | −0.121034 |

Draw 2,000 subsets without replacement at each size, using `numpy.random.default_rng(2026)` and sizes in order 1,024, 2,048, 4,096. Use the same subset across conditions and metrics. Compute SE as sample SD (`ddof=1`) divided by square root of finite region count. Mean paired SEs were nearly identical for both contrasts:

| Regions | Profile Pearson SE | Corr-of-corr SE | Wasserstein SE (bp) |
|---|---:|---:|---:|
| 1,024 | 0.00107–0.00108 | 0.00133 | 0.081 |
| 2,048 | 0.00075–0.00076 | 0.00094 | 0.057 |
| 4,096 | 0.00053–0.00054 | 0.00066 | 0.040 |

Fraction of subsets whose normal-approximation 95% paired CI excluded zero in the improvement direction:

| Comparison | 1,024 | 2,048 | 4,096 |
|---|---:|---:|---:|
| 32→64 Wasserstein | 99.8% | 100% | 100% |
| 32→64 corr-of-corr | 59% | 91% | 99.9% |
| 64→96 Wasserstein | 32% | 55% | 94% |

These are **stability frequencies within the existing validation dataset**, not prospective test-set power. Subsets share a finite parent dataset; SEs above intentionally use the ordinary region-level convention, without finite-population correction. Genomic dependence, new sampling seeds, and subsequent metric implementation changes are not modeled.

## Submitted run

Sampling array **47026014**, tasks 0–4 = 16/32/64/96/128 steps, submitted on `akundaje`, one L40S, four CPUs, 24 GB RAM, eight-hour limit per task. Deterministic region selection in `src/smf_net/data.py` uses seed 1231241; sampling seed is 0. Checkpoint: `$SCRATCH/smf_models/epoch=381-step=173810-converted.ckpt`. Output root: `$SCRATCH/SMFNet_sample_step_eval_1024_samples_test_2048_regions`. Submission was confirmed; completion was not checked in this note. Metrics were not submitted.

Exact submission command (do not rerun while the existing array is active):

```bash
cd /home/users/diamant/repos/SMFNet
sbatch --parsable --output="$SCRATCH/SMFNet_sample_step_eval_1024_samples_test_2048_regions/logs/%A_%a.out" scripts/sample_quality_array.sbatch
```

The batch script activates the environment with `load_smf`, sets `QUALITY_MAX_REGIONS=2048` and `QUALITY_OUTDIR`, and calls `scripts/sample_quality.sh sample`. Shell syntax and `sbatch --test-only` were checked before submission. Temporary compiler caches use `$L_SCRATCH`. The sampling script's metrics stage disables positional, flank, shuffle, and model-ablation baselines.

# Related Notes

- [Run-length missingness and gap bridging](run-length-missingness-gap-bridging.md): Related interpretation of run-length metrics; metric-definition changes matter when extrapolating older validation results.

# Open Questions

- Confirm actual region identity across completed test conditions, then report paired uncertainty alongside the appendix trend.
- The existing plotting script shows per-condition means/SEs; paired-difference plotting remains to be added if desired.

# Sources

- User-provided test table: `/scratch/users/diamant/evaluate_SMFNet_test_set/old_full_model/metrics_v3/metrics.tsv`.
- Validation sweep: `/scratch/users/diamant/SMFNet_sample_step_eval_1024_samples/steps-{32,64,96}/metrics/metrics.tsv`.
- Local code: `/home/users/diamant/repos/SMFNet/scripts/{sample_quality.sh,sample_quality_array.sbatch,plot_sampling_results.py}`. Analysis used `load_smf`, pandas and NumPy; procedure and seed are recorded above. Scratch outputs are temporary, not archival copies.
