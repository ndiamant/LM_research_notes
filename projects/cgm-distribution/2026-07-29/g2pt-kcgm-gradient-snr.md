---
title: G2PT kCGM gradient SNR across batch sizes and kernels
date: 2026-07-29
project: cgm-distribution
agent: codex
status: complete
sources:
  - /scratch/users/diamant/g2pt_snr_experiments/gradient_snr_comparison.tsv
  - /scratch/users/diamant/g2pt_snr_experiments/global_gradient_snr_wide.tsv
  - /scratch/users/diamant/g2pt_snr_experiments/descriptors_energy/run_config.json
  - /scratch/users/diamant/g2pt_snr_experiments/descriptors_energy_no_loo/run_config.json
  - /scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median/run_config.json
  - /scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median_no_loo/run_config.json
tags:
  - g2pt
  - kcgm
  - gradient-snr
  - mmd
  - energy-distance
  - rbf
  - loo-baseline
---

# Summary

We measured the signal-to-noise ratio (SNR) of the first-batch kCGM
score-function gradient at the unchanged pretrained G2PT weights. The
experiment used scaled molecular descriptor features and compared the energy
kernel with an RBF mixture whose bandwidths were selected from the target
median-distance heuristic.

Global gradient SNR increased consistently with batch size. With the LOO
baseline, it reached 0.1108 for energy and 0.1206 for median-heuristic RBF at
batch size 192. Without LOO, the corresponding values were 0.0557 and 0.0591.
Across batch sizes, LOO reduced gradient noise power by a factor of about
1.6--2.5 and substantially increased gradient SNR. The scalar MMD-squared
estimate was identical between paired LOO and no-LOO runs, as expected because
the baseline only changes the gradient estimator.

# Key Points

- Each gradient was evaluated at
  `xchen16/g2pt-guacamol-small-bfs`; no optimizer step changed the model.
- Seven batch sizes were tested: 8, 16, 32, 64, 96, 128, and 192.
- Each batch-size/kernel/baseline combination used 50 independent gradient
  batches.
- Energy, RBF, LOO, and no-LOO runs used identical target samples and replicate
  seeds, making all four conditions paired with respect to model sampling
  randomness.
- The target was the populated committed version of `G2PT/abx_smiles.csv`.
  The working-tree copy was empty, so it was not used directly.
- The requested target size was 500, but the ABX target contained 174
  molecules; all 174 were used.
- Each kernel was evaluated both with and without leave-one-out correction.
  Kernel auto-scaling, gradient clipping, KL regularization, and optimizer
  updates were omitted.
- Four backward chunks were used. An initial unchunked attempt reached about
  56 GB at batch size 96, while the chunked run used about 23 GB at batch size
  192.
- All 1,400 final gradient replicates and all summary statistics were finite.

# Details

## Estimators

For \(R=50\) independent batch-gradient estimates, let
\(\bar g_i\) and \(s_i^2\) be the sample mean and unbiased sample variance of
parameter coordinate \(i\). The reported signal and noise components were

\[
\widehat S = \sum_i \left(\bar g_i^2 - \frac{s_i^2}{R}\right),
\qquad
\widehat N = \sum_i s_i^2.
\]

Both components are unbiased for squared mean-gradient norm and total gradient
variance, respectively. The reported global SNR,
\(\widehat S/\widehat N\), is a plug-in ratio and is not itself generally
unbiased.

Coordinate estimates used

\[
\widehat{\mathrm{SNR}}_i =
\frac{\bar g_i^2 - s_i^2/R}{s_i^2}.
\]

Coordinates with zero estimated noise were excluded from parameter-SNR
quantiles. Negative coordinate estimates were retained rather than clamped,
because clamping would introduce positive bias.

## Configuration

- Model: `xchen16/g2pt-guacamol-small-bfs`
- Trainable parameters: 10,929,408
- Device and precision: NVIDIA H100, CUDA BF16
- Feature: 11 RDKit/QED/SA descriptors scaled by target standard deviation
- Target size: 174 molecules
- Base seed: 0
- Replicate seed: `seed + replicate + 1`
- Gradient batches per size: 50
- Backward chunks: 4
- MMD coefficient correction: evaluated with LOO and with `--no_loo`
- Kernel scale: 1.0
- Gradient clipping: none

The RBF kernel was a five-bandwidth mixture using target median-distance
factors `[0.25, 0.5, 1.0, 2.0, 4.0]`. The resulting sigmas were:

```text
[0.916054, 1.832108, 3.664216, 7.328432, 14.656864]
```

The commands were run from
`/home/users/diamant/repos/cgm_distribution`. The four final runs used this
command shape with the following variable settings:

| Condition | `KERNEL` | `OUT_DIR` | `LOO_FLAG` |
| --- | --- | --- | --- |
| Energy LOO | `energy` | `/scratch/users/diamant/g2pt_snr_experiments/descriptors_energy` | empty |
| Energy no-LOO | `energy` | `/scratch/users/diamant/g2pt_snr_experiments/descriptors_energy_no_loo` | `--no_loo` |
| RBF LOO | `rbf` | `/scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median` | empty |
| RBF no-LOO | `rbf` | `/scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median_no_loo` | `--no_loo` |

```bash
git show HEAD:G2PT/abx_smiles.csv |
  MPLCONFIGDIR=/home/users/diamant/repos/cgm_distribution/.mplconfig \
  HF_HUB_DISABLE_PROGRESS_BARS=1 \
  /home/groups/btrippe/diamant/miniforge/envs/g2pt/bin/python \
  G2PT/estimate_g2pt_gradient_snr.py \
    --out_dir "$OUT_DIR" \
    --feature descriptors \
    --kernel "$KERNEL" \
    --batch_sizes 8 16 32 64 96 128 192 \
    --num_batches 50 \
    --batch_chunks 4 \
    --n_hstar 500 \
    --target_csv /dev/stdin \
    --model_name xchen16/g2pt-guacamol-small-bfs \
    --device cuda \
    --bf16 \
    $LOO_FLAG
```

## Global Gradient SNR

| Batch size | Energy LOO | Energy no-LOO | RBF LOO | RBF no-LOO |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 0.002581 | -0.000203 | 0.002571 | 0.000095 |
| 16 | 0.006321 | 0.003051 | 0.006550 | 0.003362 |
| 32 | 0.015687 | 0.005718 | 0.016095 | 0.006224 |
| 64 | 0.032777 | 0.010149 | 0.034933 | 0.011309 |
| 96 | 0.053390 | 0.022924 | 0.057511 | 0.024871 |
| 128 | 0.068528 | 0.032112 | 0.073559 | 0.034549 |
| 192 | 0.110753 | 0.055741 | 0.120598 | 0.059127 |

The global gradient SNR was approximately monotone in batch size for both
kernels and both baseline settings. The negative energy no-LOO estimate at
batch size 8 is permitted by the unbiased finite-replicate signal-power
estimator and indicates that signal was not distinguishable from noise there.
The RBF advantage remained modest relative to the much larger effect of LOO.

## Scalar MMD-Squared SNR

| Batch size | Energy | RBF |
| ---: | ---: | ---: |
| 8 | 9.314 | 8.634 |
| 16 | 12.730 | 12.237 |
| 32 | 31.858 | 30.794 |
| 64 | 70.482 | 71.680 |
| 96 | 69.663 | 72.865 |
| 128 | 98.332 | 101.520 |
| 192 | 142.948 | 149.497 |

The paired LOO and no-LOO runs had exactly equal per-replicate MMD-squared
estimates and therefore identical scalar MMD-squared SNR. This was an
experiment validation check, not merely agreement within a tolerance.

## LOO Noise Reduction

The table reports no-LOO gradient noise power divided by LOO gradient noise
power. Values above one measure the variance reduction from the baseline.

| Batch size | Energy | RBF |
| ---: | ---: | ---: |
| 8 | 1.927 | 1.567 |
| 16 | 2.218 | 2.144 |
| 32 | 2.240 | 2.140 |
| 64 | 2.471 | 2.422 |
| 96 | 2.210 | 2.213 |
| 128 | 2.312 | 2.314 |
| 192 | 2.240 | 2.268 |

Except for RBF at batch size 8, LOO reduced gradient noise power by more than a
factor of two.

## Parameter-SNR Quantiles

Median coordinate estimates:

| Batch size | Energy LOO | Energy no-LOO | RBF LOO | RBF no-LOO |
| ---: | ---: | ---: | ---: | ---: |
| 8 | -0.01012 | -0.01060 | -0.01016 | -0.01051 |
| 16 | -0.00939 | -0.01000 | -0.00944 | -0.00991 |
| 32 | -0.00894 | -0.00967 | -0.00888 | -0.00957 |
| 64 | -0.00786 | -0.00931 | -0.00767 | -0.00921 |
| 96 | -0.00638 | -0.00856 | -0.00619 | -0.00842 |
| 128 | -0.00578 | -0.00801 | -0.00553 | -0.00782 |
| 192 | -0.00354 | -0.00645 | -0.00306 | -0.00622 |

99th-percentile coordinate estimates:

| Batch size | Energy LOO | Energy no-LOO | RBF LOO | RBF no-LOO |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 0.1208 | 0.1179 | 0.1209 | 0.1177 |
| 16 | 0.1489 | 0.1264 | 0.1515 | 0.1276 |
| 32 | 0.2201 | 0.1515 | 0.2212 | 0.1544 |
| 64 | 0.3560 | 0.1860 | 0.3742 | 0.1935 |
| 96 | 0.5590 | 0.2586 | 0.5627 | 0.2605 |
| 128 | 0.7088 | 0.3264 | 0.7290 | 0.3353 |
| 192 | 1.0475 | 0.4908 | 1.1088 | 0.5091 |

The median debiased coordinate estimate remained slightly negative at every
batch size. This does not imply negative true SNR; it indicates that at least
half of parameter coordinates had signal too weak to separate from noise with
50 replicates. With LOO, the upper tail improved substantially and the 99th
percentile exceeded 1 at batch size 192. Without LOO, it remained near 0.5.

## Execution History

1. A two-replicate energy smoke test at batch size 4 verified model loading,
   descriptor extraction, backward propagation, finite outputs, and unchanged
   model parameters.
2. A two-replicate RBF smoke test verified median-heuristic bandwidth
   selection and serialization.
3. The first full energy attempt used one backward chunk. It was stopped during
   batch size 96 after memory reached about 56 GB.
4. Both kernels were restarted from the beginning with four backward chunks.
   All final batch sizes and replicates completed successfully.
5. Validation checked seven summary rows and 350 replicate rows per kernel,
   finite numeric values, identical target snapshots, matching seed schedules,
   and the expected RBF bandwidth metadata.
6. Matched energy and RBF runs were repeated with `--no_loo`, producing 700
   additional gradients.
7. Four-condition validation checked 28 summary rows and 1,400 replicate rows,
   finite outputs, identical target snapshots, matching RBF bandwidths, and
   exact equality of paired LOO/no-LOO scalar MMD-squared estimates.

## Artifacts

Combined summary:

```text
/scratch/users/diamant/g2pt_snr_experiments/gradient_snr_comparison.tsv
```

Wide global-gradient SNR table for plotting:

```text
/scratch/users/diamant/g2pt_snr_experiments/global_gradient_snr_wide.tsv
```

Per-condition directories:

```text
/scratch/users/diamant/g2pt_snr_experiments/descriptors_energy
/scratch/users/diamant/g2pt_snr_experiments/descriptors_energy_no_loo
/scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median
/scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median_no_loo
```

Each directory contains `gradient_snr.tsv`, `batch_metrics.tsv`,
`run_config.json`, `hstar_smiles.csv`, and `run.log`.

# Related Notes

- No earlier CGM-distribution notes were present when this note was created.

# Open Questions

- Add confidence intervals or repeat the experiment with more than 50 batches
  before treating the small energy-versus-RBF differences as conclusive.
- Restore or regenerate the working-tree `G2PT/abx_smiles.csv`; it currently
  contains only the header and would make a direct rerun fail.
- Decide whether future production runs should standardize on four backward
  chunks for memory predictability.
- Extend the same paired design to other kernels and target feature sets.

# Sources

- Combined results:
  `/scratch/users/diamant/g2pt_snr_experiments/gradient_snr_comparison.tsv`
- Wide plotting table:
  `/scratch/users/diamant/g2pt_snr_experiments/global_gradient_snr_wide.tsv`
- Energy LOO configuration and logs:
  `/scratch/users/diamant/g2pt_snr_experiments/descriptors_energy/`
- Energy no-LOO configuration and logs:
  `/scratch/users/diamant/g2pt_snr_experiments/descriptors_energy_no_loo/`
- RBF LOO configuration and logs:
  `/scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median/`
- RBF no-LOO configuration and logs:
  `/scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median_no_loo/`
- Local implementation:
  `/home/users/diamant/repos/cgm_distribution/G2PT/estimate_g2pt_gradient_snr.py`
  and `/home/users/diamant/repos/cgm_distribution/cgm/gradient_snr.py`
- [CGM-distribution repository](https://github.com/ndiamant/cgm_distribution)
- [G2PT calibration submission configuration](https://github.com/ndiamant/cgm_distribution/blob/main/G2PT/submit_calibrate_g2pt_abx_lambda_sweep_v2.sh)
