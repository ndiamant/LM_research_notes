---
title: G2PT kCGM gradient SNR across batch sizes and kernels
date: 2026-07-29
project: cgm-distribution
agent: codex
status: complete
sources:
  - /scratch/users/diamant/g2pt_snr_experiments/gradient_snr_comparison.tsv
  - /scratch/users/diamant/g2pt_snr_experiments/descriptors_energy/run_config.json
  - /scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median/run_config.json
tags:
  - g2pt
  - kcgm
  - gradient-snr
  - mmd
  - energy-distance
  - rbf
---

# Summary

We measured the signal-to-noise ratio (SNR) of the first-batch kCGM
score-function gradient at the unchanged pretrained G2PT weights. The
experiment used scaled molecular descriptor features and compared the energy
kernel with an RBF mixture whose bandwidths were selected from the target
median-distance heuristic.

Global gradient SNR increased consistently with batch size. At batch size 192,
it was 0.1108 for energy and 0.1206 for median-heuristic RBF. The RBF kernel
gave modestly higher global gradient SNR for every batch size except 8. The
scalar MMD-squared estimate had much higher SNR than its parameter gradient,
reaching 142.9 for energy and 149.5 for RBF at batch size 192.

# Key Points

- Each gradient was evaluated at
  `xchen16/g2pt-guacamol-small-bfs`; no optimizer step changed the model.
- Seven batch sizes were tested: 8, 16, 32, 64, 96, 128, and 192.
- Each batch-size/kernel combination used 50 independent gradient batches.
- Energy and RBF runs used identical target samples and replicate seeds, making
  the kernel comparison paired with respect to model sampling randomness.
- The target was the populated committed version of `G2PT/abx_smiles.csv`.
  The working-tree copy was empty, so it was not used directly.
- The requested target size was 500, but the ABX target contained 174
  molecules; all 174 were used.
- Leave-one-out correction was enabled. Kernel auto-scaling, gradient clipping,
  KL regularization, and optimizer updates were omitted.
- Four backward chunks were used. An initial unchunked attempt reached about
  56 GB at batch size 96, while the chunked run used about 23 GB at batch size
  192.
- All 700 final gradient replicates and all summary statistics were finite.

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
- MMD coefficient correction: leave-one-out
- Kernel scale: 1.0
- Gradient clipping: none

The RBF kernel was a five-bandwidth mixture using target median-distance
factors `[0.25, 0.5, 1.0, 2.0, 4.0]`. The resulting sigmas were:

```text
[0.916054, 1.832108, 3.664216, 7.328432, 14.656864]
```

The final runs used the following command shape, once with `KERNEL=energy` and
once with `KERNEL=rbf`:

```bash
git show HEAD:G2PT/abx_smiles.csv |
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
    --bf16
```

## Global Results

| Batch size | Energy gradient SNR | RBF gradient SNR | Energy MMD-squared SNR | RBF MMD-squared SNR |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 0.002581 | 0.002571 | 9.314 | 8.634 |
| 16 | 0.006321 | 0.006550 | 12.730 | 12.237 |
| 32 | 0.015687 | 0.016095 | 31.858 | 30.794 |
| 64 | 0.032777 | 0.034933 | 70.482 | 71.680 |
| 96 | 0.053390 | 0.057511 | 69.663 | 72.865 |
| 128 | 0.068528 | 0.073559 | 98.332 | 101.520 |
| 192 | 0.110753 | 0.120598 | 142.948 | 149.497 |

The global gradient SNR was approximately monotone in batch size for both
kernels. The RBF advantage was small at low batch sizes and reached about 9%
at batch size 192. The scalar MMD-squared SNR was noisier as a function of
batch size, including a small decrease from 64 to 96.

## Parameter-SNR Quantiles

| Batch size | Energy median | RBF median | Energy q99 | RBF q99 |
| ---: | ---: | ---: | ---: | ---: |
| 8 | -0.01012 | -0.01016 | 0.1208 | 0.1209 |
| 16 | -0.00939 | -0.00944 | 0.1489 | 0.1515 |
| 32 | -0.00894 | -0.00888 | 0.2201 | 0.2212 |
| 64 | -0.00786 | -0.00767 | 0.3560 | 0.3742 |
| 96 | -0.00638 | -0.00619 | 0.5590 | 0.5627 |
| 128 | -0.00578 | -0.00553 | 0.7088 | 0.7290 |
| 192 | -0.00354 | -0.00306 | 1.0475 | 1.1088 |

The median debiased coordinate estimate remained slightly negative at every
batch size. This does not imply negative true SNR; it indicates that at least
half of parameter coordinates had signal too weak to separate from noise with
50 replicates. The upper tail improved substantially, with the 99th percentile
exceeding 1 only at batch size 192.

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

## Artifacts

Combined summary:

```text
/scratch/users/diamant/g2pt_snr_experiments/gradient_snr_comparison.tsv
```

Per-kernel directories:

```text
/scratch/users/diamant/g2pt_snr_experiments/descriptors_energy
/scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median
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
- Energy configuration and logs:
  `/scratch/users/diamant/g2pt_snr_experiments/descriptors_energy/`
- RBF configuration and logs:
  `/scratch/users/diamant/g2pt_snr_experiments/descriptors_rbf_median/`
- Local implementation:
  `/home/users/diamant/repos/cgm_distribution/G2PT/estimate_g2pt_gradient_snr.py`
  and `/home/users/diamant/repos/cgm_distribution/cgm/gradient_snr.py`
- [CGM-distribution repository](https://github.com/ndiamant/cgm_distribution)
- [G2PT calibration submission configuration](https://github.com/ndiamant/cgm_distribution/blob/main/G2PT/submit_calibrate_g2pt_abx_lambda_sweep_v2.sh)
