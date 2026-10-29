---
title: Frozen AlphaGenome encoder as SMFNet's DNA encoder (ONT mESC), design and val results vs baseline and ATAC
date: 2026-10-29
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - SMFNet_public commit 803eb27 on branch nd_alpha_genome (not pushed)
  - SMFNet_public/ALPHAGENOME_IMPLEMENTATION_PLAN.md (untracked)
  - https://genomicsxai.github.io/blogs/2026-002/
  - https://github.com/genomicsxai/alphagenome-pytorch (commit 72268c0)
  - https://huggingface.co/gtca/alphagenome_pytorch
tags: [alphagenome, dna-encoder, pretraining, transfer, atac, ont-mesc, sampling-metrics]
---

# Summary

Replacing SMFNet's learned DNA encoder with AlphaGenome's pretrained
convolutional encoder, **frozen**, plus a light trainable decoder to 1 bp, gives
almost the same validation loss as the baseline (0.29373 vs 0.29436). But the
sampled profiles are clearly better on the ONT mESC validation active regions:
`pearson_matched` 0.760 vs 0.731, paired per-region difference +0.028 ± 0.002.
That closes roughly 70–80% of the gap between the baseline and the
ATAC-conditioned model (0.772), without any extra input. Training is also ~17% faster per step.
Resampling the baseline with the current sampling code reproduced its old
samples within noise, so the gain is not an artifact of the sampling code
(update below). Remaining caveat: one seed per model.

# Key Points

- **Design.** AlphaGenome's encoder (stem + 6 down blocks, 90.0M params) is
  frozen and outputs features at 1, 2, 4, …, 64 bp. Each resolution goes
  through `CNNLayerNorm1D` → 1x1 conv to 512 channels → `repeat_interleave`
  back to 1 bp. The seven are summed, then one `CNNBlock` (k=9), `TakeMiddle`
  and FiLM as before. 8.8M trainable params, against 21.2M for the baseline.
  There is no dilated DNA body.
- **Validation loss** barely separates the models: baseline 0.29436,
  AlphaGenome 0.29373, ATAC 0.29328 (best epoch 377 for all three).
- **Sampled-SMF metrics** separate them clearly. The ranking on profile
  metrics and on zeros run lengths is baseline < AlphaGenome < ATAC.
  Covariation corr-of-corr improves only slightly.
- **Ones run lengths** are slightly *worse* than the baseline for AlphaGenome
  (+0.13 bp) and more so for ATAC (+0.31 bp).
- **Inferred explanation (not tested directly).** Validation loss averages
  over mask rates and is dominated by prediction from the read's own visible
  calls at low mask rates. Sampling starts from fully masked reads, so it
  depends entirely on the DNA (and ATAC) conditioning, which is where the
  pretrained encoder helps. Freezing showed no sign of underfitting.
- **Speed.** At full batch size (256 regions × 32 reads, bf16, compiled) a
  step takes 576 ms against 695 ms for the baseline, with peak memory 58.5 vs
  63.3 GiB. 400 epochs took ~33 h against 1 d 15 h for the baseline.

# Details

## Why not the blog's recipe directly

The genomicsxai post adapts AlphaGenome to 100–300 bp MPRA sequences. It keeps
the whole conv encoder, drops the transformer and decoder, and puts an MLP on
globally pooled 128 bp embeddings. It reports both frozen probing (which
generalized better out of distribution) and progressive unfreezing (1 epoch of
head-only training first). SMFNet needs a per-base embedding for FiLM, so this
note's design rebuilds 1 bp resolution from the encoder's multi-resolution
outputs instead of pooling.

An earlier plan used only the 1 bp stem in place of SMFNet's first conv layer.
It was dropped because only 3M of the encoder's ~22M parameters would have
been pretrained, with a 19 bp receptive field.

## AlphaGenome encoder facts used

- Stem: `Conv1d(4, 768, 15)` plus a residual `ConvBlock(768, 768, 5)`.
  `ConvBlock` = `RMSBatchNorm` → GELU (`sigmoid(1.702x)·x`) → weight-standardized
  conv. Each down block adds 128 channels, with max-pooling by 2 in between.
- `RMSBatchNorm` uses stored statistics and never updates them, so it is a
  fixed affine map. Batch size and batch composition don't matter.
- The encoder needs a length divisible by 64 (128 including the final pool).
  The multiple-of-2048 constraint comes only from the transformer's pair path.
  SMFNet's 2048 bp DNA input fits.
- Weights: HF `gtca/alphagenome_pytorch`, `model_all_folds.safetensors`
  (921 MB; 80 `encoder.*` tensors, kernels in bf16). They are subject to
  DeepMind's AlphaGenome Model Terms, a constraint on publishing derived
  checkpoints. That was deferred because this is an internal experiment.
- AlphaGenome's folds hold out genomic intervals, not chromosomes, so it has
  seen SMFNet's val and test chromosomes. The comparison measures transfer,
  not generalization of the pretrained part to unseen sequence.
- By resolution, the frozen feature RMS on random DNA was 4.8, 4.7, 4.3, 4.7,
  5.4, 7.1 and 11.4.

## Implementation (SMFNet_public `nd_alpha_genome`, commit 803eb27)

- `src/smf_net/alphagenome.py`: encoder modules copied from alphagenome-pytorch
  (Apache 2.0) under the same weight names. Copying keeps SMFNet checkpoints
  independent of upstream renames. `load_alphagenome_encoder_weights` loads
  strictly from the `encoder.` keys.
- `src/smf_net/cnn.py`: `AlphaGenomeDNAEncoder`, with
  `CNNConfig.dna_encoder_type` (`"cnn"` default or `"alphagenome"`) and
  `alphagenome_frozen` (default true). The frozen part runs under `no_grad`
  unless the input DNA requires grad.
- `scripts/train.py`: the `alphagenome_weights=` option, required when frozen.
  Checkpoints contain the AlphaGenome weights (465 MB each), so `load_model`
  doesn't need the weights file.
- `tests/test_alphagenome.py`: covers shapes and gradients, weight loading,
  and an exact match against upstream `SequenceEncoder` at every resolution.
- Pitfall for later fine-tuning: `configure_optimizers` skips parameters with
  `requires_grad=False`. Unfreezing later would need the pretrained parameters
  kept in the optimizer at lr 0 instead.

## Runs and commands

Environment: Sherlock job 49826490 (btrippe partition, 1× H100), with
`module load devel pixi/0.53.0` and `PIXI_CACHE_DIR=$SCRATCH/pixi-cache`, run
from `~/repos/SMFNet_public`.

Training: 400 epochs, wandb run `1x9j2po8`, output in
`/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_frozen`.

```bash
pixi run python scripts/train.py \
    hydra.run.dir=/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_frozen \
    wandb_logger.save_dir=/scratch/users/diamant/wandb wandb_logger.name=ONT_mESC_alphagenome_frozen \
    assay_type=cg_gc_accessibility \
    "data.cfg.h5_path=[/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions.h5,/scratch/groups/btrippe/ndiamant/ONT_mESC_bgRegions.h5]" \
    data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/mm10.fa \
    model.cfg.cnn_config.dna_encoder_type=alphagenome \
    alphagenome_weights=/scratch/users/diamant/alphagenome_pytorch/model_all_folds.safetensors
```

Sampling: the defaults are 1,024 reads, 32 steps, seed 0 and bf16 (1 h 26 m).
Then metrics with plots (~20 min).

```bash
R=/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_frozen
pixi run python scripts/sample.py hydra.run.dir=$R/samples_val_ac_epoch377 \
    "checkpoint='$R/epoch=377-step=171990.ckpt'" regions.source=split regions.split=val \
    data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/mm10.fa \
    "data.cfg.h5_path=[/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions.h5]"
pixi run python scripts/metrics.py hydra.run.dir=$R/metrics_val_ac_epoch377_plots \
    samples_path=$R/samples_val_ac_epoch377/samples.h5 plot.enabled=true
```

Comparison: `pixi run -e cpu python scripts/compare.py`, using
[scripts/compare.py](frozen-alphagenome-encoder/scripts/compare.py). Its output
is in [comparison.txt](frozen-alphagenome-encoder/comparison.txt). The step
benchmark is in
[scripts/bench.py](frozen-alphagenome-encoder/scripts/bench.py).

Models compared:

- Baseline: `smf_models/SMFNet_full_run_2026-08-11/epoch=377-step=171990.ckpt`,
  trained with an older repo version but the same config. Its metrics are in
  `/scratch/users/diamant/evaluate_SMFNet_full_run_2026-08-11/metrics_val_ac_current_code_plots`.
- ATAC: `smf_models/ONT_mESC_ATAC_TKO_pooled`, with metrics in
  `metrics_val_ac_epoch377_plots`.

## Results (val split, active regions, 7,559 regions)

| | baseline | AlphaGenome frozen | ATAC |
| --- | --- | --- | --- |
| Best val/loss (ac + bg) | 0.29436 | 0.29373 | 0.29328 |
| `pearson_matched` (fraction of gap) | 0.731 (0.888) | 0.760 (0.925) | 0.772 (0.941) |
| `pearson_flank_residual_matched` | 0.702 (0.850) | 0.736 (0.900) | 0.747 (0.916) |
| `run_length_zeros_wasserstein` (bp) | 3.32 | 3.16 | 3.11 |
| `run_length_ones_wasserstein` (bp) | 7.07 | 7.19 | 7.43 |
| Pooled site Pearson (all reads) | 0.866 | 0.904 | 0.911 |

The "fraction" is the share of the gap between the null and the floor that the
model closes, as in `metrics_summary.tsv`.

Paired per-region differences from the baseline (mean ± SE, and the share of
regions where the model is better):

| metric | AlphaGenome | ATAC |
| --- | --- | --- |
| `pearson_matched` | +0.028 ± 0.002 (61%) | +0.036 ± 0.002 (62%) |
| `pearson_flank_residual_matched` | +0.029 ± 0.002 (61%) | +0.039 ± 0.002 (62%) |
| `covariation_corr_of_corr_pearson` | +0.002 ± 0.0005 (52%) | +0.004 ± 0.0005 (54%) |
| `covariation_residual_corr_of_corr_pearson` | +0.003 ± 0.001 (51%) | +0.008 ± 0.001 (53%) |
| `run_length_zeros_wasserstein` | −0.33 ± 0.03 bp (54%) | −0.46 ± 0.03 bp (56%) |
| `run_length_ones_wasserstein` | +0.13 ± 0.03 bp (48%) | +0.31 ± 0.03 bp (45%) |

## Caveats

- One seed per model. The per-region SEs capture variation across regions, not
  across training runs.
- (Resolved, see the update below.)
- The baseline samples were produced by older sampling code (sample-file
  schema v1; the AlphaGenome and ATAC samples are v2). The samples themselves
  are unchanged; only a `data_config` attribute was added on 2026-10-27. A
  sampling-code difference could explain part of the gap.
- The baseline checkpoint's weight names carry `model._orig_mod.` (it was
  compiled the old way) and don't load with the current code. A
  prefix-stripped copy is at
  `SMFNet_full_run_2026-08-11/epoch=377-step=171990-converted.ckpt` and loads
  with `load_model`.

## Update: baseline resampled with current code (2026-10-29)

Job 50155720 (btrippe, 1x H100, COMPLETED in 1 h 38 m, MaxRSS 10 GB) ran
[scripts/sample_metrics_baseline_val.sbatch](frozen-alphagenome-encoder/scripts/sample_metrics_baseline_val.sbatch):

```bash
sbatch $SCRATCH/smfnet_ag_smoke/sample_metrics_baseline_val.sbatch
```

It resampled the baseline (`epoch=377-step=171990-converted.ckpt`) on the same
val active regions and settings as above, then ran metrics with plots. Outputs:
`SMFNet_full_run_2026-08-11/samples_val_ac_epoch377_current_code` and
`metrics_val_ac_epoch377_current_code_plots`. The comparison used
[scripts/compare_resampled.py](frozen-alphagenome-encoder/scripts/compare_resampled.py),
with output in [compare_resampled.txt](frozen-alphagenome-encoder/compare_resampled.txt).

- **Old vs new baseline samples** agree within sampling noise (paired, mean ±
  SE): `pearson_matched` +0.0012 ± 0.0012 (51% of regions), zeros run length
  +0.02 ± 0.02 bp, ones run length −0.007 ± 0.03 bp, pooled site Pearson
  0.86616 vs 0.86621.
- **Against the resampled baseline**, the gains are unchanged:

| metric | AlphaGenome | ATAC |
| --- | --- | --- |
| `pearson_matched` | +0.029 ± 0.001 (62%) | +0.037 ± 0.001 (64%) |
| `pearson_flank_residual_matched` | +0.031 ± 0.001 (62%) | +0.040 ± 0.002 (65%) |
| `covariation_corr_of_corr_pearson` | +0.003 ± 0.0004 (53%) | +0.005 ± 0.0004 (55%) |
| `covariation_residual_corr_of_corr_pearson` | +0.005 ± 0.0007 (52%) | +0.010 ± 0.0008 (55%) |
| `run_length_zeros_wasserstein` | −0.31 ± 0.02 bp (55%) | −0.44 ± 0.03 bp (56%) |
| `run_length_ones_wasserstein` | +0.12 ± 0.02 bp (47%) | +0.31 ± 0.02 bp (43%) |

The resampled baseline's `pearson_matched` is 0.7296 (old samples: 0.7310).
AlphaGenome closes 71% of the gap between the baseline and ATAC using the
summary values ((0.7601 − 0.7296) / (0.7723 − 0.7296)), or 79% using the paired
per-region means (0.0293 / 0.0369).

# Related Notes

- [Flank baseline for profile metrics](../2026-10-05/flank-baseline-logistic-regression.md):
  defines the flank baseline behind `pearson_flank_residual_matched`, the
  metric showing that the gain goes beyond local sequence composition.
- [Positional baselines and noise-matched corr-of-corr](../2026-09-29/positional-covariation-baselines-and-noise-matched-corr-of-corr.md):
  background on the covariation metrics, which move only slightly here.
- [What ends observed SMF runs](../2026-10-08/run-length-missingness-gap-bridging.md):
  background on the run-length Wasserstein metrics, including the ones run
  lengths that got worse.

# Open Questions

- (Resolved) Does the gain hold when the baseline is resampled with the
  current sampling code? Yes; see the update below.
- Are the AlphaGenome and ATAC gains additive? Train with both.
- Does unfreezing after 1 epoch, or a trainable dilated body on top, improve
  further, or does it overfit?
- Is the improvement concentrated at high mask rates, as inferred? Check with
  `log_mask_rate_bins` or post hoc binned val loss.
- Is it concentrated in the CpG (endogenous methylation) or GpC
  (accessibility) channel? AlphaGenome likely has no methylation training
  targets.
- Why do ones run lengths get worse for both conditioned models?
- Test-split results are reserved for the selected configuration.

# Sources

- Blog post: https://genomicsxai.github.io/blogs/2026-002/
- alphagenome-pytorch source (commit 72268c0): https://github.com/genomicsxai/alphagenome-pytorch
- Weights: https://huggingface.co/gtca/alphagenome_pytorch (the license file
  points to https://deepmind.google.com/science/alphagenome/model-terms)
- SMFNet_public branch `nd_alpha_genome`, commit 803eb27, and
  `ALPHAGENOME_IMPLEMENTATION_PLAN.md` (untracked, in that checkout)
- [comparison.txt](frozen-alphagenome-encoder/comparison.txt),
  [compare.py](frozen-alphagenome-encoder/scripts/compare.py),
  [bench.py](frozen-alphagenome-encoder/scripts/bench.py)
