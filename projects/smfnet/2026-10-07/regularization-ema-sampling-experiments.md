---
title: GCH overfitting, regularization, EMA, and sampling-length experiments
date: 2026-10-07
project: smfnet
agent: Claude Code
status: draft
sources:
  - regularization-ema-sampling-experiments/RUN_LOG.md (full timeline of the proxy runs)
  - regularization-ema-sampling-experiments/summary.tsv, loss_curves.png, param_norms.png
  - regularization-ema-sampling-experiments/windowed_sampling_comparison.tsv
  - regularization-ema-sampling-experiments/scripts/ (run_one.sh, run_queue.sh, summarize.py, compare_windows.py)
  - wandb runs jnk4c93j (ES), 6dwjbqab (ES_jitter128), 3ekr0isw (ONT full model), project diamant-stanford-university/SMFNet
  - SMFNet repo, scripts/train.py and src/smf_net/model.py (branch nd/more_assays, uncommitted at time of writing)
  - User discussion, 2026-10-06 to 2026-10-07
tags: [smfnet, gch, overfitting, dropout, weight-decay, ema, jitter, diffusion, sampling, checkpoint-selection, ont]
---

# Summary

The GCH models overfit early: training loss keeps falling while validation
loss flattens and creeps up after epoch ~30–70. On a small-data proxy (ES, 5,000
training regions), **EMA of the weights helped most** (validation −0.0045 vs
the baseline), **dropout 0.2 helped less** (−0.0025) but closed most of the
train/validation gap, and **together they matched the gain from tripling the
training data** (−0.0054, gap 0.024 → 0.007). Weight decay 0.1, a half-width
model and positional jitter (64 bp in the proxy; 128 bp on full data) did
little or nothing. Sampling only a short window, to match the short training
reads, did not change scores (dense full-window sampling is not the problem).
A full-data ES run with dropout 0.2 + EMA 0.999 (`ES_dropout20_ema999`) was
started, with early stopping still on validation loss.

This follows the finding in
[the covariation diagnosis](../2026-10-06/covariation-differences-diagnostics.md)
that the GCH models fall well short of their noise ceilings.

# Key Points

- **The proxy works:** ES capped at 5,000 training regions, validated on the
  same 15,000 regions in every run, ~90 epochs in 45–60 min on an H100. Two seeds
  of the baseline agree to 0.0004 in mean validation loss over epochs 70–89.
- **Where memorizing happens:** the DNA encoder's weights keep growing (~2x by
  epoch 80, 3x with 15,000 regions) long after validation stops improving. The
  rest of the network barely changes. DNA-encoder-only dropout matches dropout
  everywhere.
- **Overfitting is real but not the main limit:** dropout closes most of the gap
  but lowers validation loss by only ~40% of what tripling the data does.
- **Jitter can't help this architecture much:** the network is fully
  convolutional and close to shift-equivariant, so shifting the input mostly
  shifts the output. Full-data `ES_jitter128` tracked the original ES run within
  0.0006 at every epoch window.
- **Validation loss is a flat average over masking levels:** with the cosine
  schedule's weighting, each mask-rate bin contributes equally. (An earlier claim
  in the discussion that heavily, then lightly, masked reads dominate was wrong
  both times.)
- **Sampled scores keep moving after the validation minimum, but noisily:**
  covariation (correlation of correlations) kept improving after the minimum in
  one EMA run. Repeat-seed noise of 0.005–0.008 per evaluation (256 regions ×
  256 reads) makes most later differences 1–2x noise, so the decision was to
  keep early stopping on validation loss and use sampled scores as a monitor.
- **ONT is a different regime:** the ONT model's validation loss improves to
  epoch ~377, its gap is only 0.009 at epoch 400 and its noise is ~4x lower than
  GCH ES. Dropout is unlikely to help it, and EMA likely adds little at a fully
  annealed checkpoint.

# Details

## Proxy results

All runs: ES, 5,000 training regions, same 15,000 validation regions. Validation
and gap are means over epochs 70–89 (`summary.tsv`, `loss_curves.png`).

| Run | Change | Validation | vs baseline (0.4500) | Gap |
|---|---|---|---|---|
| `drop02_ema098` | dropout 0.2 + EMA 0.98 | 0.4446 | −0.0054 | 0.007 |
| `ema098` | EMA 0.98 | 0.4455 | −0.0045 | 0.019 |
| `dropout02` | dropout 0.2 everywhere | 0.4475 | −0.0025 | 0.010 |
| `dropout03` | dropout 0.3 everywhere | 0.4479 | −0.0021 | 0.007 |
| `dnadrop03` | dropout 0.3, DNA encoder only | 0.4480 | −0.0020 | 0.012 |
| `dropout01` | dropout 0.1 everywhere | 0.4484 | −0.0016 | 0.015 |
| `dnadrop05` | dropout 0.5, DNA encoder only | 0.4486 | −0.0014 | 0.008 |
| `wd01` | weight decay 0.1 | 0.4488 | −0.0012 | 0.023 |
| `base5k_s0` / `_s1` | none, two seeds | 0.4497 / 0.4502 | | 0.024 |
| `jitter64` | 64-bp jitter | 0.4503* | +0.0003 | 0.024 |
| `narrow` | half width | 0.4507 | +0.0007 | 0.023 |

\* scored on a 992-bp window, so only roughly comparable. With 15,000 training
regions the baseline reached ~0.443, so tripling the data gives about −0.0054.

EMA decay 0.98 in the proxy (20 steps per epoch) stands in for 0.999 on full
data (~615 steps per epoch); both average over ~1.6–2.5 epochs.

Weight decay barely matters because AdamW decay is scaled by the learning rate:
0.1 × 1e-3 is 1e-4 per step. Early stopping ended `narrow` at epoch 151, so
runs are compared on a matched epoch window.

## Sampled scores during training (EMA runs)

Medians over 256 fixed validation regions × 256 reads, every 10 epochs.
Replicate ceilings: profile 0.726, covariation 0.739.

| Epoch | Profile, EMA | Profile, dropout + EMA | CoC, EMA | CoC, dropout + EMA |
|---|---|---|---|---|
| 49 | 0.446 | 0.408 | 0.625 | 0.618 |
| 59 | 0.447 | 0.430 | 0.626 | 0.629 |
| 69 | 0.419 | 0.461 | 0.633 | 0.622 |
| 89 | 0.438 | 0.440 | 0.631 | 0.629 |
| 99 | 0.433 | 0.458 | 0.639 | 0.632 |

In `ema098` the profile dip after epoch 59 turned out to be churn, not a
uniform decline: the 25th-percentile region *r* rose while the median fell.
Predicted profile spread plateaued at ~2/3 of observed. Mean bias drifted from
+0.8 to ~+2.5 pp late in training. The spread of predicted pair correlations
kept growing toward the observed level (dependence still being learned). The
non-EMA runs with sampling evaluation (`base5k_sampled`, `dropout02_sampled`)
were stopped to free the GPU for the full-data run, so the effect of EMA alone
on sampled scores is not isolated.

## Sampling only the training read length

ES model, 2,000 validation regions, 1,024 reads, same seed. Dense (whole
1,024-bp window per read) vs windowed (`sampling.sampling_window=[384,640]`),
both scored on the central 256 bp (`windowed_sampling_comparison.tsv`).
Windowed minus dense, with paired region-bootstrap 95% CIs: profile *r* +0.006
pooled and +0.009 within regions, covariation −0.001. So generating reads far
longer than the ~165-bp training reads does not distort site-level or
short-range predictions.

## Full-data comparisons (wandb)

- **ES vs ES_jitter128:** single best epochs 33 vs 74 are noise on a flat
  plateau (epoch-to-epoch SD 0.0012). Best 5-epoch means are 0.4323 (epoch 35)
  and 0.4321 (epoch 45). Window means are within 0.0006, and gaps grow
  identically.
- **ONT (3ekr0isw):** validation 0.3017 (epochs 20–39) → 0.2953 (380–399), gap
  +0.003 → +0.009, noise 0.0003–0.0004. That run did not log the learning rate,
  so the annealing argument assumes it followed its 400-epoch cosine config.

## Decisions and code

- Full-data runs use dropout 0.2 + EMA 0.999, early stopping on validation loss
  (patience 100), validation loss by mask rate, and checkpoints plus sampling
  evaluation every 10 epochs. `ES_dropout20_ema999` was started by the user;
  `scripts/train_gch_dropout_ema.sbatch` does the same for MEL, NP and C2C12.
- Kept in SMFNet (`scripts/train.py`, `src/smf_net/model.py`,
  `src/smf_net/evaluation/sampling.py`): `ema_decay`, `log_mask_rate_bins`,
  `checkpoint_every_n_epochs`, `sample_eval_every_n_epochs`, the CSV logger and
  `sampling.sampling_window`.
- Removed after the experiments: `data.cfg.max_train_regions_per_h5`,
  `model.cfg.cnn_config.dna_dropout` and `log_layer_norms`. The copied
  `scripts/run_one.sh` uses the first and last of these, so it reproduces the
  runs only against the 2026-10-07 code, not the current branch.

# Related Notes

- [Covariation differences diagnosis](../2026-10-06/covariation-differences-diagnostics.md): Showed the models fall short of their ceilings and predict cell types too alike, which motivated these experiments.
- [Differential evaluation with held-out molecules](../2026-10-06/differential-evaluation-held-out-molecules.md): The evaluation the new models will be scored with.

# Open Questions

- Do the proxy gains carry over to full data? Compare `ES_dropout20_ema999`
  with the original ES run (wandb `jnk4c93j`) on validation loss, sampled
  scores and the differential strata.
- Does EMA alone change sampled scores? This needs the unfinished non-EMA runs
  with sampling evaluation.
- Would one model trained on all four cell types, with shared DNA encoder and a
  cell-type input, reduce memorization more than any regularizer here?

# Sources

- [Run log](regularization-ema-sampling-experiments/RUN_LOG.md), [summary table](regularization-ema-sampling-experiments/summary.tsv), [loss curves](regularization-ema-sampling-experiments/loss_curves.png), [parameter norms](regularization-ema-sampling-experiments/param_norms.png), [windowed sampling comparison](regularization-ema-sampling-experiments/windowed_sampling_comparison.tsv)
- Scripts: [run_one.sh](regularization-ema-sampling-experiments/scripts/run_one.sh), [run_queue.sh](regularization-ema-sampling-experiments/scripts/run_queue.sh), [summarize.py](regularization-ema-sampling-experiments/scripts/summarize.py), [compare_windows.py](regularization-ema-sampling-experiments/scripts/compare_windows.py)
- Run outputs (purged after 90 days unmodified): `/scratch/users/diamant/smf_models/overnight_20261006/`, `/scratch/users/diamant/windowed_sampling_20261006/`
