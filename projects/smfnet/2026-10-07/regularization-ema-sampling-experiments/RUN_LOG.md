# Overnight regularization runs, 2026-10-06/07

Goal: compare ways to reduce early overfitting of the GCH models (train loss
keeps falling while validation loss rises after epoch ~30–70), using fast
small-data proxies. Run by Claude Code in Slurm job 46855199 (sh04-03n07,
1 H100, ends 2026-10-07 08:40).

## Setup

- Model/data: ES, `gch_accessibility`, same settings as the production ES run
  (`smf_models/ES`), RC augmentation off.
- Proxy: `data.cfg.max_regions_per_h5_per_split=15000` (a fixed, seeded subset
  of about 10% of training regions; validation is capped at 15,000 too, the same
  regions in every run). All reads kept.
- Length: `trainer.max_time=1h20m` wall time per run. `max_epochs` stays 400
  so the cosine learning-rate schedule matches production.
- Logging: wandb offline (nothing uploaded; runs in
  `$SCRATCH/smf_models/overnight_20261006/wandb`), plus a local CSV log of
  every metric, plus per-layer-group gradient norms (every 50 steps) and
  parameter norms (each epoch).
- Runner: `analyses/overnight_regularization_20261006/run_one.sh NAME [overrides]`.
  Outputs: `$SCRATCH/smf_models/overnight_20261006/<NAME>/` (checkpoints,
  `csv/metrics.csv`), logs in `.../logs/`.

## Code changes (all opt-in; revert by removing)

| File | Change | Default behavior |
|---|---|---|
| `scripts/train.py` | CSV logger next to wandb (`<run dir>/csv/metrics.csv`) | adds a local file only |
| `scripts/train.py` | `log_layer_norms` flag and `LayerNormMonitor` callback | off |
| `src/smf_net/cnn.py` | `CNNConfig.dna_dropout` (dropout inside the DNA encoder only) | `None` = use `dropout`, as before |

## Planned runs

1. `baseline_s0`: no changes, seed 0
2. `baseline_s1`: seed 1, to measure run-to-run noise
3. `jitter64`: `data.cfg.model_input_length=1984 model.cfg.input_length=1984`
4. `dropout02`: `model.cfg.cnn_config.dropout=0.2` (all blocks)
5. `wd01`: `model.cfg.wd=0.1`
6. `dnadrop03`: `model.cfg.cnn_config.dna_dropout=0.3`
7. `narrow`: `model.cfg.cnn_config.n_filters=64 model.cfg.cnn_config.dna_encoder_filters=256`
8. If time remains: `jitter128` or read thinning 0.1

## Log

- 22:45 Smoke test (512 regions, 3 min) passed: time limit, CSV log and layer
  norms all work. Even 512 regions overfit within ~55 epochs.
- 22:53 Started `baseline_s0` (59 steps per epoch, ~42 s per training epoch, ~1 min per epoch with validation).
- 22:58 GPU: the baseline alone uses 74 GB and 100% of the H100, so runs go one at a time via `run_queue.sh` (queue file `$SCRATCH/smf_models/overnight_20261006/queue.txt`, editable while running). Order: dropout02, dnadrop03, wd01, baseline_s1, narrow, jitter64. jitter64 is last because a full-data `ES_jitter128` run (started by the user) is already going in another job.
- 23:30 `baseline_s0` at epoch 35: validation loss 0.4424 (best, epoch 30), training
  0.435, gap 0.006. Not overfitting yet; validation is still flat to slightly falling.
  Epoch-to-epoch validation noise is about ±0.002–0.003, so differences between
  runs smaller than about 0.003 won't be meaningful without the second seed. If the
  baseline hasn't overfit by the end (~75 epochs), the remaining runs will use a
  smaller region cap so they overfit within their time.
- 00:12 `baseline_s0` (15,000 training regions) finished 77 epochs: validation
  plateaued at ~0.443 from epoch 30 on (best 0.4403 at epoch 50, within noise),
  while training loss kept falling 0.436 → 0.425 (gap 0.006 → 0.018). Mild
  overfitting, but no clear validation rise in 80 min, so this cap is too
  large for comparing regularizers in the time available.
- 00:13 Change of plan: remaining runs train on 5,000 regions and validate on the
  same 15,000 regions as `baseline_s0` (new opt-in `data.cfg.max_train_regions_per_h5`
  in `src/smf_net/data.py`), 45 min each (~75 epochs). A 5,000-region baseline
  (`base5k_s0`, `base5k_s1`) replaces the 15,000-region one for comparisons.
  Queue: base5k_s0, dropout02, dnadrop03, wd01, base5k_s1, narrow, jitter64,
  dnadrop05, wd03. Also fixed `run_one.sh` so a failed run still logs its end.
  The queue runner keeps its startup cutoff of 85 min left in the job, so the
  last run starts by ~07:15.
- 00:17 Started `base5k_s0` (20 steps per epoch, ~35 s per epoch with validation).
- 01:07 `base5k_s0` (5,000 training regions) finished 93 epochs: validation plateaued
  at ~0.449 from epoch ~45 on (last-20-epoch mean 0.4496), training kept falling to
  ~0.423, gap 0.024 (vs 0.018 at 15,000 regions). So 5,000 regions overfit more,
  but validation still plateaus rather than rising, probably because the cosine
  learning rate (over 400 epochs) is still near its peak. Same schedule for every
  run, so comparisons are fair. Comparison measures added to `summarize.py`:
  mean validation loss over the last 20 epochs (`val_last20`), lowest 5-epoch
  rolling mean, and the last-20 train/validation gap.
- 01:07 Parameter norms: the DNA encoder's weights keep growing steadily (2.3x the
  first-epoch norm at 5,000 regions, 3.1x at 15,000) long after validation stops
  improving; read blocks, DNA conditioning and the input layer stay at 1.0–1.25x.
  The output layer grows ~1.5–2x, and jumped abruptly at epoch 50 in `baseline_s0`.
  Consistent with memorization in the DNA encoder, which favors `dna_dropout`
  and stronger weight decay (both queued).
- 01:04 Started `dropout02`.
- 01:59 `dropout02` (dropout 0.2 in every block) finished 88 epochs: last-20
  validation 0.4475 vs 0.4496 for `base5k_s0` (−0.002, about the size of the
  epoch-to-epoch noise; needs `base5k_s1` to judge). The gap shrank from 0.024 to
  0.010, mostly because training loss stayed higher (0.437 vs 0.421) rather than
  because validation fell. DNA-encoder weight growth slowed (1.77x vs 2.17x at
  epoch 85); other parts similar. Reading so far: dropout stops much of the
  memorization but buys little generalization at this data size.
- 01:51 Started `dnadrop03`.
- 02:50 `dnadrop03` (dropout 0.3 in the DNA encoder only) finished 91 epochs:
  last-20 validation 0.4480 (−0.0016 vs `base5k_s0`), gap 0.0125, DNA-encoder
  growth 1.65x at epoch 80 (baseline 2.09x, `dropout02` 1.68x). DNA-only dropout
  gives nearly the same effect as dropout everywhere, consistent with the
  memorizing being concentrated in the DNA encoder. Both gains (~0.002) are at
  the noise level until `base5k_s1` finishes.
- 02:50 Timing: with the queue's 85-min cutoff, the last slot is probably lost, so
  `wd03` (last) likely won't run. Will choose between it and `dnadrop05` after `wd01`.
- 02:38 Started `wd01`.
- 03:31 `wd01` (weight decay 0.1, 10x the default) finished 93 epochs: last-20
  validation 0.4486 (−0.001 vs `base5k_s0`), gap 0.023 (unchanged), DNA-encoder
  growth still 1.93x. Little effect: with AdamW the decay per step is lr × wd =
  1e-4, too weak to stop the growth. Dropped `wd03` from the queue.
- 03:31 New queue after `base5k_s1`: narrow, dnadrop05, jitter64, drop02_dna05
  (dropout 0.2 in the read blocks plus 0.5 in the DNA encoder). Four slots fit
  before the 85-min cutoff.
- 03:24 Started `base5k_s1`.
- 04:17 `base5k_s1` (seed 1) finished: last-20 validation 0.4500 vs 0.4496 for seed 0
  (difference 0.0004; within-run standard error of the last-20 mean ~0.0003), gap
  0.0245 vs 0.0242. Against the mean of the two baselines (0.4498): `dropout02`
  −0.0023, `dnadrop03` −0.0018 (both ~5x the seed difference, so real), `wd01`
  −0.0012 (borderline). For scale, tripling the training regions (5,000 → 15,000)
  lowers validation by 0.0054, so dropout recovers about 40% of that.
- 04:11 Started `narrow`.
- 05:03 `narrow` (n_filters 64, DNA encoder 256; half width) ran faster and early
  stopping ended it at epoch 151 (best at epoch 50 + patience 100). Its last-20
  window is therefore epochs 131–150, so I added a matched window (epochs 70–89)
  to `summarize.py`. Matched: validation 0.4507 vs 0.4500 for the two baselines,
  gap 0.023 vs 0.024. Halving the width does not help.
  Matched epochs 70–89, validation / gap: base5k_s0 0.4497 / 0.024, base5k_s1
  0.4502 / 0.024, dropout02 0.4475 / 0.010, dnadrop03 0.4480 / 0.012, wd01
  0.4488 / 0.023, narrow 0.4507 / 0.023.
- 04:49 Started `dnadrop05`.
- 05:47 `dnadrop05` (DNA-encoder dropout 0.5) finished: epochs 70–89 validation
  0.4486, gap 0.008, DNA-encoder growth 1.54x. Versus `dnadrop03` (0.4480, gap
  0.012): stronger DNA dropout narrows the gap further but no longer improves
  validation, so the benefit peaks around 0.3. Replaced `drop02_dna05` (likely
  over-regularized) with `dropout01` and `dropout03` (dropout everywhere at 0.1 and
  0.3), which map the curve around the best run so far, `dropout02`. Both fit
  before the cutoff. `summarize.py` now accepts a 70–89 window with >=15 epochs
  (`dropout02` stopped at 88).
- 05:36 Started `jitter64`.
- 06:31 `jitter64` (input 1,984 bp, up to 64 bp of shift) finished 94 epochs: epochs
  70–89 validation 0.4503 (baselines 0.4497/0.4502; scored on a 992-bp window, so
  only roughly comparable), gap 0.024 (unchanged), DNA-encoder growth 2.1x
  (unchanged). No effect. Likely reason: the network is fully convolutional, so
  it is close to shift-equivariant. Shifting the input mostly shifts the output;
  only the window edges see different context, so the model can still memorize
  each region's sequence. The user's full-data `ES_jitter128` run will probably
  show the same, though 128 bp changes the edges a bit more.
- 06:23 Started `dropout01`.
- 07:17 `dropout01` (dropout 0.1 everywhere) finished 91 epochs: epochs 70–89
  validation 0.4484, gap 0.015. So far dropout everywhere: 0.1 → 0.4484 (gap 0.015),
  0.2 → 0.4475 (gap 0.010); `dropout03` will show whether 0.3 helps further.
- 07:09 Started `dropout03`, the last run; ends ~07:56.
- 07:56 `dropout03` (dropout 0.3 everywhere) finished 90 epochs: epochs 70–89
  validation 0.4479, gap 0.007. Queue empty at 07:56; the final wrap-up check did
  not run before the job ended, so these results were collected the next morning.

## Final results (5,000 training regions, same 15,000 validation regions)

Validation and gap are means over epochs 70–89. DNA-encoder growth is the
parameter norm at epoch 80 relative to epoch 0.

| Run | Change | Validation | vs baseline mean (0.4500) | Gap | DNA-encoder growth |
|---|---|---|---|---|---|
| `dropout02` | dropout 0.2 everywhere | 0.4475 | −0.0025 | 0.010 | 1.68x |
| `dropout03` | dropout 0.3 everywhere | 0.4479 | −0.0021 | 0.007 | 1.72x |
| `dnadrop03` | dropout 0.3, DNA encoder only | 0.4480 | −0.0020 | 0.012 | 1.65x |
| `dropout01` | dropout 0.1 everywhere | 0.4484 | −0.0016 | 0.015 | 1.64x |
| `dnadrop05` | dropout 0.5, DNA encoder only | 0.4486 | −0.0014 | 0.008 | 1.54x |
| `wd01` | weight decay 0.1 (10x) | 0.4488 | −0.0012 | 0.023 | 1.93x |
| `base5k_s0` | none, seed 0 | 0.4497 | | 0.024 | 2.09x |
| `base5k_s1` | none, seed 1 | 0.4502 | | 0.024 | 1.94x |
| `jitter64` | 64-bp jitter | 0.4503* | +0.0003 | 0.024 | 2.10x |
| `narrow` | half width | 0.4507 | +0.0007 | 0.023 | 1.63x |

\* scored on a 992-bp window, so only roughly comparable.
Reference: `baseline_s0` with 15,000 training regions has best 5-epoch validation
0.4428 (it stopped before epoch 89), so tripling the data improves validation by
about 0.005–0.006.

## Conclusions

1. Dropout is the only lever that clearly helps: 0.2–0.3 everywhere, or 0.3 in
   the DNA encoder alone, all lower validation by ~0.002 (5x the seed difference
   of 0.0004). The differences among these settings (≤0.001) are within ~2–3x the
   seed difference, so treat them as a tie; dropout 0.2 everywhere is the simplest
   choice.
2. Dropout closes most of the train/validation gap (0.024 → 0.007–0.012) but
   lowers validation only modestly: about 40% of what tripling the training data
   does. So overfitting is real but is not the main limit; the amount of
   information in the training data is.
3. Weight decay 0.1 (decoupled, lr × wd = 1e-4 per step), half width, and 64-bp
   jitter do little or nothing. Jitter is weak because a fully convolutional
   network is close to shift-equivariant.
4. The DNA encoder is where weights keep growing after validation stops
   improving (2x in 80 epochs; 3x with 15,000 regions). Dropout slows it, and DNA
   dropout alone matches dropout everywhere, so the memorizing is concentrated
   there.

Caveats: these are small-data proxies (5,000 training regions, ~90 epochs, the
learning rate still near its peak), and lower validation loss does not guarantee
better profile or differential scores. Next step: a full-data ES run with
`model.cfg.cnn_config.dropout=0.2`, evaluated with `metrics.py` and the
differential pipeline.

## EMA follow-up, 2026-10-07 (job 46675939, same node, H100)

Question: does EMA help, alone or with dropout, and do sampled metrics keep
improving after the validation minimum? Four runs, same proxy as above (5,000
training regions, same 15,000 validation regions) but 60 min each, so the new
sampling evaluation (256 regions x 256 reads every 10 epochs) still leaves ~90
epochs. All log validation cross-entropy by mask rate.

EMA decay 0.98 (horizon ~50 steps = ~2.5 proxy epochs) stands in for 0.999 on
full data (~1,000 steps = ~1.6 epochs at 615 steps per epoch). With 0.999 the
proxy's average would lag ~50 epochs and still hold ~16% of the initial
weights at the end.

Queue (`queue_ema.txt`): ema098, drop02_ema098, base5k_sampled, dropout02_sampled.
The last two repeat last night's baseline and dropout 0.2 with sampling
evaluation, so sampled scores can be compared across all four.

- 10:16 Started `ema098`.
- 11:21 `ema098` (EMA decay 0.98) finished 103 epochs. Epochs 70–89 validation
  0.4455 (−0.0045 vs the baseline mean 0.4500, ~11x the seed difference; dropout
  0.2 gave −0.0025, tripling the data −0.0054), gap 0.019 (baseline 0.024,
  dropout 0.010). EMA lowers validation more than dropout, while dropout closes
  more of the gap.
  Sampled scores (256 validation regions x 256 reads, medians), with sampling
  noise from a repeat seed at epoch 9 (profile 0.005, covariation 0.003; the
  matched profile Pearson is much noisier, 0.038, so ignored):

  | Epoch | Profile Pearson | Corr. of corr. |
  |---|---|---|
  | 9 | 0.115 | 0.517 |
  | 29 | 0.310 | 0.602 |
  | 49 | 0.446 | 0.625 |
  | 59 | 0.447 | 0.626 |
  | 69 | 0.419 | 0.633 |
  | 89 | 0.438 | 0.631 |
  | 99 | 0.433 | 0.639 |

  Replicate ceilings 0.726 (profile) and 0.739 (covariation). Covariation keeps
  improving after the validation minimum (~epoch 50): +0.014 by epoch 99, ~5x
  noise. The profile peaks at epochs 49–59 and slips ~0.01–0.015 afterwards. So
  the best-by-loss checkpoint suits profiles, while covariation favors longer
  training.
  Validation cross-entropy by mask rate: every bin improves until epoch ~40–60,
  then stays flat; the 0.4–0.8 bins creep up ~0.001–0.002. Heavily masked reads
  (0.8–1.0) are hardest (0.625 vs 0.347 for 0–0.2). No masking level keeps
  falling, even though sampled covariation keeps improving.
  Sampling evaluation took 48 s each (99 s the first time, with the repeat seed).
- 11:18 Started `drop02_ema098`.
- 12:23 `drop02_ema098` (dropout 0.2 + EMA 0.98) finished 101 epochs: epochs 70–89
  validation 0.4446 (−0.0054 vs baseline mean; EMA −0.0045 and dropout −0.0025
  combine sub-additively, matching the gain from tripling the data), gap 0.007,
  best epoch 79 (later than EMA alone, 50). Sampled profile Pearson keeps rising
  to ~0.46 by epoch 69 and stays there; covariation levels off at ~0.63. So the
  profile decline in `ema098` was likely noise: repeat-seed noise here is 0.008
  for both metrics (vs 0.005/0.003 in `ema098`), so sampled-score differences
  of ~0.01–0.02 after epoch 50 are only 1–2x noise.
  Decision (with the user): keep early stopping on validation loss; sampled
  scores are a monitor only. Full-data run started by the user:
  `ES_dropout20_ema999` (dropout 0.2, EMA 0.999, mask-rate bins, checkpoints
  and sampling evaluation every 10 epochs).
- 12:42 Stopped the queue: `base5k_sampled` was cut short at ~20 min and
  `dropout02_sampled` never ran, to free the shared H100 for `ES_dropout20_ema999`
  in the same job. So there is no non-EMA run with sampling evaluation yet.
