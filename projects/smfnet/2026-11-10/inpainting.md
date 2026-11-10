---
title: Span-inpainting metric; the U-Net's long-gap weakness and its drift over training; inpainting-variant sweep; k=9 + dilated fine convs as the new best (val loss 0.29229)
date: 2026-11-10
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - SMFNet_public branch nd_unet_inpaint (worktree ~/repos/SMFNet_inpaint), commits 8833c79, ec73658; runs from the detached checkout ~/repos/SMFNet_inpaint_runs at ec73658
  - inpainting/scripts/ (scorer, runners, summaries, outputs)
tags: [evaluation, inpainting, unet, dilation, normalization, attention, sampling, muon, ema]
---

# Summary

**Why a new metric.** The k=9 attention U-Net improved val loss
(0.29314 → 0.29259), but its sampled-read metrics did not move. All models
share the same short-range errors in their sampled reads, from the original
SMFNet to the U-Nets:

- sampled reads are too weakly correlated at 0–50 bp, by ~0.04;
- they have too many short accessible runs.

Doubling the sampling steps shrinks these errors by ~25%. That suggests the
parallel-decoding sampler, not the models, limits fine-scale structure in
the samples. Confidence-ordered sampling was rejected earlier because it
distorts the distribution.

**The metric.** A sampling-free metric that sees fine-scale structure
directly. One contiguous span (20, 50 or 150 bp) of each held-out read is
hidden and scored, given the rest of the read and the DNA:

- **joint:** the exact chain-rule likelihood of the span's true calls,
  revealed one site at a time in random order. The main inpainting number.
- **one-shot:** all of the span in one pass. A diagnostic of how far the
  model generalizes beyond its training masks; the sampler never visits this
  state.

**Findings.**

1. **The dilated-CNN-body models inpaint long spans better than every
   U-Net,** despite worse val loss. Frozen AlphaGenome against the k=5
   attention U-Net at 150 bp: joint −0.0064, one-shot −0.034 nats per call.
2. **The U-Net's one-shot inpainting gets steadily worse over training**
   while val loss improves. This happens for the attention U-Net with the
   per-read norm, not for the plain U-Net (per-position norm, no attention).
   The drift looks tied to read-level mechanisms that reward specializing in
   scattered masks (inferred), not to capacity.
3. **50-epoch architecture sweep** (linear upsampling, dilated fine convs,
   local attention, a 16 bp bottleneck, a mask-aware read norm): only
   **dilated fine convs** helped. At epoch 49 they improved 150 bp joint by
   −0.0020 ± 0.0002 against a run-to-run difference of ~0.0003. The
   mask-aware norm was worse.
4. **Full run, k=9 + dilated fine convs (dilations 4, 4, 2 at the three
   finest levels):**
   - **new best val loss, 0.29229** (k=9: 0.29259), at no extra compute;
   - a modest 32-step sampled-metric gain over k=9;
   - its inpainting advantage disappeared by epoch 399, because the drift
     continued.

   Adopt it as the default. Long-gap inpainting would need training-side
   changes, such as span masks, if it ever matters. Inpainting is a metric
   here, not an objective, so that's not planned.

# Key Points

- **The sampler limits fine-scale structure in samples.** Pooled over 7,559
  val active regions, model minus observed:

  | | observed | original SMFNet | k=5 | k=9 |
  | --- | --- | --- | --- | --- |
  | mean pair correlation, sites 0–10 bp apart | 0.832 | −0.042 | −0.042 | −0.043 |
  | share of accessible runs 0–10 bp long | 0.347 | +0.022 | +0.023 | +0.024 |

  - Same model, same 1,000 regions, 32 → 64 steps: the 0–10 bp correlation
    error goes −0.041 → −0.030, and the short-run excess +0.024 → +0.017.
    The share of short inaccessible runs got worse (−0.008 → −0.016).
  - The loss weights low mask rates very heavily (π/2·tan(π/2·(1−t)), up to
    100). That is where architecture gains such as k=9 appear.
- **The seven-model pilot** (1,000 val active regions; `span_infill.py`;
  nats per hidden call; reference = predicting each call from the region's
  other reads, ~0.62):

  | model (SMF body) | 150 bp joint | 150 bp one-shot | 50 bp joint | 20 bp joint |
  | --- | --- | --- | --- | --- |
  | original SMFNet (dilated) | 0.2370 | 0.3946 | 0.2272 | 0.2077 |
  | frozen AlphaGenome (dilated) | **0.2357** | **0.3907** | **0.2259** | 0.2061 |
  | ATAC model (dilated) | **0.2356** | 0.3910 | 0.2260 | **0.2057** |
  | plain U-Net (per-position norm, no attention, AdamW) | 0.2403 | 0.4085 | 0.2294 | 0.2091 |
  | attention U-Net (per-read norm, AdamW) | 0.2428 | 0.4272 | 0.2291 | 0.2084 |
  | attention U-Net, Muon+EMA, k=5 | 0.2452 | 0.4393 | 0.2299 | 0.2082 |
  | attention U-Net, Muon+EMA, k=9 | 0.2432 | 0.4292 | 0.2287 | 0.2078 |
  | patched U-Net, Muon+EMA | 0.2443 | 0.4332 | 0.2295 | 0.2088 |

  These are all on the L40S, so the spans are identical. Paired SEs against
  k=5 are ≤ 0.0004 for joint and ≤ 0.0015 for one-shot at 150 bp.
- **Drift over training** (150 bp, the same spans within each column pair):

  | checkpoint | plain U-Net one-shot / joint | attention U-Net (AdamW) one-shot / joint | k=9 Muon+EMA one-shot / joint |
  | --- | --- | --- | --- |
  | epoch 24 | 0.4029 / 0.2463 | **0.3999** / 0.2441 | 0.3992 / 0.2405 |
  | epoch 99 | 0.4026 / 0.2419 | 0.4207 / 0.2441 | 0.4141 / 0.2414 |
  | epoch 199 | 0.4022 / 0.2407 | 0.4231 / 0.2432 | 0.4158 / 0.2414 |
  | epoch 299/399 | 0.4060 / **0.2401** (399) | 0.4286 / 0.2439 (399) | 0.4211 / 0.2419 (299) |

  - The plain U-Net holds its one-shot score and keeps improving its joint
    score.
  - The attention U-Net with the per-read norm starts ahead and then
    drifts (+0.029). Its joint score never improves, and it ends worse
    than the plain U-Net, despite better val loss (0.2937 vs 0.2946).
  - Optimizer, encoder and k are the same in the first two columns; they
    differ only in the per-read norm and bottleneck attention.
- **50-epoch sweep** (each on the k=9 U-Net, Muon+EMA, 256 × 32, the first
  50 epochs of the 400-epoch schedule; epoch-49 checkpoints, all scored on
  H100s; Δ = paired against the control, with SE):

  | variant | fixed-mask val loss | 50 bp joint Δ | 150 bp joint Δ | 150 bp one-shot |
  | --- | --- | --- | --- | --- |
  | control (k=9 rerun) | 0.29629 | — | — | 0.4123 |
  | original k=9 run, epoch 49 (same config) | 0.29636 | +0.0003 (0.0002) | −0.0004 (0.0001) | 0.4104 |
  | **dilated** (block-2 dilations 4, 4, 2, 1, 1, 1, 1) | **0.29624** | **−0.0012 (0.0003)** | **−0.0020 (0.0002)** | **0.4002** |
  | linear upsampling | 0.29630 | +0.0004 (0.0002) | +0.0001 (0.0001) | 0.4135 |
  | local attention, 4 bp level | 0.29629 | −0.0004 (0.0003) | +0.0001 (0.0002) | 0.4168 |
  | bottleneck at 16 bp (5 levels) | 0.29661 | +0.0011 (0.0003) | +0.0004 (0.0002) | 0.4109 |
  | mask-aware read norm | 0.29638 | +0.0009 (0.0002) | +0.0012 (0.0002) | 0.4179 |
  | windowed norm (±32 positions) | not run (waited >20 h for an H200 on owners) | | | |

  The in-training val curves (fixed masks) agree with these numbers: the
  control and the original run are within 0.0001 at every epoch, and the
  dilated run is slightly ahead throughout.
- **Full runs** (256 × 32, Muon+EMA, 400 epochs; fixed-mask val loss on
  H100s):

  | model | val loss (final = best) | 32-step profile Pearson / corr-of-corr (median) |
  | --- | --- | --- |
  | k=5 attention U-Net | 0.29314 | 0.7598 / 0.6444 |
  | k=9 | 0.29259 | 0.7611 / 0.6444 |
  | **k=9 + dilated** | **0.29229** | 0.7606 / 0.6449 |

  - **Paired, dilated minus k=9, 32 steps** (7,559 regions):
    - profile Pearson +0.0015 (0.0007);
    - residual corr-of-corr +0.0018 (0.0005);
    - corr-of-corr +0.0005 (0.0003);
    - run lengths within noise.
  - **Paired, 64 steps** (1,000 regions): all within noise.
  - **Val loss by mask rate,** dilated minus k=9: −0.0004 to −0.0005 at
    0.05–0.45, −0.0001 at 0.95 and 1.0. The gain is in fine-scale read
    context.
  - **Memorization** (train − val gap): identical to k=9's.
  - **Final inpainting** (identical CPU-drawn spans; Δ vs k=9, with SE):

    | model | 50 bp joint | 150 bp joint | 150 bp one-shot |
    | --- | --- | --- | --- |
    | k=9 + dilated | −0.0006 (0.0003) | −0.0003 (0.0002) | +0.0028 (0.0004) |
    | dilated body (frozen AlphaGenome) | −0.0020 (0.0005) | −0.0055 (0.0003) | −0.0286 (0.0013) |

    The dilated run's own checkpoints drift like k=9's: 150 bp one-shot /
    joint 0.4148 / 0.2413 (epoch 99), 0.4195 / 0.2421 (199), 0.4278 /
    0.2432 (299), 0.4352 / 0.2444 (399).

# Details

## The inpainting scorer (`span_infill.py`)

- **Spans.** For each read and span length S, the start is drawn uniformly
  so the span lies inside the read's covered extent; reads too short for
  the span are skipped. The hidden calls are the read's valid calls inside
  the span; all other calls stay visible.
- **One-shot:** one forward pass.
- **Joint:** the span's calls are revealed one at a time in a random order
  with the true values (teacher forcing), and their log-likelihoods summed.
  That takes as many passes as the longest span has calls. Reads are scored
  in parallel, ~4 minutes per model for 1,000 regions on an L40S.
- **Data.** 1,000 val active regions (`max_regions_per_h5_per_split`),
  32 reads per region, chosen deterministically. Per-region sums are saved,
  so paired SEs can be computed.
- **GPU-type caveat.**
  - `span_infill.py` and `rescore_val.py` draw their spans and masks with a
    CUDA generator, which gives different draws on different GPU types.
    Compare only runs scored on the same type. An early L40S rescore of the
    reference checkpoint looked 0.00055 worse for this reason; the H100
    rescore matches.
  - `span_infill_cpu_draws.py` draws on the CPU, so its spans are identical
    on any GPU. Use it from now on; `rescore_val.py` needs the same change.

## Architecture options (`CNNConfig`, branch nd_unet_inpaint; all default off, channels-last U-Net)

- `unet_norm=read_masked`: per-read statistics weighted by each position's
  local density of observed calls (±8 positions at its level), relative to
  the read's mean density, capped at 1 and floored at 0.01. Reads with no
  observed calls get uniform weights. Plain PyTorch, not fused.
- `unet_norm=window` with `unet_norm_window`: statistics over ±W positions,
  from prefix sums. Plain PyTorch.
- `unet_upsample=linear`: matches
  `F.interpolate(mode="linear", align_corners=False)`, written as
  elementwise ops.
- `unet_down_dilations`: a dilation for each level's second down conv. cuDNN
  runs dilated k=5 convs as fast as undilated ones (1–16, L40S), so this is
  free.
- `unet_local_attention_levels` / `_heads` / `_window`: pre-norm
  self-attention with rotary positions and a zero-initialized output, after
  the up block at the listed levels. Full attention by default; with a
  window, exact halo-block local attention by SDPA.
- **Tests:** `tests/test_inpaint_options.py`.
  - The norm special cases (whole read = GroupNorm(1, C), window 0 =
    LayerNorm, and zero-weight positions don't shift the statistics).
  - Observed weights.
  - Linear upsampling against `F.interpolate`.
  - Halo attention against SDPA with a band mask.
  - fp32 equivalence of channels-last and channels-first with the new
    options.
  - With the existing tests: 304 passed, 1 skipped.

## Commands and jobs

```bash
module load devel pixi/0.53.0; export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
export PYTHONPATH=$HOME/repos/SMFNet_inpaint_runs/src
# inpainting score of one checkpoint (CPU-drawn spans)
$PIXI python inpainting/scripts/span_infill_cpu_draws.py --checkpoint <ckpt> --out <json>
# pilot (7 models) and the two extra U-Nets
bash inpainting/scripts/run_pilot.sh; bash inpainting/scripts/run_pilot_unets.sh
# 50-epoch variants on owners (requeued on preemption; resumes from last.ckpt)
sbatch --job-name=inpaint_<variant> inpainting/scripts/train_variant.sbatch <variant>
# full dilated run + evaluation on btrippe (job 51737387)
sbatch inpainting/scripts/train_eval_k9_dilated.sbatch
```

- **Variant jobs** (owners, H100):
  - 51719888 control, 51719890 linear, 51719893 dilated, 51719896 local
    attention (preempted once, resumed as intended), 51719898 bottleneck16
    and 51719905 read_masked.
  - The windowed norm (51719907) needed an H200 for memory (~72 GiB) and
    never started.
  - The reference epoch-49 scoring job on an H100 was 51734594.
- **Full run:** `ONT_mESC_alphagenome_unet_attn_muon_ema_k9_dilated` in
  `/scratch/users/diamant/smf_models/`.
  - Training took 13.3 h (2026-11-09 20:55 → 2026-11-10 10:10); the
    evaluation finished at 11:16.
- **Analysis:** `fine_scale_compare.py` and `steps_short_range.py` (sampled
  reads vs observed, by distance and run length), and
  `compare_k9_dilated.py` (sampled metrics).
- **Problems on the way:**
  - A 2026-11-09 pixi manifest change in the main checkout (new plotting
    dependencies) blocked `pixi run` until the shared environment was
    updated. The lock diff touched only plotting packages; core imports and
    the GPU tests were checked before submitting.
  - Running the environment's Python directly skips pixi's activation, and
    then Triton cannot find `ptxas` for the fused kernels. Use `pixi run`.

# Related Notes

- [Kernel width vs depth](../2026-11-09/kernel-width.md): the k=9 model all
  variants build on, and the finding that fine-scale context, not receptive
  field alone, drives the low-mask gains. The dilation result here extends
  it.
- [Memorization](../2026-11-08/memorization.md): the train/val diagnostics
  reused here. The dilated run memorizes no more than k=9.
- [Patched U-Net](../2026-11-07/patched-unet.md): its inpainting is scored in
  the pilot; it has the same low-mask deficit.
- [Frozen AlphaGenome encoder](../2026-10-29/frozen-alphagenome-encoder.md):
  the dilated-body models that inpaint best.

# Open Questions

- **Can the U-Net keep its early inpainting** without training on span masks?
  Candidates: the windowed norm (not run), or removing the read-level paths
  only at fine levels.
- **Is the dilated U-Net's lead in val loss and sampled metrics robust?** A
  second seed would tell. The training-run noise in val loss at epoch 49 is
  ~0.0001; at epoch 399 it is unmeasured.
- **Better samplers** that are distribution-preserving: a spread-out
  reveal order, a reallocated step schedule, random-block Gibbs refinement.
  They should shrink the errors all models share, and would show whether
  architecture gains reach the samples.
