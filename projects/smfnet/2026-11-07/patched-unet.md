---
title: Patched U-Net (4 bp patches, 1 bp head) matches the attention U-Net early in training at 0.6-0.8x the step and 0.7x the memory; Muon+EMA vs AdamW+EMA final results
date: 2026-11-07
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - SMFNet_public branch nd_unet_patch (worktree ~/repos/SMFNet_patch; not pushed), commits c8a8db5, 7889c57, on nd_unet_fast
  - src/smf_net/cnn.py (CNN._forward_patched, UNetConvBlock.forward_channels_last), scripts/train.py (ForeachEMAWeightAveraging, StopAfterEpochs)
  - patched-unet/scripts/
tags: [unet, efficiency, roofline, long-range, ema, muon, adamw, evaluation]
---

# Summary

**Patched U-Net.** Most of the attention U-Net's per-read memory traffic is
at 1 and 2 bp, where the activations are largest and the channel counts too
low for good tensor-core use. The new `unet_patch_size` option runs the
U-Net on 4 bp patches, so its first level is at 4 bp. A light head then
returns to 1 bp, adding:

- the U-Net's patch features, unpatched;
- a projection of the 1 bp AlphaGenome features, computed once per region;
- a projection of the input calls;
- optional residual k=9 conv blocks.

Results, against the 1 bp attention U-Net:

- **The best variant** (`p4_wide_h2`: channels 128→256 at 4–64 bp, plus two
  32-channel head blocks) matched it on val loss at every epoch of the first
  30, and at every mask level. Its step is 151 vs 183 ms and its peak memory
  25 vs 36 GiB.
- **Without the head blocks** (`p4_wide`, 135 ms), it trails by ~0.0004.
  The gap sits at low mask rates, where 1 bp detail of nearby observed calls
  matters.
- **Cost scales with positions only.** An 8 kb window at the same number of
  positions per batch costs the same as 1 kb. This is the property needed for
  long-range modeling.

**Optimizer.** With both 400-epoch attention U-Net + EMA runs now evaluated,
AdamW+EMA is **not** better than Muon+EMA:

- **Fixed-mask val loss of the final weights:** Muon+EMA 0.29314,
  AdamW+EMA 0.29361.
- **32-step sampling:** AdamW+EMA is slightly worse on profile Pearson
  (−0.0015 ± 0.0007).
- **64-step sampling:** the two tie on every metric.

The impression that AdamW+EMA was better came from its mid-training lead
(for example, val loss 0.2953 vs 0.2961 at epoch 100). The Muon+EMA run
also logged noisy in-training val loss, because it started before
validation masks were made deterministic.

**EMA speed.** Lightning's `EMAWeightAveraging` updates one tensor at a time,
costing ~14 ms per step with the frozen AlphaGenome encoder's tensors.
`ForeachEMAWeightAveraging` does the same update in ~3 ms, so real training
of `p4_wide_h2` went from 173 to 163 ms per step.

# Key Points

- **Patched variants**, 30 epochs of the 400-epoch schedule, AdamW (lr as in
  the baseline) + EMA 0.999, and the same deterministic val masks as the
  baseline run:

  | variant | U-Net channels (levels) | step (bench) | peak | train s/epoch | val loss, epoch 9 | 19 | 29 |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | 1 bp attention U-Net (`..._unet_attn_adamw_ema`) | 48→144 (1–64 bp) | 183 ms | 36.3 GiB | 173* | 0.3055 | 0.3001 | **0.2986** |
  | p4_narrow | 80→144 (4–64 bp) | 103 ms | 16.9 GiB | 66 | 0.3071 | 0.3020 | 0.3001 |
  | p4_mid | 96→224 | 118 ms | 21.0 GiB | 73 | 0.3062 | 0.3010 | 0.2994 |
  | p4_wide | 128→256 | 135 ms | 24.6 GiB | 82 | 0.3060 | 0.3006 | 0.2989 |
  | p2 | 64→144 (2–64 bp) | 135 ms | 25.6 GiB | 81 | 0.3061 | 0.3006 | 0.2991 |
  | **p4_wide_h2** (+2 head blocks, 32 ch) | 128→256 | 151 ms | 25.2 GiB | 89 (85†) | 0.3057 | 0.3003 | **0.2986** |

  \* The baseline ran on the old code without the speed options. With them
  its steady-state step is ~189 ms (see
  [the speed note](../2026-11-06/unet-training-speed.md)), so its epoch would
  be roughly 95 s. That figure is estimated, not measured.

  † 85 s with the foreach EMA.

  Bench steps are compiled, bf16, 256 regions × 32 reads, with AdamW and
  without EMA.
- **Where the patch variants lose** (val cross-entropy by mask rate, epoch 29):

  | run | 0–0.2 | 0.2–0.4 | 0.4–0.6 | 0.6–0.8 | 0.8–1.0 |
  | --- | --- | --- | --- | --- | --- |
  | 1 bp baseline | 0.2088 | 0.2318 | 0.2686 | 0.3323 | 0.5122 |
  | p4_narrow | 0.2111 | 0.2339 | 0.2702 | 0.3335 | 0.5126 |
  | p4_wide | 0.2092 | 0.2325 | 0.2691 | 0.3326 | 0.5123 |
  | p4_wide_h2 | 0.2089 | 0.2320 | 0.2687 | 0.3323 | 0.5121 |

  The narrow-model gap is ~0.002 at low masking and ~0 when nearly
  everything is masked. Inference: what is lost is 1 bp detail of nearby
  observed calls. The head blocks restore it.
- **Per-position norm** (channels-last, compiled), as an alternative to the
  fused per-read norm: not faster (190 vs 183 ms at 1 bp, 121 vs 118 ms for
  p4_mid), and more memory (61 vs 36 GiB at 1 bp), because inductor saves
  its intermediates. The fused per-read norm stays.
- **Long range** (8,192 bp SMF window × 32 regions × 32 reads, the same
  positions per batch as 1,024 bp × 256 regions):

  | variant | 1 kb step / peak | 8 kb step / peak |
  | --- | --- | --- |
  | 1 bp attention U-Net | 183 ms / 36.3 GiB | 185 ms / 36.3 GiB |
  | p4_narrow | 103 / 16.9 | 101 / 16.5 |
  | p4_wide | 135 / 24.6 | 134 / 24.0 |

  - **Attention stays cheap.** It runs at the 64 bp bottleneck, which is 128
    tokens per read at 8 kb, so cost is linear in positions for every
    variant. The patched U-Net holds about 1.5–2× more positions per GPU.
  - **Caveat: the fused-norm forward path.** At 8 kb × 32 regions there are
    only 1,024 reads, below the single-pass threshold of 2,048, so the fused
    norm uses its chunked forward (as designed).
- **Muon+EMA vs AdamW+EMA**, attention U-Net, both 400 epochs:

  | | fixed-mask val loss, best ckpt | final (epoch 399) | profile Pearson (32 steps) | corr-of-corr (32) | Δ AdamW−Muon at 64 steps, 1,000 regions |
  | --- | --- | --- | --- | --- | --- |
  | Muon+EMA (best: epoch 377) | 0.29318 | **0.29314** | **0.7598** | 0.6444 | — |
  | AdamW+EMA (best: epoch 396) | 0.29361 | 0.29361 | 0.7556 | 0.6452 | Pearson +0.0006 ± 0.0018; corr-of-corr +0.0003 ± 0.0007; run lengths within 1 SE |
  | AdamW, no EMA (epoch 366) | 0.29365 | | 0.7568 | 0.6456 | |

  - **Paired per-region at 32 steps** (7,559 val active regions),
    AdamW+EMA minus Muon+EMA:
    - profile Pearson −0.0015 (0.0007);
    - corr-of-corr −0.0001 (0.0003);
    - run-length zeros +0.035 (0.019), ones −0.013 (0.017).
  - **Fixed-mask** means the same masks for every model (`rescore_val.py`),
    over all 13,922 val regions.
  - **Medians are noisy.** Paired differences are the reliable comparison.

# Details

## The patched U-Net (`CNNConfig.unet_patch_size`, `unet_head_channels`, `unet_head_blocks`)

Implemented in `CNN._forward_patched`. It needs the channels-last U-Net.

1. **Stem.**
   - Each read's (L, 3) one-hot calls are reshaped to (L/p, 3p).
   - The input conv runs on patches. Its kernel is
     `input_kernel_size // p`, rounded to odd: 21 bp → 5 patches.
2. **U-Net** levels at p, 2p, … bp.
   - With the AlphaGenome encoder, the encoder returns
     `len(unet_channels) + log2(p)` levels from 1 bp. The U-Net uses those
     from p bp, so it still has FiLM at every level from AlphaGenome features
     at the matching resolution, plus the full-window 64 bp tokens for
     cross-attention.
   - AlphaGenome's maps stop at 64 bp, so p=4 allows at most 5 U-Net levels.
   - With the CNN encoder, the 1 bp map is average-pooled by p.
3. **Head at 1 bp**, `unet_head_channels` D wide (16 by default; 32 in h2):
   - `Linear(C0 → p·D)` on the U-Net output, unpatched to (L, D);
   - plus `Linear(512 → D)` of the 1 bp AlphaGenome map, once per region
     and broadcast over reads;
   - plus `Linear(3 → D)` of the input calls;
   - then `unet_head_blocks` residual `UNetConvBlock`s (fused per-read norm
     → GELU → weight-standardized k=9 conv, residual scale 0.15);
   - then GELU and the k=75 output conv as Toeplitz matmuls.
4. **Old checkpoints** load unchanged. All three options default to
   patch 1, D=16 and 0 blocks, which is the original model.
5. **Tests** (`tests/test_fast_unet.py`):
   - patched forward with both norms: shapes, `dna`/`dna_embed`
     equivalence, gradients reach the stem and every head part, and each
     read's output depends only on that read;
   - the channels-last per-position U-Net against the channels-first one in
     fp32.

   The full suite gives 293 passed, 1 skipped.

## Foreach EMA (`scripts/train.py: ForeachEMAWeightAveraging`)

- **Why it was slow.** `EMAWeightAveraging(decay=…)` hard-codes `avg_fn`,
  so `AveragedModel.update_parameters` loops over every parameter and buffer.
  That is 316 tensors and 101M values, 94M of them the frozen AlphaGenome
  encoder's.
- **The subclass** passes `multi_avg_fn=get_ema_multi_avg_fn(decay)`, one
  `torch._foreach_lerp_` per dtype group: 13.9 → 2.8 ms per update
  (`bench_ema.py`).
- **Checkpoints.** It reports `state_key = "EMAWeightAveraging"`, so
  existing EMA runs resume with it. On a toy model it gives the same average
  as Lightning's callback (`test_ema_callback.py`).
- **Real training of p4_wide_h2:** 173 → 163 ms per step (tqdm timestamps,
  epochs 1–2; the bench is 151 ms). About 12 ms per step of Lightning or
  data overhead is still unexplained for the patched models.

## Commands

```bash
module load devel pixi/0.53.0; export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
cd ~/repos/SMFNet_patch; export PYTHONPATH=$PWD/src
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
# step time and memory: <variant> [smf_length regions]
$PIXI python patched-unet/scripts/bench_variants.py p4_wide_h2
$PIXI python patched-unet/scripts/bench_variants.py p4_wide 8192 32
# 30-epoch runs (one after another on this GPU), then curves
bash patched-unet/scripts/run_variants.sh p4_narrow p4_wide p2 p4_mid p4_wide_h2
python3 patched-unet/scripts/compare_curves.py; python3 patched-unet/scripts/compare_mask_bins.py 29
```

- **Run directories:**
  `/scratch/users/diamant/smf_models/ONT_mESC_unet_attn_adamw_ema_<variant>`.
  Each run is `scripts/train.py` with the common attention U-Net arguments
  (`$SCRATCH/smfnet_unet/common_args.txt`), the three speed options,
  `ema_decay=0.999 stop_after_epochs=30`, and the variant's options (see
  `run_variants.sh`). wandb runs are named `ONT_mESC_unet_attn_adamw_ema_<variant>`.
- **AdamW+EMA evaluation:**
  [eval_adamw_ema.sh](patched-unet/scripts/eval_adamw_ema.sh) runs the
  fixed-mask rescoring, `eval_run.sh` (32 steps, all val active regions,
  metrics, mask-rate diagnostics) and the 64-step 1,000-region subset. Then
  [after_eval.sh](patched-unet/scripts/after_eval.sh) rescores the epoch-399
  weights, tests and benchmarks h2, and trains it.
- **Comparison outputs:** `$SCRATCH/smfnet_ag_smoke/compare_optim_ema.txt`
  and `$SCRATCH/smfnet_unet/compare_steps_with_adamw_ema.txt`.
- **Jobs.** The first four patched runs ran in Sherlock job 51222511 and
  all finished before it ended. The evaluation, h2 and the EMA timing ran in
  job 51468982.

## Batch shape: regions × reads per region

- **Benchmark** (`bench_variants.py <variant> 1024 <regions> <reads>`):
  - More regions per batch doesn't raise throughput: 512 × 32 is only 3%
    faster per read than 256 × 32, because the GPU is already saturated.
  - More reads per region shares the per-region work across more reads
    (AlphaGenome encoder, FiLM projections, cross-attention keys and values).
    For p4_wide_h2, reads per second rise 54k → 61k at 128 × 64, 65k at
    64 × 128 and 67k at 128 × 128. The 1 bp model gains about the same
    (+12% to +22%).
- **Data.** Training regions have 138 reads (median; 5–95%: 111–171), so
  an epoch at 32 reads sees ~23% of the reads. At 128 reads per region, 25%
  of regions would have to repeat reads.
- **Matched test.** p4_wide_h2 at 128 regions × 64 reads, with
  `trainer.max_epochs=200 model.cfg.warmup_epochs=5`, stopped after 15
  epochs. That gives the same steps (13,650), reads per step, warmup length
  in steps and LR schedule as the 256 × 32 run's first 30 epochs.
  `warmup_epochs` is a new `SMFNetConfig` option (commit 05f0aea); its
  default keeps the old rule.

  | batch | real step | train s/epoch | fixed-mask val loss (32 reads/region) at step 13,650 |
  | --- | --- | --- | --- |
  | 256 × 32 | 163 ms* | 85* (455 steps) | 0.29861 |
  | 128 × 64 | 149 ms | 150 (910 steps) | 0.29868 |

  \* With the foreach EMA, measured in the 3-epoch timing run.

  The two are equal at this point in training, and 128 × 64 is ~9% faster
  in real training. The in-training val loss of the 128 × 64 run uses 64
  reads per region, so it is not comparable to the other runs; hence the
  rescoring (`run_r64.sh`).

# Related Notes

- [Training-step speed optimizations](../2026-11-06/unet-training-speed.md):
  the speed options all variants here build on, and the roofline argument
  that the 1 bp levels are bandwidth-bound.
- [U-Net SMF body](../2026-11-03/unet-smf-body.md): the U-Net design
  (AlphaGenome blocks, per-level FiLM, residual scale 0.15) that the patched
  variant keeps from 4 bp up.
- [Cheri blocks and Muon in SMFNet](../2026-11-03/cheri-muon-smfnet.md):
  where Muon was introduced. Its learning-rate findings are behind the
  Muon+EMA run compared here.

# Open Questions

- **Does p4_wide_h2 still match the 1 bp model at 400 epochs?** That needs
  a full run (~11 h at ~100 s per epoch including validation) and sampling.
  Muon+EMA is the better optimizer on the final evaluation, so a full run
  should probably use it, or both.
- **Long-range data and levels.** Training at 8–10 kb needs data regions of
  that length. Deeper U-Nets (128–256 bp levels) would need pooled
  AlphaGenome features past 64 bp.
- **Remaining overhead.** What is the ~12 ms per step between bench and real
  training for the patched models (real data, Lightning, loader)?
- **Head size.** Is two 32-channel head blocks the right size? One block or
  16 channels might be enough, and cheaper.

# Sources

- SMFNet_public `nd_unet_patch` commits c8a8db5 (patched U-Net,
  channels-last per-position norm, `stop_after_epochs`) and 7889c57 (head
  blocks, foreach EMA). Worktree `~/repos/SMFNet_patch`, not pushed.
- Scripts: [patched-unet/scripts/](patched-unet/scripts/).
- Lightning `pytorch_lightning/callbacks/weight_averaging.py` and
  `torch/optim/swa_utils.py` (PyTorch 2.5) in the pixi environment.
