---
title: Kernel width vs depth in the attention U-Net, roofline and measured costs, and a k=9 one-block-per-down-level run at matched receptive field
date: 2026-11-09
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - SMFNet_public branch nd_unet_patch (worktree ~/repos/SMFNet_patch; not pushed), commits c05eafc, 4258aa2
  - kernel-width/scripts/ (benchmarks, training and evaluation scripts, outputs)
tags: [unet, kernel-size, receptive-field, roofline, efficiency, muon, ema]
---

# Summary

Wider conv kernels and more conv layers both widen the U-Net's receptive
field.

**Roofline.** A wider kernel adds compute but no memory traffic, and the
U-Net's convs are mostly memory-bound. So widening should be nearly free,
while an extra layer always costs another full memory pass, its
norm/FiLM/GELU kernels and stored activations.

**Measured.** cuDNN does not make wider kernels free: conv time grows close
to linearly in k at the 1 bp level. Per unit of receptive field, though,
widening still beats depth:

| 1 bp attention U-Net change | step time | receptive field | memory |
| --- | --- | --- | --- |
| k=9 (from 5) | +10% | 2× | +0 |
| one extra block per down level | +11% | 1.3× | +4 GiB |

**Matched-receptive-field run.** One k=9 block per down level, with k=5 on
the up path, keeps exactly the default receptive field (1,772 bp).

- **Cost:** 8–9% faster in the benchmark and 17% faster per epoch in real
  training (Muon has fewer matrices to orthogonalize). It uses 4 GiB less
  memory.
- **Quality:** at 128 regions × 64 reads it reached fixed-mask val loss
  0.29356, against 0.29318 for the default.
- **Where it loses:** at low mask rates, 0.1935 vs 0.1927 at mask 0.05.
  Inference: the second nonlinear stage per level matters for read context,
  rather than receptive field.
- **Sampled metrics:** within noise.

**Not finished yet.** A full k=9 run with two down blocks (2× receptive
field) is training as btrippe job 51631889.

# Key Points

- **Roofline arithmetic.** For a "same" conv of N reads × L positions × C
  channels, kernel k, bf16:
  - each of forward, input-gradient and weight-gradient does 2·N·L·C²·k FLOPs;
  - memory traffic is ~4·N·L·C bytes per pass, independent of k;
  - arithmetic intensity ≈ C·k/2 FLOP/byte, against an H100 ridge point of
    ~990 TFLOP/s ÷ 3.35 TB/s ≈ 295.

  So convs are memory-bound for C·k < ~590. At the 1 bp model's level 0
  (C=48) that means k ≲ 12; at the patched model's level 0 (C=128), k ≲ 4.
- **Per-layer measurements** (`bench_kernel_layers.py`; 8,192 reads,
  forward + backward, channels-last conv2d as the model runs them):

  | level | k=3 | k=5 | k=9 | k=15 | k=25 | roofline at k=5 / 9 / 15 |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 bp U-Net, L=1024, C=48 | 2.4 | 3.3 | 5.8 | 10.5 | 16.7 ms | 1.4 / 1.4 / 1.8 |
  | 1 bp U-Net, L=512, C=64 | 1.2 | 1.4 | 2.3 | 4.9 | 7.9 | 1.0 / 1.0 / 1.6 |
  | patched, L=256, C=128 | 1.2 | 1.6 | 2.7 | 4.6 | 8.1 | 1.0 / 1.9 / 3.1 |

  - **cuDNN runs these convs at 15–20% of peak**, because so few channels
    leave its implicit GEMMs too narrow. So time grows roughly linearly with
    k even where the roofline says it shouldn't.
  - **Same receptive field, one wide block vs two narrow ones** (norm + conv
    blocks): at the 1 bp level 0, one k=9 block takes 8.4 ms against
    2 × 5.9 ms for two k=5 blocks. At the patched model's level 0 it is
    4.5 ms against 2 × 3.8 ms.
- **Full model** (`bench_variants.py`; 1 bp attention U-Net, frozen
  AlphaGenome encoder, AdamW, no EMA). The conv-path receptive field (RF),
  ignoring attention, is (k_down−1)·(down blocks + extra)·127 +
  (k_up−1)·189 bp: 127 is the sum of the down levels' resolutions
  (1 + 2 + … + 64 bp) and 189 comes from the up path.

  | variant | RF | RF at 1 bp | 256 × 32 step / peak | 128 × 64 step / peak |
  | --- | --- | --- | --- | --- |
  | k=5, 2 down blocks (default) | 1,772 bp | 12 bp | 184 ms / 36.3 GiB | 164 / 32.8 |
  | k=3 | 886 | 6 | 176 / 36.3 | |
  | k=7 | 2,658 | 18 | 193 / 36.3 | |
  | k=9 | 3,544 | 24 | 202 / 36.3 | 183 / 32.8 |
  | k=15 | 6,202 | 42 | 238 / 36.4 | |
  | k=5, +1 block per down level | 2,280 | 16 | 204 / 40.4 | |
  | k=5, +2 blocks per down level | 2,788 | 20 | 225 / 44.3 | |
  | **k=9, 1 down block, up k=5** | **1,772** | 12 | **168 / 32.4** | **149 / 28.9** |
  | k=9, 1 down block, up k=9 | 2,528 | 16 | 176 / 32.4 | 157 / 28.9 |
  | k=13, 1 down block, up k=5 | 2,280 | 16 | 174 / 32.4 | 155 / 28.9 |

  In the patched U-Net, k=9 and one extra down block each cost ~+9%
  (165 vs 151 ms), but the extra block adds 2.4 GiB.
- **The bottleneck attention already sees the whole read** (at 64 bp), so
  the conv RF mainly sets the fine-scale context each level integrates
  before pooling.

# Details

## Training run: k=9, one down block, up k=5, at 128 × 64

- **Run.** `ONT_mESC_alphagenome_unet_attn_muon_ema_r64_k9d1u5` with
  `model.cfg.cnn_config.unet_kernel_size=9 unet_down_blocks=1
  unet_up_kernel_size=5`.
  - Otherwise like `ONT_mESC_alphagenome_unet_attn_muon_ema_r64`: Muon
    lr 0.01, EMA 0.999, 128 regions × 64 reads, `trainer.max_epochs=200
    model.cfg.warmup_epochs=5` (182,000 steps).
  - Ran in Sherlock job 51468982 with
    [train_eval_attn_muon_ema_r64_k9d1u5.sh](kernel-width/scripts/train_eval_attn_muon_ema_r64_k9d1u5.sh),
    then the standard evaluation.
- **Time.** Training took 10.1 h, against 11.5 h for the k=5 run at the
  same batch shape (180 vs 218 s per epoch).

| | k=5, 2 down blocks (128 × 64) | k=9, 1 down block, up k=5 (128 × 64) |
| --- | --- | --- |
| fixed-mask val loss (final = best, epoch 199) | **0.29318** | 0.29356 |
| val cross-entropy at mask 0.05 / 0.45 / 0.95 | 0.1927 / 0.2514 / 0.5176 | 0.1935 / 0.2519 / 0.5177 |
| train − val gap, mask 0.05 / 1.0 | +0.0048 / +0.0013 | +0.0048 / +0.0013 |

- **Paired, k=9 one-block minus k=5, 32 steps** (7,559 val active regions):
  - profile Pearson −0.0007 (0.0006);
  - corr-of-corr −0.0001 (0.0003);
  - run-length zeros +0.009 (0.017), ones +0.021 (0.016).
- **Paired, 64 steps** (1,000 regions): Pearson −0.0027 (0.0017),
  corr-of-corr +0.0002 (0.0007).
- **Reading.**
  - The loss is concentrated where the read's own nearby calls matter,
    and is zero when nearly everything is masked. That is the same signature
    as the patched U-Net's deficit. Fine-scale read context seems to want
    nonlinear depth more than receptive field.
  - Memorization is unchanged.
  - The step-time saving (12% wall clock) costs about as much val loss as
    the AdamW-vs-Muon difference.

## Commands

```bash
module load devel pixi/0.53.0; export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
cd ~/repos/SMFNet_patch; export PYTHONPATH=$PWD/src
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
$PIXI python kernel-width/scripts/bench_kernel_layers.py
$PIXI python kernel-width/scripts/bench_variants.py base_k9_d1_u5 1024 128 64   # variant [length regions reads]
bash $SCRATCH/smfnet_kernel/train_eval_attn_muon_ema_r64_k9d1u5.sh
sbatch $SCRATCH/smfnet_kernel/train_eval_attn_muon_ema_k9.sbatch            # job 51631889
```

The new options (defaults leave models unchanged) are `unet_kernel_size`,
`unet_extra_blocks`, `unet_down_blocks` (1 or 2) and `unet_up_kernel_size`.
They come with fp32 equivalence tests in `tests/test_fast_unet.py`.

# Related Notes

- [Patched U-Net](../2026-11-07/patched-unet.md): the same low-mask-rate
  deficit when fine-scale capacity is cut, and the batch-shape results
  (128 × 64) used here.
- [Training-step speed optimizations](../2026-11-06/unet-training-speed.md):
  the roofline framing and the channels-last conv path benchmarked here.
- [Memorization](../2026-11-08/memorization.md): the diagnostics and the
  128 × 64 k=5 run this run is compared with.

# Open Questions

- **Does more receptive field help?** That is the k=9 run with two down
  blocks (2× RF; btrippe job 51631889, 256 × 32, 400 epochs, compared with
  0.29314). Its evaluation runs automatically after training.
- **Could a custom conv kernel make wide kernels nearly free at C=48–64?**
  For example, a Triton direct conv that loads a tile once for all k taps.
  The roofline says k≈12 should cost the same as k=5 at the 1 bp level.
- **Speed without losing fine-scale depth:** cheaper second blocks, for
  example depthwise convs, instead of removing them.
