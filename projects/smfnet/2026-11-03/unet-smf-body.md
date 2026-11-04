---
title: U-Net SMF body with AlphaGenome-style blocks and per-level FiLM from frozen AlphaGenome features
date: 2026-11-03
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - SMFNet_public branch nd_unet, commit d8cc2dc (worktree ~/repos/SMFNet_unet; not pushed)
  - src/smf_net/cnn.py (UNetConvBlock, UNetDownBlock, UNetUpBlock, UNetBody, AlphaGenomeDNAEncoder multiscale_levels)
  - alphagenome-pytorch commit 72268c0, src/alphagenome_pytorch/convolutions.py (DownResBlock, UpResBlock)
  - unet-smf-body/scripts/bench_unet.py, unet-smf-body/scripts/train_unet.sh
tags: [unet, smf-body, efficiency, alphagenome, film, initialization]
---

# Summary

After the frozen AlphaGenome encoder proved both fast and effective (see
[the AlphaGenome note](../2026-10-29/frozen-alphagenome-encoder.md)), the goal
shifted toward efficiency. The per-read SMF body dominates SMFNet's compute and
memory: it runs on 256 regions × 32 reads = 8,192 reads per batch.

This note records the design of a U-Net SMF body built from AlphaGenome's
down and up residual blocks, conditioned at every level by the frozen
AlphaGenome features at the same resolution. At about the same parameter
count as the current dilated body (2.4M vs 2.2M), a training step is about
1.8× faster (329 vs 578 ms) and uses less memory (48.9 vs 58.4 GiB).

The first smoke test was unstable. Starting each residual branch's learned
scale at 0.15 instead of AlphaGenome's 1 fixed it, and the fixed U-Net beat
the dilated body on the same smoke test. A full 400-epoch run started
2026-11-03; results are pending.

# Key Points

- **Why a U-Net should be cheaper.** A conv layer costs about parameters ×
  positions. Putting most layers at reduced resolution cuts compute at equal
  parameter count. Estimated per read, conv FLOPs and stored layer outputs
  only:

  | SMF body | params | forward FLOPs | stored layer outputs |
  | --- | --- | --- | --- |
  | current (8 dilated k=9 blocks, 128 channels at 1 bp) | 1.2M | 1× (2.45 GFLOP) | 1× |
  | U-Net, channels 48 → 144 at 1 → 64 bp | 1.3M | 0.14× | 0.53× |
  | U-Net, channels 64 at 1 bp then 128 → 384 | 7.8M | 0.54× | 0.94× |

  Memory falls less than compute, because the full-resolution levels hold the
  largest activations. The width at 1 bp is the main memory knob.
- **Measured** (full batch, bf16, compiled, H100, frozen AlphaGenome encoder):

  | SMF body | body params | trainable | step | peak memory |
  | --- | --- | --- | --- | --- |
  | dilated (current) | 2.24M | 8.8M | 578 ms | 58.4 GiB |
  | U-Net 48 → 144 | 2.43M | 6.6M | **329 ms** | **48.9 GiB** |
  | U-Net 64, 128 → 384 | — | — | out of memory | — |

  The wide variant runs out of memory because each U-Net level stores
  several activations, and doubling channels as length halves keeps
  per-level memory roughly constant.
- **Why the design suits SMF beyond efficiency (inferred, untested):**
  - Reads look like a segmentation problem: accessible stretches, footprints
    and nucleosomes. Skip connections keep 1 bp precision, while coarse
    levels see nucleosome-scale (~150 bp) structure and spacing.
  - Coarse levels that see the whole read are a plausible route to better
    within-molecule covariation. No conditioning change so far has closed
    more than ~13% of the covariation gap.

# Details

## Design choices

1. **Blocks: AlphaGenome's `ConvBlock`** (norm → GELU → conv). This is
   `UNetConvBlock`:
   - **Weight-standardized k=5 convs** (`StandardizedConv1d`, from the
     copied AlphaGenome module); a plain 1×1 conv for the pointwise blocks.
   - **AlphaGenome's GELU** approximation, `sigmoid(1.702x)·x`.
   - **SMFNet's channel LayerNorm (`CNNLayerNorm1D`) in place of
     AlphaGenome's `RMSBatchNorm`.** `RMSBatchNorm` divides by stored
     statistics that are never updated, so trained from scratch it is only
     an affine map and normalizes nothing.
2. **Down path: AlphaGenome's `DownResBlock`** (`UNetDownBlock`). Two k=5
   conv blocks per level. The residual is the input zero-padded to the new
   width, so channels must not decrease with depth.
   - **Downsampling is 2× average pooling,** not AlphaGenome's max pooling.
     The input is sparse one-hot calls with mask tokens, and max pooling
     after the k=21 input conv could let lone sites dominate. Strided convs
     were judged equally reasonable; only average pooling is implemented.
3. **Up path: AlphaGenome's `UpResBlock`** (`UNetUpBlock`): a conv block
   plus a channel-sliced residual, 2× upsampling by `repeat_interleave`, a
   learned residual scale (init 1), the skip connection added through a 1×1
   conv block, then a final k=5 conv block. There is no bottleneck module;
   AlphaGenome's transformer sits there.
4. **Levels:** 7 levels at 1, 2, 4, …, 64 bp, with channels
   `[48, 64, 80, 96, 112, 128, 144]` (`CNNConfig.unet_channels`). Level 0's
   width replaces `n_filters` for the input conv (k=21) and the output conv
   (k=75). The read length (1,024) must be divisible by 64.
5. **Conditioning: FiLM at every level, at the matching resolution.**
   - **Where:** each down level's first conv block and each up level's last
     conv block modulate their normalized input (adaLN-style placement):
     `h·(1 + s) + b`.
   - **The projection:** `CNNLayerNorm1D` → 1×1 conv to 2C, initialized to
     zero so the U-Net starts unconditioned. It runs once per region and is
     broadcast over the region's reads with the existing `FiLM` class.
   - **With the AlphaGenome encoder,** `AlphaGenomeDNAEncoder(multiscale_levels=7)`
     returns each level's existing projection (`CNNLayerNorm1D` → 1×1 conv
     to 512), each cropped to the middle half at its own resolution,
     instead of their sum at 1 bp. Level i of the U-Net is conditioned on
     AlphaGenome's 2ⁱ bp features. The `mix` block is not built in this mode.
   - **With the CNN DNA encoder,** the U-Net average-pools the 1 bp
     embedding to each level.
6. **Initialization fix: residual branches start small.** The first smoke
   test spiked as the learning rate warmed up, with step losses of 1.87 and
   1.97 and val/loss 2.36 after one epoch. Diagnosis: about 30 residual
   blocks of weight-standardized convs at output scale 1 make the residual
   stream grow quickly. SMFNet's default-initialized convs add much less,
   and AlphaGenome only ever runs its blocks with trained scales.
   - **Fix:** `StandardizedConv1d.scale` starts at
     `UNetConvBlock.RESIDUAL_SCALE_INIT = 0.15`, Cherimoya's residual scale.
   - **Smoke test** (2,000 regions per h5 per split, 4 epochs, batch 64,
     1-epoch warmup):

     | | epoch 0 | 1 | 2 | 3 |
     | --- | --- | --- | --- | --- |
     | U-Net, scale init 1 | 2.360 | 0.447 | 0.367 | 0.352 |
     | U-Net, scale init 0.15 | 0.503 | 0.352 | 0.338 | 0.332 |
     | dilated body, frozen AlphaGenome | 0.530 | 0.366 | 0.344 | 0.337 |

     One step spike (1.14) remained during the smoke test's short warmup.
     The full run warms up over 10 epochs.
7. **Plumbing.**
   - `CNN.forward` accepts either a single DNA-embedding tensor or the
     per-level list. It casts each to the SMF dtype and checks the finest
     level's length.
   - PRISM gets the batch size of either form via `_embed_batch_size`.
   - `CNNConfig.smf_body_type` (`"dilated"` by default, or `"unet"`) and
     `unet_channels` are new fields; their defaults leave existing models
     unchanged.
   - Tests in `tests/test_unet.py` cover shapes, the dna/dna_embed
     equivalence, per-level AlphaGenome embeddings with frozen weights, and
     zero-initialized FiLM that broadcasts over reads.

## Full run (started 2026-11-03)

Run `/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_frozen_unet`, wandb
`j7k71gli`. It runs inside Sherlock interactive job 51001983 on
`sh04-03n07`, launched detached with
[scripts/train_unet.sh](unet-smf-body/scripts/train_unet.sh). It runs from the
`nd_unet` worktree with the main checkout's pixi environment:

```bash
module load devel pixi/0.53.0; export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
cd ~/repos/SMFNet_unet
PYTHONPATH=$PWD/src pixi run --frozen --manifest-path ~/repos/SMFNet_public/pyproject.toml \
    python scripts/train.py hydra.run.dir=/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_frozen_unet \
    ... model.cfg.cnn_config.dna_encoder_type=alphagenome model.cfg.cnn_config.smf_body_type=unet \
    alphagenome_weights=/scratch/users/diamant/alphagenome_pytorch/model_all_folds.safetensors \
    log_mask_rate_bins=true checkpoint_every_n_epochs=25
```

The data and training settings match the frozen-AlphaGenome run. It trains
at about 2.8 steps/s (≈2.7 min per epoch, ≈19 h for 400 epochs).
Validation loss at epochs 0, 1 and 2: 0.3444, 0.3191, 0.3161. The
frozen-AlphaGenome dilated run had 0.3670, 0.3250, 0.3219, and the
baseline 0.3623, 0.3297, 0.3233. These are early numbers only.

The benchmark is [scripts/bench_unet.py](unet-smf-body/scripts/bench_unet.py),
run with `PYTHONPATH=~/repos/SMFNet_unet/src`.

# Related Notes

- [Frozen AlphaGenome encoder](../2026-10-29/frozen-alphagenome-encoder.md):
  the DNA encoder and the five-model comparison (baseline, frozen and
  fine-tuned AlphaGenome, ATAC) this run will be judged against. It also
  documents the multi-resolution AlphaGenome features reused here.
- [Cheri blocks and Muon in SMFNet](cheri-muon-smfnet.md): the alternative
  efficiency direction tried just before. It ran out of memory without
  checkpointing and was about 2× slower than the baseline. It is also the
  source of the 0.15 residual scale used in the fix above.
- [Positional baselines and noise-matched corr-of-corr](../2026-09-29/positional-covariation-baselines-and-noise-matched-corr-of-corr.md):
  defines the covariation metrics that the U-Net is hoped to improve.

# Open Questions

- Does the U-Net match or beat the frozen AlphaGenome + dilated body run on
  val/loss, sampled profiles and covariation, at about 0.57× the step time?
- Does it improve within-molecule covariation (corr-of-corr), as
  hypothesized?
- If it matches: how far can it be widened at about the current compute?
  The 64 → 384 variant needs less 1 bp width or activation checkpointing to
  fit in memory.
- Average pooling vs strided convs, and the number of levels (stopping at
  32 or 128 bp), are untested.
- Is the fixed 0.15 residual scale init close to optimal? Alternatives are
  zero init (SkipInit/ReZero) or scaling by depth.

# Sources

- SMFNet_public `nd_unet` commit d8cc2dc: `src/smf_net/cnn.py`,
  `src/smf_net/evaluation/prism.py`, `tests/test_unet.py`.
- alphagenome-pytorch commit 72268c0:
  https://github.com/genomicsxai/alphagenome-pytorch (`convolutions.py`,
  `model.py`).
- Cherimoya residual scale: https://github.com/jmschrei/cherimoya (commit
  5729b57, `cherimoya/cheri.py`).
- [scripts/bench_unet.py](unet-smf-body/scripts/bench_unet.py),
  [scripts/train_unet.sh](unet-smf-body/scripts/train_unet.sh).
