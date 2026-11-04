---
title: Cheri blocks and Muon in SMFNet, the fused Triton kernel on PyTorch 2.5, memory, and Muon learning rate
date: 2026-11-03
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - SMFNet_public branch nd_cheri, commits de82785 and 2379432 (worktree ~/repos/SMFNet_cheri; not pushed)
  - https://github.com/jmschrei/cherimoya (commit 5729b57; cherimoya/cheri.py, cherimoya/training.py)
  - torch.optim.Muon in PyTorch 2.9.1 (torch/optim/_muon.py)
  - cheri-muon-smfnet/scripts/ (bench_cheri.py, train_cheri_muon.sh, train_cheri_muon_lr0.01.sbatch, early_curves.py)
tags: [cherimoya, cheri-block, muon, triton, memory, learning-rate, film]
---

# Summary

Both SMFNet bodies (the DNA encoder and the per-read SMF body) were rebuilt
from Cherimoya's Cheri blocks, with FiLM conditioning in the SMF body, and
the blocks' MLP weights were trained with Muon.

- **The kernel works on PyTorch 2.5 / Triton 3.1.** Cherimoya's fused
  conv + norm training kernel, copied in, matches its pure-PyTorch reference
  to about 1e-6 in fp32. The inference megakernel was not ported.
- **Memory is the binding constraint.** The plain Cheri model ran out of
  memory even at half the batch. Per-block activation checkpointing fits
  the full batch (48.9 GiB), at about 2× the baseline's step time.
- **Muon's learning rate.** With Cherimoya's default (0.025), training loss
  stalled as the rate warmed up to its peak, and the baseline caught up by
  epoch ~16. At 0.01 there was no stall: an early lead in val/loss (−0.0037
  vs the baseline at epoch 14) that had converged to the other models by
  epoch ~40–60.

Given that, and the frozen AlphaGenome encoder's speed, the focus moved to
efficiency; see the [U-Net note](unet-smf-body.md).

# Key Points

- **Cheri block in SMFNet (`smf_net.cheri.CheriBlock`):**
  - 3-tap dilated depthwise conv + per-read layer norm over the whole
    (L, C) plane (the fused Triton op), then a bias-free MLP
    (C → 2C → C, tanh-GELU), then `x + 0.15·MLP(...)`.
  - **Conditioning is adaLN-style:** `h·(1 + s) + b` between the norm and
    the MLP. `s` and `b` come from `LayerNorm(dna) → Linear(512 → 2C)`,
    zero-initialized, computed once per region and broadcast over its 32
    reads. This is SMFNet's existing FiLM formula, moved after the
    normalization.
  - **10 blocks with dilations 1–512,** to keep a receptive field of about
    2 kb with 3-tap kernels. The stem convs are followed by a GELU, as in
    Cherimoya.
- **Muon (`smf_net.muon.MuonAdamW`):**
  - One optimizer: Muon for the Cheri blocks' `linear1`/`linear2` weights
    (lr 0.025, decoupled weight decay 0.03, Nesterov momentum 0.95, 5
    Newton–Schulz steps, lr × √max(1, rows/cols), all following
    `torch.optim.Muon` in 2.9.1), and AdamW for everything else with
    SMFNet's existing rules.
  - **FiLM projections stay on AdamW,** because zero-initialized
    modulation and Muon's fixed-size updates seemed a risky mix. That was
    judged, not tested.
  - **Why one combined optimizer:** Lightning's automatic optimization,
    gradient clipping and `LinearWarmupCosineAnnealingLR` work unchanged.
    The scheduler's scalar warmup start and floor (lr/100 = 1e-5) apply to
    the Muon group too.
- **Why the code was copied in rather than installed:**
  - Cherimoya requires torch ≥ 2.9 and Triton ≥ 3.5.1, plus heavy
    dependencies (macs3, modisco, …).
  - SMFNet pins torch 2.5.1. conda-forge has torch 2.9.1 builds for
    glibc 2.17, but they need CUDA 12.9, while the `sh04-03n07` driver
    (550.163.01) reports CUDA 12.4.
  - `cheri.py`'s training kernels use only basic Triton operations and
    worked unchanged under Triton 3.1.
- **Two PyTorch 2.5 workarounds:**
  - The fused op runs under `torch.compiler.disable`: 2.5 does not trace
    these autotuned kernels.
  - The residual scale is applied *before* `linear2`, an equivalent form.
    Multiplying `linear2`'s output by a Python float matched inductor's
    weight-only-quantization pattern and crashed compiled `no_grad`
    forwards (validation and sampling) with
    `AttributeError: 'float' object has no attribute 'meta'`
    (`torch/_inductor/fx_passes/quantization.py`).

# Details

## Kernel check (GPU, Triton 3.1)

| dtype | shape (N, L, C), dilation | forward max abs | dx rel | dw rel |
| --- | --- | --- | --- | --- |
| fp32 | (64, 1024, 128), 1 | 1.4e-6 | 1.1e-7 | 3.3e-7 |
| fp32 | (64, 1024, 128), 512 | 1.4e-6 | 9.4e-8 | 2.9e-7 |
| fp32 | (8, 2048, 512), 64 | 2.9e-6 | 1.3e-7 | 2.4e-7 |
| bf16 | same three | 3–4e-2 | ~2e-5 | ~1.7e-3 |

The bf16 differences are at the level of bf16 rounding of normalized
values.

## Memory and speed (full batch: 256 regions × 32 reads, bf16, compiled, H100)

| variant | params | step | peak memory |
| --- | --- | --- | --- |
| baseline CNN (AdamW), 256 regions | 21.2M | 701 ms | 63.3 GiB |
| Cheri + Muon, 256 regions | 12.5M | out of memory | — |
| Cheri + Muon, 128 regions | 12.5M | out of memory | — |
| Cheri + Muon + per-block checkpointing, 256 regions | 12.5M | 1,367 ms | 48.9 GiB |

In real training the run took about 12 minutes per epoch, against the
baseline's 5.9. Reasons, mostly inferred rather than profiled:

1. **Checkpointing** recomputes each block's forward (≈ +33%).
2. **Cherimoya's speed claims are mainly for inference** with its
   megakernel, which isn't used here.
3. **Fewer FLOPs, more memory traffic.** A Cheri block costs about 8C² FLOPs
   per position against about 18C² for a dense k=9 conv, but it makes many
   bandwidth-limited passes: norm statistics, normalization, FiLM, GELU on
   2C, transposes, and an fp32 scratch buffer in the fused backward.
4. **The comparison baseline differs.** Cherimoya's gains are against
   ChromBPNet's full-width convs; SMFNet's per-read body is already narrow
   (128 channels).
5. **Smaller costs:** 10 blocks instead of 8, graph breaks around the fused
   op, and kernels tuned for newer Triton.

Fewer reads per region would cut memory roughly proportionally, but it
changes the training signal. Checkpointing at the full batch was used to
keep the comparison with the baseline exact.

## Runs and the Muon learning rate

Both runs: Cheri DNA encoder (the CNN type, not AlphaGenome) + Cheri SMF
body, `cheri_checkpoint=true`, `optimizer=muon`, `log_mask_rate_bins=true`,
`checkpoint_every_n_epochs=25`, and otherwise the baseline's
data and training settings (no EMA).

- **Muon lr 0.025:** `ONT_mESC_cheri_muon`, wandb `g9gh7vni`, ran in job
  50721347 ([scripts/train_cheri_muon.sh](cheri-muon-smfnet/scripts/train_cheri_muon.sh)).
  Interrupted (SIGINT) in epoch 37; checkpoints for epochs 24 and 35 are
  kept.
- **Muon lr 0.01:** `ONT_mESC_cheri_muon_lr0.01`, job 50810614
  ([scripts/train_cheri_muon_lr0.01.sbatch](cheri-muon-smfnet/scripts/train_cheri_muon_lr0.01.sbatch)),
  still running at the time of writing.

Validation loss (both region sets). The baseline's curve comes from its
local wandb file (run `3ekr0isw`), read offline with
[scripts/early_curves.py](cheri-muon-smfnet/scripts/early_curves.py):

| epoch | baseline | AG frozen | ATAC | Cheri lr 0.025 | Cheri lr 0.01 |
| --- | --- | --- | --- | --- | --- |
| 4 | 0.3126 | 0.3127 | 0.3123 | 0.3082 | 0.3112 |
| 8 | 0.3100 | 0.3090 | 0.3086 | 0.3068 | 0.3053 |
| 12 | 0.3092 | 0.3076 | 0.3055 | 0.3047 | 0.3031 |
| 14 | 0.3062 | 0.3070 | 0.3048 | 0.3047 | 0.3025 |
| 18 | 0.3037 | 0.3058 | 0.3032 | 0.3042 | 0.3013 |
| 30 | 0.3019 | 0.3033 | 0.3012 | 0.3032 | 0.3003 |
| 40 | 0.3003 | 0.3011 | 0.2998 | — | 0.2994 |
| 60 | 0.2988 | 0.3000 | 0.2978 | — | 0.2993 |

- **At lr 0.025,** training loss stalled (0.3040 → 0.3046) over epochs 6–9,
  exactly as Muon's lr warmed up from 0.0175 to 0.025. The val/loss lead
  then disappeared, and with cosine decay over 400 epochs the rate stays
  near its peak for a long time. AdamW-only runs trained smoothly over the
  same epochs.
- **At lr 0.01** there was no stall and a clear early lead. By epoch 60
  it's level with the AlphaGenome models and slightly behind ATAC.
- These are single seeds, and one-epoch differences are about 0.0005–0.001.
- EMA (Cherimoya uses decay 0.999) was not used, for comparability, and was
  deferred.

Commands, from the `nd_cheri` worktree with the main checkout's
environment (`module load devel pixi/0.53.0`,
`PIXI_CACHE_DIR=$SCRATCH/pixi-cache`):

```bash
cd ~/repos/SMFNet_cheri
PYTHONPATH=$PWD/src pixi run --frozen --manifest-path ~/repos/SMFNet_public/pyproject.toml pytest -q tests  # incl. GPU kernel test
cd $SCRATCH/smfnet_cheri && PYTHONPATH=~/repos/SMFNet_cheri/src pixi run --frozen \
    --manifest-path ~/repos/SMFNet_public/pyproject.toml python bench_cheri.py
sbatch $SCRATCH/smfnet_cheri/train_cheri_muon_lr0.01.sbatch
```

# Related Notes

- [U-Net SMF body](unet-smf-body.md): the efficiency-focused follow-up. It
  reuses Cherimoya's 0.15 residual scale as an initialization.
- [Frozen AlphaGenome encoder](../2026-10-29/frozen-alphagenome-encoder.md):
  the comparison models (baseline, AlphaGenome frozen and fine-tuned, ATAC),
  and the conditioning-placement reasoning reused for adaLN-style FiLM.
- [GCH regularization, EMA and sampling](../2026-10-07/regularization-ema-sampling-experiments.md):
  EMA 0.999 helped most on the GCH proxy. It's relevant to the deferred
  Cheri EMA test.

# Open Questions

- Does Cheri at lr 0.01 keep any advantage by the end of the run, on
  sampled profiles and covariation as well as val/loss?
- With EMA 0.999 (Cherimoya's recipe), does lr 0.025 work? That needs an EMA
  baseline too.
- Would Muon on the FiLM projections help or destabilize?
- How much of the 2× slowdown is checkpointing, the missing megakernel, or
  bandwidth-limited operations? A one-step profile is untested.
- Would torch ≥ 2.9 (compiling through the fused op) close the speed gap on
  a node whose driver supports CUDA 12.9?

# Sources

- Cherimoya: https://github.com/jmschrei/cherimoya (commit 5729b57), MIT
  License; `README.md`, `cherimoya/cheri.py`, `cherimoya/training.py`
  (`_split_parameters`, `configure_optimizers`).
- PyTorch 2.9.1 `torch/optim/_muon.py`:
  https://github.com/pytorch/pytorch/blob/v2.9.1/torch/optim/_muon.py
- SMFNet_public `nd_cheri`: `src/smf_net/cheri.py`, `src/smf_net/muon.py`,
  `src/smf_net/cnn.py`, `src/smf_net/model.py`, `tests/test_cheri.py`.
- [scripts/](cheri-muon-smfnet/scripts/): `bench_cheri.py`,
  `train_cheri_muon.sh`, `train_cheri_muon_lr0.01.sbatch`, `early_curves.py`.
