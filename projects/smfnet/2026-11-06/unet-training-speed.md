---
title: Attention U-Net training step 292 → 184 ms with a fused Triton norm kernel, channels-last activations and a Toeplitz output conv (H100)
date: 2026-11-06
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - SMFNet_public branch nd_unet_fast (worktree ~/repos/SMFNet_fast; not pushed), commits 98416ce, ecadace, c4641d9, 6b1707b, f85020f
  - src/smf_net/fused_norm.py, src/smf_net/channels_last.py, src/smf_net/cnn.py, src/smf_net/alphagenome.py, tests/test_fast_unet.py
  - unet-training-speed/scripts/ (benchmark and profiling scripts)
tags: [efficiency, triton, kernel-fusion, channels-last, unet, profiling, h100]
---

# Summary

The model sped up here is the attention U-Net SMFNet: a frozen AlphaGenome
encoder, per-read norm, 2 bottleneck attention layers and cross-attention
(see [the U-Net note](../2026-11-03/unet-smf-body.md)). Its training step was
cut from 292 ms and 56 GiB to 184 ms and 38 GiB. Measurements are at full
batch (256 regions × 32 reads, 2,048 bp DNA, 1,024 bp SMF window), bf16,
compiled, on one H100.

The changes are exact rewrites:

- outputs match the original in fp32 to ~1e-6;
- a smoke training run follows the reference loss trajectory.

All changes are behind `CNNConfig` options, off by default:
`unet_fused_norm`, `unet_channels_last` and `output_conv_matmul`. Old
checkpoints load unchanged.

In real Lightning training on the full data, the steady-state step is
188–190 ms. The same configuration took ~300 ms before the changes, so the
real-training speedup is about 1.6×.

# Key Points

- **Cumulative effect** (plain training loop: forward, backward, gradient
  clipping, AdamW; `cudnn.benchmark=True`):

  | change (cumulative) | step | peak memory |
  | --- | --- | --- |
  | none (current attention U-Net) | 292 ms | 56 GiB |
  | fused per-read norm → affine → FiLM → GELU Triton kernel | 281 ms | 33 GiB |
  | channels-last (N, L, C) activations in the U-Net | 256 ms | |
  | output conv (k=75, 1 output channel) as a matmul | 234 ms | |
  | … as blockwise Toeplitz matmuls instead | 216 ms | |
  | native column sum for channels-last bias gradients | 193 ms | |
  | channels-last AlphaGenome encoder | 183 ms | 38 GiB |
  | single-pass fused-norm forward (sustained 200-step loop) | 186 → 184 ms | |

- **Real training matches the plain loop.** A Lightning callback (below)
  profiled steps 60–69 of a real run:
  - CUDA kernel time was 184 ms per step;
  - host time per step was 188–190 ms;
  - CUDA sync-debug mode found no host–device syncs in the step.

  An earlier estimate of 220–230 ms per step came from tqdm's cumulative it/s,
  which includes start-up. Differencing tqdm timestamps (steps 112 → 192)
  gives ~200 ms for that commit (before the channels-last encoder) and
  ~300 ms for the reference. Validation and DataLoader start-up add a
  per-epoch cost on top; that cost was not measured separately.
- **Where the remaining 184 ms goes** (compiled-step profile): convs and
  matmuls ~69 ms, the fused norm ~40 ms (backward ~27 ms), inductor pointwise
  kernels ~36 ms, reductions ~21 ms, other ~13 ms. cuDNN's layout transposes
  are gone. Most remaining kernels are bandwidth-bound. Further large gains
  would probably need model changes (for example, a narrower level 0 or a
  cheaper norm) rather than rewrites. This is inferred, not tested.

# Details

## 1. Fused per-read norm kernel (`src/smf_net/fused_norm.py`)

- **What it computes.** `unet_norm=read` is a GroupNorm(1, C) per read,
  followed by the channel affine, the FiLM `h·(1 + s) + t`, and AlphaGenome's
  GELU `x·sigmoid(1.702x)`.
- **Why fuse it.** Compiled PyTorch runs this as several kernels and saves
  intermediates for the backward pass.
- **The Triton version** is registered as `torch.library.custom_op`s
  (`smfnet::read_norm_film_gelu_fwd` / `_bwd`), so torch.compile treats it as
  opaque. It supports both layouts:
  - The forward stores only per-read mean and rstd.
  - The backward recomputes the normalized values. Each program loops over a
    region's R reads, so the FiLM gradients (dst) and dγ/dβ accumulate
    without atomics.
- **Peak memory.** Most of the saving (56 → 33 GiB) comes from not saving
  intermediates.
- **Single-pass forward** (f85020f). When there are at least
  `SINGLE_PASS_MIN_READS = 2048` reads, one program per read computes the
  statistics and then normalizes. The second read of x then comes mostly
  from L2.
  - Forward time per level fell 18–50%, for example 1.06 → 0.87 ms at
    L=1024, C=48, and 0.12 → 0.06 ms at L=16, C=144.
  - Block 2048 with 8 warps is used for reads of at least 16k elements,
    otherwise 1024 with 4 warps.
  - Outputs agree to 1 bf16 ulp.
  - With fewer reads (for example a few long 10 kb reads), the chunked
    two-kernel forward is kept, because one program per read would not fill
    the GPU.
  - Only 2 ms per step overall, because the backward dominates.

## 2. Channels-last activations (`src/smf_net/channels_last.py`, `cnn.py`)

- **The problem.** cuDNN's fastest bf16 kernels on Hopper are NHWC. For
  (N, C, L) inputs it inserted NCHW↔NHWC transposes around every conv, about
  a fifth of the step.
- **The fix.** The U-Net, and in 6b1707b the AlphaGenome encoder, keep
  (N, L, C) activations. `conv_channels_last` runs each conv as `conv2d` on a
  zero-copy (N, C, 1, L) `channels_last` view; k=1 convs are `F.linear`.
- **What didn't work.** Swapping conv1d for conv2d alone, without the
  channels-last layout, gave no gain. The earlier observation that conv2d is
  faster in the original SMFNet holds only together with the layout.
- **Bias gradients.** Over (N, L, C), these are column sums. Inductor's
  compiled version (PyTorch 2.5) reads with a stride and was ~15× slower. The
  bias add is an `autograd.Function` whose backward calls the native sum
  through the opaque op `smfnet::channel_sum` (c4641d9).
- **Average pooling and upsampling** become `view` + `mean` and
  `repeat_interleave` on dim 1.

## 3. One-output-channel output conv (`conv_to_one_channel`, ecadace)

- cuDNN is very slow on the 48 → 1 channel, k=75 conv.
- A matmul formulation brought the step from 256 to 234 ms.
- The blockwise Toeplitz version brought it to 216 ms. It builds a
  (C·(block+k−1)) × block Toeplitz matrix from the weight (block 256), so the
  conv is one dense GEMM over 256-position tiles, about 10× faster than cuDNN
  for this layer.
- A first attempt as a direct Triton conv took 56 ms and was dropped.

## 4. Pitfalls found on the way (PyTorch 2.5)

- **`zip` over unequal-length lists in `forward`** caused graph breaks. With
  them the cache limit of 8 was hit, so the model fell back to eager, with 2×
  time and memory and an OOM. Fixed by slicing the lists to equal length.
- **An mm followed by a float scale** crashed an inductor quantization pattern.
- **Concurrent GPU jobs fail.** Sherlock btrippe GPUs are exclusive-process,
  so a second job gets "CUDA-capable device(s) is/are busy". Benchmarks must
  run one at a time.

## 5. Rejected or not worth it

- **Triton column sum for bias gradients** (`bench_channel_sum.py`).
  - The native sum is within 1.5–2× of the bandwidth ideal at the large
    levels.
  - The Triton kernel was slower at several mid-size levels.
  - At most ~1–2 ms total could be saved.
- **Other host-side causes, ruled out:** the progress bar refresh rate,
  `log_every_n_steps`, DataLoader throughput (25 ms per batch, below the step
  time), `get_loss` overhead (3 ms) and GPU clock throttling. The sustained
  loop stayed at 1815 MHz and ~640 W.

## Verification

- `tests/test_fast_unet.py`, GPU tests that are skipped without CUDA:
  - the fused kernel against a PyTorch reference in both layouts, with and
    without FiLM, and with both forward paths;
  - the fast U-Net against the original in fp32 with TF32 off, outputs and
    all parameter gradients;
  - `conv_to_one_channel` against conv1d, including partial blocks and the
    weight gradient;
  - the channels-last AlphaGenome encoder against the original.

  With `tests/test_unet.py`, 13 tests pass.
- A smoke training run with the options on matched the reference trajectory.

## Commands

The environment: the worktree runs with the main checkout's pixi env.

```bash
module load devel pixi/0.53.0; export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
cd ~/repos/SMFNet_fast
PYTHONPATH=$PWD/src pixi run --frozen --manifest-path ~/repos/SMFNet_public/pyproject.toml \
    python <script>
```

The scripts, all in [unet-training-speed/scripts/](unet-training-speed/scripts/),
originals in `$SCRATCH/smfnet_fast/`:

| script | what it measures |
| --- | --- |
| `bench_all_model.py [--cudnn-benchmark]` | full-model equivalence and step time |
| `sustained_loop.py` | 200-step loop with clock and power readings |
| `bench_fused_norm.py`, `bench_single_pass.py` | per-level kernel timings |
| `toeplitz_bench.py` | the output conv |
| `bench_channel_sum.py` | bias-gradient column sums |
| `profile_triton_names.py compiled --fused --cl` | per-kernel profile |

`prof_cb.py` profiles real Lightning training. It was loaded through a
temporary 6-line hook in `scripts/train.py` keyed on the environment
variable `SMFNET_PROFILE_STEPS`. The hook was reverted and is not committed.

To turn the options on in training:

```bash
python scripts/train.py ... model.cfg.cnn_config.unet_fused_norm=true \
    model.cfg.cnn_config.unet_channels_last=true model.cfg.cnn_config.output_conv_matmul=true
```

# Related Notes

- [U-Net SMF body](../2026-11-03/unet-smf-body.md): the architecture sped up
  here, and its original step time (329 ms before the attention, Muon and
  cudnn.benchmark changes).
- [Cheri blocks and Muon in SMFNet](../2026-11-03/cheri-muon-smfnet.md): the
  first fused Triton kernel on PyTorch 2.5. The custom-op approach and the
  per-region read loop in the backward follow from that work.
- [Frozen AlphaGenome encoder](../2026-10-29/frozen-alphagenome-encoder.md):
  the encoder now also run channels-last.

# Open Questions

- Is the fused-norm backward (~27 ms per step) close to its bandwidth bound?
  A per-read single-pass backward conflicts with accumulating FiLM gradients
  across a region's reads.
- How do the options scale to 10 kb windows with fewer reads per region?
  The single-pass forward switches off below 2,048 reads; the Toeplitz block
  size and memory headroom are untested there.
- Should the freed memory (56 → 38 GiB) go to a wider level 0 or more reads
  per region?
- Is the per-epoch overhead (validation, DataLoader start-up) significant
  relative to ~1,800 steps per epoch?

# Sources

- SMFNet_public `nd_unet_fast` commits 98416ce, ecadace, c4641d9, 6b1707b,
  f85020f (worktree `~/repos/SMFNet_fast`, not pushed).
- Benchmark and profiling scripts in
  [unet-training-speed/scripts/](unet-training-speed/scripts/). Logs are in
  `$SCRATCH/smfnet_fast/`: `lightning_torchprof.log`, `lightning_sync.log`,
  `throughput_ref.log`, `throughput_fast.log`.
