---
title: A frozen AlphaGenome encoder prior improves likelihood and the sequence prior, but not posterior ChIP agreement
date: 2026-11-06
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - smf_hmm repo, branch nd/experiments, commits 4c239f7 (AlphaGenome encoder and prior) and babfb52 (LayerNorm, clipping, weight EMA, bf16)
  - AlphaGenome weights google/alphagenome-all-folds (Hugging Face), downloaded to /scratch/users/diamant/models/alphagenome-all-folds; code google-deepmind/alphagenome_research (cloned at /home/users/diamant/repos/alphagenome_research, commit 3eecf48)
  - Runs /scratch/users/diamant/2026-11-05_fixed_blocks_global_kl3_map_targets_alphagenome_clip_ema_bf16 (in this allocation) and /scratch/users/diamant/2026-11-05_fixed_blocks_global_kl3_map_targets_clip_ema (Slurm job 51281138, owners, L40S)
  - Stopped runs ..._alphagenome (unnormalized features) and ..._alphagenome_ln (diverged)
  - Analyses /scratch/users/diamant/mesc_validation/covariation_windows/{,footprint_}{chr8,chr16}-512-valid90-reads10/analysis_clip_ema; plots .../covariation_windows/comparative_posteriors_2026-11-06/
  - Precision benchmark /scratch/users/diamant/alphagenome_precision_benchmark (Slurm job 51262269)
  - SMFNet's AlphaGenome encoder, /home/users/diamant/repos/SMFNet_unet (branch nd_unet), for the LayerNorm approach
  - User discussion, 2026-11-05 and 2026-11-06
tags: [smf-hmm, alphagenome, cnn-prior, sequence-prior, transfer, tf-chip, ctcf, nfy, oct4, training-stability, global-normalization, bf16, ema]
---

# Summary

To give the HMM's sequence prior more power, we replaced the prior's CNN trunk with the frozen
early layers of AlphaGenome: the DNA embedder and the first 3 residual blocks. They run at 1 bp
through an "à trous" rewrite that removes the max pools. Per-scale LayerNorm and trainable
adapters follow, then the existing conv layers and heads. Everything else matches the
lower-TF-prior `neural_global` configuration (`2026-11-04_..._map_targets`). A CNN control was
trained with the same gradient clipping and weight EMA.

**Results:**
- **Likelihood:** validation NLL per read 12.208 with the AlphaGenome prior, against 12.360 for
  the CNN control and 12.351 for the original CNN run. The 0.15-nat (1.2%) gain comes from the
  features, since clip and EMA alone change nothing.
- **Sequence prior:** the AlphaGenome prior alone tracks ChIP better than the CNN prior for
  NF-Y (+0.10, both chromosomes) and OCT4 (+0.06), and the same for CTCF. Metric: protected-read
  partial Spearman given accessibility.
- **Posterior ChIP agreement:** unchanged. AlphaGenome vs CNN: CTCF 0.64 vs 0.66, NF-Y 0.31 vs
  0.30, OCT4 0.17 vs 0.20. For NF-Y and OCT4, the AlphaGenome *prior alone* beats every
  posterior, including its own.

**Getting it to train** took three fixes:
1. **LayerNorm before the adapters.** AlphaGenome features have RMS about 4.5 with large
   per-channel offsets.
2. **Gradient clipping and weight EMA**, after a run whose globally normalised prior potentials
   blew up to about 3 × 10⁵. That made the float32 likelihood a cancellation artefact that the
   optimiser then exploited.
3. **A bf16 frozen encoder**, which is 3.9× faster.

# Key Points

- **The port matches DeepMind's Haiku modules to float32 rounding** (relative error about
  1e-7 at every scale) on real mm10 sequence with the real weights
  (`evaluation/alphagenome/check_encoder_port.py`).
- **Validation NLL and checkpoint health:**

| run | validation NLL | max \|prior potential\| | log Z (median) |
|---|---|---|---|
| CNN, lower TF prior (2026-11-04) | 12.351 | 27 | 26 |
| CNN + clip/EMA (control) | 12.360 | 21 | 16 |
| **AlphaGenome + clip/EMA, bf16 encoder** | **12.208** | 41 | 19 |
| AlphaGenome + LayerNorm, no clip (diverged, stopped) | (−9.9 log p, artefact) | ~3 × 10⁵ | ~1.1 × 10⁶ |

- **ChIP agreement:** protected-read partial Spearman given accessibility, mean of chr8 and
  chr16, with chr8 / chr16 in brackets:

| | CTCF | NF-Y | OCT4 | MAX::MYC | ZFP281 | ESRRB | KLF4 |
|---|---|---|---|---|---|---|---|
| phase_type | 0.60 | 0.21 | 0.10 | 0.04 | 0.06 | 0.04 | 0.02 |
| CNN + clip/EMA | 0.66 (0.66/0.66) | 0.30 (0.45/0.15) | 0.20 (0.14/0.25) | 0.04 | 0.06 | 0.03 | 0.01 |
| AlphaGenome + clip/EMA | 0.64 (0.64/0.64) | 0.31 (0.45/0.17) | 0.17 (0.12/0.22) | 0.06 | 0.05 | 0.04 | 0.00 |
| CNN + clip/EMA, prior only | 0.68 | 0.31 (0.43/0.19) | 0.19 (0.17/0.22) | 0.07 | 0.02 | −0.08 | −0.01 |
| AlphaGenome + clip/EMA, prior only | 0.67 | **0.41** (0.52/0.30) | **0.25** (0.22/0.27) | 0.09 | 0.02 | −0.03 | −0.04 |
| motif score | 0.54 | 0.17 | 0.38 | 0.12 | 0.15 | 0.20 | 0.19 |

  The raw and all-site partial metrics show the same pattern (in `analysis_clip_ema`).
- **Comparative plots:** where the two kinds of prior disagree, AlphaGenome changes how
  protected stretches are split.
  - At a high-ChIP NF-Y site (`chr8:75033999`), it calls TF blocks (motif P(TF) 0.56) where
    both CNNs call `other_long` (0.24–0.26).
  - At a CTCF site (`chr8:114541315`), it prefers open (0.31) where the CNNs call
    `TF_long_strong` (0.62–0.63).
  - The CNN with and without clip/EMA is nearly identical everywhere.

# Details

**AlphaGenome encoder** (`smf_hmm/alphagenome_encoder.py`; architecture from
`alphagenome_research/model/{model,convolutions,layers}.py`):
- **Architecture:**
  - `DnaEmbedder`: a kernel-15 conv to 768 channels plus a residual kernel-5 `ConvBlock`
    (RMSBatchNorm, then GELU, then a weight-standardised conv).
  - Then `DownResBlock`s at bin sizes 2, 4, 8, …, each after a 2× max pool and each adding 128
    channels.
  - Only the embedder is at 1 bp. The organism embedding enters later, at the 128 bp trunk, so
    the encoder is organism-agnostic.
- **À trous rewrite:**
  - each pool becomes `max(x[t], x[t + b/2])`, and each block's convolutions are dilated by its
    bin size b;
  - output t then equals the pooled network's output for the bin starting at t;
  - each scale is shifted right by ⌊(b − 1)/2⌋ so it is centred on t.
- **Receptive field** (estimated from kernels and dilations): ±9, ±18, ±35 and ±69 bp for the
  embedder and blocks 0–2.
- **Weights:**
  - read straight from the Orbax OCDBT/zarr checkpoint with TensorStore, since our Orbax 0.4.4
    can't restore the newer layout;
  - cast to float32 (some are stored in bf16);
  - held as non-trainable `nnx.Variable`s and saved in our checkpoints, about 63 MB.
- **Size:** 32.8 M frozen parameters.

**Prior** (`AlphaGenomePriorArchitecture`, a subclass of `SimpleCNNPriorArchitecture`):
- the frozen encoder (3 blocks), then per scale `nnx.LayerNorm` (ε = 1e-5) and a linear layer to
  64 channels;
- the scales are concatenated (256 channels), then 2 conv layers (128 channels, kernel 9, ReLU);
- then the existing state-logit and emission heads, zero-initialised for fixed transitions.
- Config: `prior.architecture=alphagenome`, `prior.alphagenome_checkpoint`,
  `prior.alphagenome_blocks`, `prior.alphagenome_adapter_features`,
  `prior.alphagenome_compute_dtype`.

**Training stability:**
- **Run 1 (no normalisation):** losses ran 0.1–0.2 above the CNN at matched steps.
  - AlphaGenome feature RMS by scale is 4.8, 4.7, 4.3 and 4.6, with per-channel std of only
    about 0.5–1 (large offsets) and outliers near 90, matching SMFNet's measurements.
  - Each adapter sums about 1,000 such inputs, so Adam's effective step was hundreds of times
    that of the CNN, whose first layer sees one-hot input.
  - Stopped.
- **Run 2 (per-scale LayerNorm, as in SMFNet's `CNNLayerNorm1D`):** trained 0.12–0.17 *below*
  the CNN at matched steps.
  - At step 7,800 a single loss spike reached 58.7, against 11.1 for the CNN on the same batch.
  - Afterwards the training loss first sat about 2 nats above the CNN, then implausibly 2.5
    nats *below* it (8.8 at step 10,000), then 20–26.
  - The epoch-23 "best" checkpoint had prior potentials of about 3 × 10⁵ and log Z of about
    1.1 × 10⁶. With global normalisation the likelihood is a difference of numbers around 10⁶
    in float32, a resolution of about 0.1, so the "improvement" was cancellation error.
  - Training used plain `optax.adam` with no clipping. Stopped.
- **Run 3:** `optimizer.clip_norm=1.0` (`optax.clip_by_global_norm` before Adam),
  `optimizer.ema_decay=0.999` (a weight EMA used for validation and checkpoints) and the bf16
  encoder.
  - No spikes, healthy potentials, still improving near epoch 100.
  - The user chose clipping without a bound on the potentials. A soft bound,
    `B · tanh(x / B)`, was considered but not implemented.

**Precision** (Slurm job 51262269, L40S; encoder only, 128 × 512 bp, 3 blocks):

| setting | seconds per batch | error vs float32 |
|---|---|---|
| float32 | 0.163 | — |
| TF32 | 0.070 | RMS relative 3e-5 to 1.1e-4 |
| bfloat16 | 0.042 | RMS relative 2.8e-3 to 4.8e-3; correlation after LayerNorm ≥ 0.99999 |

- Training ran at 3.8 steps/s with the bf16 encoder, against 2.6 with float32.
- The forward–backward pass through the trainable prior and the HMM stays float32.
- bf16 autotuning logs 12 harmless "Results mismatch between different convolution algorithms"
  warnings.

**Interpretation (inferred):**
- **The features help the sequence model.** AlphaGenome's features give a better model of the
  reads (NLL) and a better sequence-only predictor of NF-Y and OCT4 binding.
- **The posterior doesn't keep that.** On protected reads, P(TF) is decided by "TF footprint
  versus nucleosome or `other_long`", and that read evidence ranks sites differently from ChIP.
  For NF-Y and OCT4 it even pulls the ranking away from the better prior.
- **The limit may now be the read model.** Candidates: the HMM's TF-versus-other-protection
  emissions or block structure, rather than the sequence model.
- **The likelihood gain may matter elsewhere.** Better per-molecule footprint calls, which ChIP
  can't score, could account for some of it.

**Licence:** AlphaGenome weights fall under the
[model terms](https://deepmind.google.com/science/alphagenome/model-terms): non-commercial, with
restrictions on derived models and on training with outputs. This is an unpublished early
experiment, and the user is contacting the AlphaGenome team about the use case.

# Commands

Runs used interactive Slurm job 51132065 (sh04-06n09, 1 L40S) plus owners-partition jobs
51262269 (benchmark) and 51281138 (CNN control).

```bash
cd /home/users/diamant/repos/smf_hmm
# Weights (gated; needs a Hugging Face token with access to public gated repos)
curl -sfL -H "Authorization: Bearer $(cat ~/.cache/huggingface/token)" -o <file> \
  https://huggingface.co/google/alphagenome-all-folds/resolve/main/<file>   # each checkpoint file, into /scratch/users/diamant/models/alphagenome-all-folds
# Port check against DeepMind's Haiku modules (Haiku etc. in a --no-deps --target dir, not the hmm env)
JAX_PLATFORMS=cpu PYTHONPATH=/home/groups/btrippe/diamant/pylibs/alphagenome_check:$PWD \
  python evaluation/alphagenome/check_encoder_port.py
sbatch /scratch/users/diamant/alphagenome_precision_benchmark/benchmark.sbatch                    # precision benchmark
bash /scratch/users/diamant/2026-11-05_fixed_blocks_global_kl3_map_targets_alphagenome_clip_ema_bf16/pipeline.sh   # train + posteriors
sbatch /scratch/users/diamant/2026-11-05_fixed_blocks_global_kl3_map_targets_clip_ema/train_and_infer.sbatch      # CNN control
bash /scratch/users/diamant/mesc_validation/covariation_windows/analysis_clip_ema_2026-11-06.sh    # combined analysis
```

- Comparative plots: `plot_posterior_examples` with `--site-key` on the 3 largest and 3
  smallest AlphaGenome − CNN gaps in motif P(TF), among bound sites with ≥ 10 covered reads
  (keys in `comparative_posteriors_2026-11-06/*_site_keys.txt`).
- The checkpoint health checks (potentials, log Z, mean log p on 40 val windows), the feature
  scale measurement and the site selection were one-off snippets.

# Related Notes

- [Controlled ChIP metrics and occupancy targets](../2026-11-04/controlled-chip-metrics-and-context-offsets.md): The lower-TF-prior (MAP-target) configuration this experiment builds on.
- [Protected-read ChIP metric](../2026-11-05/protected-read-chip-metric.md): The headline metric used here.
- [MAX::MYC, NF-Y and ZFP281](../2026-11-05/chip-nexus-tfs-max-myc-nfy-zfp281.md): The footprint benchmark datasets and why CTCF, NF-Y and OCT4 are the focus.

# Open Questions

- Why does combining a better prior with the reads lower NF-Y and OCT4 site agreement? Would TF
  block types or emissions that better separate TF from nucleosome or `other_long` keep the
  prior's gain?
- Does the likelihood gain correspond to better per-molecule calls, e.g. measured by posterior
  predictive checks?
- Would more AlphaGenome blocks (more context), or fine-tuning the last block, help further?
- Should the globally normalised prior get a soft bound on its potentials to rule out the
  cancellation exploit for good?

# Sources

- Code: smf_hmm `nd/experiments` `4c239f7`, `babfb52` (`smf_hmm/alphagenome_encoder.py`, `smf_hmm/prior_architectures.py`, `smf_hmm/train.py`, `smf_hmm/cli/train.py`, `evaluation/alphagenome/`)
- AlphaGenome: https://github.com/google-deepmind/alphagenome_research, https://huggingface.co/google/alphagenome-all-folds, [Nature paper](https://www.nature.com/articles/s41586-025-10014-0)
- Runs and analyses: paths in `sources` above
