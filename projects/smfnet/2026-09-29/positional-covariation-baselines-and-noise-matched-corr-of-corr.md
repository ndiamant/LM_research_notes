---
title: Train-set positional baselines and noise-matched correlation-of-correlations
date: 2026-09-29
project: smfnet
agent: Claude Code
status: draft
sources:
  - SMFNet repo, scripts/explore_positional_baseline.py (branch nd/more_assays, uncommitted)
  - SMFNet repo, src/smf_net/evaluation/{metrics,covariation,summary,metrics_plotting}.py
  - /scratch/users/diamant/positional_baseline_explore/ (summary.json, *_bandwidths.tsv, test_preview.tsv, *.png, baselines.npz)
  - /scratch/users/diamant/evaluate_SMFNet_test_set/downsample_025/{sample/samples.h5,metrics/metrics.tsv}
  - User and advisor discussion, 2026-09-23 to 2026-09-29
tags: [smfnet, evaluation, baselines, covariation, correlation-of-correlations, pseudo-replicate, sampling-noise]
---

# Summary

We checked how much of SMFNet's correlation metrics a **sequence-blind train-set
baseline** can reach. There are two baselines: a positional mean profile and a
mean pair-correlation surface f(separation, midpoint offset). They were fit on
train chromosomes, bandwidths were chosen on val, and they were evaluated on the
`downsample_025` model's 13,744 test regions.

- **Profile.** The positional baseline is weak. Median per-region Pearson is
  0.29, against the model's 0.81, and the model wins in 98% of regions.
- **Covariation.** The baseline looked stronger than the model. Its
  corr-of-corr is 0.715, above both the model's metric value (0.603) and the
  pseudo-replicate (0.657). **This is a noise artifact, not a model failure.**
  The baseline is noise-free. The metric's model and replicate estimates each
  come from ~80 reads, so they pay a sampling-noise penalty the baseline
  doesn't.
- **Scored from all 1,024 of its reads** against the same observed sample, the
  model reaches **0.758** and beats the baseline in 85% of regions.
- **Residual corr-of-corr** subtracts the baseline from every pair correlation
  before correlating. There the model keeps **0.387**, about 70% of the
  noise-free ceiling. That is region-specific coupling the generic curve
  cannot explain.

Conclusion: the matched-read corr-of-corr is not broken, but it is dominated by
generic nucleosome phasing and read-count noise. The all-reads and residual
versions are the informative ones. None of this is in `metrics.py` yet.

# Key Points

- **Why a noise-free baseline wins.** Roughly,
  corr(x, y) ≈ corr(signal_x, signal_y) · √(rel_x · rel_y), where rel is the
  fraction of variance that is signal rather than sampling noise.
  - The replicate pays the noise penalty on both sides, so it scores ≈ rel_n.
  - The noise-free baseline pays it once, so it scores ≈ corr(f, ρ) · √rel_n.
  - Noise-free ceiling: √(replicate corr-of-corr) = √0.657 ≈ 0.81.
  - Attenuation predicted the all-reads model at 0.603 / √0.657 = 0.744;
    measured 0.758. That's close; the correction assumes homogeneous noise.
- **Fractions of ceiling** (Pearson): model all-reads 0.93, baseline 0.88,
  model matched 0.92 (against the replicate itself).
- **This is consistent with the downsampling curve.** The user reports a clean
  improvement from 1/8 → 1/4 → 1/2 → full training data on this metric. Every
  model is scored with the same n-read noise, so rankings still hold. All of
  them are just capped below a noise-free generic curve.
- **Bandwidths barely matter for the profile.** Val median Pearson is
  0.355–0.356 for 0–100 bp and drops only at 400 bp (0.333). With 55k train
  regions the raw per-bp mean is already smooth. Mirroring the profile changes
  nothing, and the profile is symmetric.
- **Covariation needs almost no separation smoothing.** 2 bp was best, the
  smallest tried. Midpoint dependence helps a little: pooled gives 0.721 and a
  50 bp midpoint bandwidth gives 0.730 on val.
- **Train pair correlation shows clear ~190 bp nucleosome phasing.** The
  phasing is stronger in the region flanks than at the centre. Mean call has
  a broad central dip (about 0.72 on the flanks, 0.55 at the centre). The
  regions are selected to contain TF motifs and to be highly active under a
  ChromBPNet model, so they are loosely anchored.
- **Pooling CpG and GpC calls is fine** under `cg_gc_accessibility`: their
  profiles match. GCG sites sit about 0.08 lower everywhere, a constant offset
  rather than a different shape.

# Details

## Earlier metric changes this builds on

These changes were committed by the user as `21152e9` "more sensible metric
baselines" on `nd/prism_validation`. They are **not** on `nd/more_assays`.

- **Corr-of-corr shuffle scalars dropped.** Both within-read and read-axis
  shuffles give corr-of-corr ≈ 0 by construction: one leaves pure noise, the
  other a constant offset. Their decay curves are kept. Corr-of-corr therefore
  has no gap-closed row.
- **Matched profile correlations added:** `pearson_matched`,
  `pearson_matched_pseudo_replicate` and `pearson_matched_shuffle`, plus the
  Spearman equivalents. The model sample, replicate half and target half are
  disjoint and equal in size, with n = min(N_pred, N_obs // 2) and no
  128-read cap.
  - These samples are drawn last from the generator, so earlier columns keep
    their values.
  - The old leaky `pearson_shuffle` was removed. It overlapped its target, and
    its read count differed from the model's.
- **The matched shuffle is ≈ 0 for these long reads:** median 0.039 against a
  replicate of 0.797. It stays in the TSV for possible short-read models but is
  hidden from the plots.
- **Plotting.** Profile panels overlay Model (all reads), Model (matched reads)
  and Pseudo-replicate (matched reads). Figures use `layout="constrained"`,
  which fixed a suptitle overlap.

## Test preview numbers (medians, 13,744 regions)

| | Profile Pearson | corr-of-corr | residual corr-of-corr |
|---|---|---|---|
| Model, matched reads (metric) | — | 0.603 | 0.203 |
| Model, all 1,024 reads | 0.812 | 0.758 | 0.387 |
| Pseudo-replicate | — | 0.657 | 0.302 |
| Positional baseline | 0.294 | 0.715 | 0 (by construction) |

- **The replay matches the metrics exactly.** The preview replays each region's
  `region_generator(seed=0, index)` draw of `_prepare_{mmd,matched}_samples`
  (max 128 reads), so the observed sample is the one the metrics used. The
  replayed read counts match `n_mmd_reads_per_sample` in 100% of regions.
- **Model reads get donated missingness.** All-reads and matched model samples
  receive per-read missingness donated from observed reads, from a separate
  generator. Re-donated matched scores (0.603) agree with the metric column
  (0.604).
- **Profile scoring.** The positional baseline is scored against the
  full-depth observed mean, like the headline `pearson`.

## Performance note

Per-region tensors are small: median 194 reads × 134 sites on val. torch's
default thread pool made the per-region loops 4–17× slower than
`torch.set_num_threads(1)`. The first full run was killed after ~40 min for
this reason. The script now pins one thread, and the full fit plus preview
takes about 20 minutes on 8 CPUs.

## Reproduce

These ran interactively on a Sherlock compute-node allocation (8 CPUs, 125 GB),
not through `sbatch`. Environment: `$GROUP_HOME/uv-envs/SMFNet-cu124`, CPU only.

```bash
cd ~/repos/SMFNet   # branch nd/more_assays; script is uncommitted
D=$SCRATCH/evaluate_SMFNet_test_set/downsample_025
O=$SCRATCH/positional_baseline_explore
CUDA_VISIBLE_DEVICES="" MPLBACKEND=Agg \
  $GROUP_HOME/uv-envs/SMFNet-cu124/bin/python -u scripts/explore_positional_baseline.py \
  --samples-path $D/sample/samples.h5 \
  --metrics-tsv $D/metrics/metrics.tsv \
  --output-dir $O > $O/run.log 2>&1
# Rerun only the test preview from the saved fit (skips the ~8 min train load):
#   ... same command ... --preview-only
# Smoke run: --max-regions-per-split 300 --max-test-regions 200
```

Defaults are ONT_mESC_acRegions.h5, mm10.fa, `cg_gc_accessibility`, and the
default chromosome splits (train 16 chromosomes, val chr8/chr16, test
chr1/chr3/chr6). Fitting uses full-depth train reads, not the downsampled
training fraction. The script's import accepts either helper name, so it runs
on both branches.

# Related Notes

- [Fixed-HMM posterior-sampling baseline comparison](../../smf-hmm/2026-09-04/fixed-hmm-posterior-sampling-comparison.md):
  The fixed geometric and phase-type HMMs in `smf_hmm/baselines/fixed_hmm` are
  close to the Markov-1 and semi-Markov generative baselines proposed below.
  They might be reusable as SMFNet read generators instead of writing new ones.

# Open Questions

- **Promote to `metrics.py`?** Candidates are all-reads model corr-of-corr and
  residual corr-of-corr. The residual version needs the train pair grid
  (`baselines.npz`) as an input.
- **Generative sequence-blind baselines.** A ladder of Markov-1
  P(call | previous call, gap), then a semi-Markov hazard model
  P(run continues | state, span so far, gap), then SMFNet. Output would be
  `samples.h5`, scored by `metrics.py` like any model, so every metric
  (including run-length Wasserstein) is noise-matched.
  - The hazard form handles runs censored by missing calls or window edges.
    A simple renewal sampler of train run spans would drop those runs, biasing
    it towards short runs.
  - A pooled train run-length distribution compared directly with test regions
    is **not** fair: run spans depend on each region's site layout and
    missingness.
  - Not implemented.
- **Profile baseline refinements** (speculative): add sequence context (site
  type, local CpG/GC density), or a GCG indicator. These might raise the 0.29
  profile baseline, but the profile gap to the model is already large.
- **Replicate the check** on the other downsample levels and the full model, to
  see whether the residual metric tracks training data as cleanly as the raw
  one.

# Sources

- Script: `SMFNet/scripts/explore_positional_baseline.py` (uncommitted, branch `nd/more_assays`).
- Outputs: `/scratch/users/diamant/positional_baseline_explore/`
  - `summary.json`: all medians, ceilings, fractions
  - `profile_bandwidths.tsv`, `covariation_bandwidths.tsv`: val bandwidth sweeps
  - `test_preview.tsv`: per-region baseline and model scores joined to `metrics.tsv`
  - `profile.png`, `decay_by_midpoint.png`, `test_preview.png`
  - `baselines.npz`: fitted profile and pair grid, with raw sums and counts
- Model evaluated: `downsample_025` checkpoint
  `/scratch/users/diamant/smf_models/SMFNet_downsample/downsample_025/epoch=47-step=21840.ckpt`,
  1,024 sampled reads per test region.
- Metric changes: SMFNet commit `21152e9` on `nd/prism_validation`.
