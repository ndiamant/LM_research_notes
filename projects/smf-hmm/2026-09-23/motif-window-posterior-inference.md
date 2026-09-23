---
title: Motif-window posterior inference pipeline
date: 2026-09-23
project: smf-hmm
agent: Codex
status: draft
sources:
  - /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/prepare_dataset.py
  - /home/users/diamant/repos/smf_hmm/smf_hmm/posterior_inference.py
  - /home/users/diamant/repos/smf_hmm/smf_hmm/cli/infer_posterior.py
  - /scratch/users/diamant/2026-09-08_fixed_phase_type_kernel_9/checkpoints
tags:
  - smf-hmm
  - tf-chip
  - posterior-marginals
  - fixed-hmm
  - inference
---

# Summary

A model-independent preparation step now selects motif windows using cached
ChIP and CpG/GpC SMF coverage, and a shared inference CLI writes posterior
marginals for neural or fixed HMMs. The design uses 512-bp motif-centered
windows, computes sequence-conditioned prior terms once per window, and
broadcasts them across read batches. It stores all named element probabilities
at every window base, including bases with missing SMF observations.

The new code passed synthetic tests. The supplied 138-state neural checkpoint
also loaded and completed a small synthetic GPU run. No real chromosome-scale
preparation or inference has been run yet.

# Key Points

- Site selection uses each TF/chromosome's cached pooled ChIP sum, with a
  configurable quantile (default 0.5) and strict `>` selection. The threshold
  is computed over all cached motif sites before SMF filtering; ties can make
  the selected fraction smaller than one half. This is an upper-half signal
  subset, not a matched genomic background or calibrated binding threshold.
- SMF rows are read from the TSV produced by `evaluation/datasets/smf/scripts/h5_to_calls.py`.
  Assayable positions are plus-strand C bases with G on either neighbor (CpG
  or GpC); a C satisfying both is counted once. Context outside the 512-bp
  window is used to classify edge cytosines.
- Default read eligibility requires the row interval to span the full window
  and nonzero calls at >=90% of its reference CpG/GpC positions. Default
  window eligibility requires at least 10 eligible rows. These thresholds are
  configurable and should be assessed from the emitted candidate summary.
- `prepare_dataset` writes `dataset.h5`, retained window metadata, all
  candidate selection diagnostics, per-TF thresholds, and provenance. It
  streams the SMF TSV once and appends selected observations to HDF5.
- `infer_posterior` supports a neural checkpoint and the fixed geometric and
  phase-type baselines using the same loop. Its `posterior` dataset has shape
  `(total window-read pairs, window_size, n_element_groups)` and float32 dtype.
  Window offsets, group names, molecule IDs, strands, metadata, and provenance
  accompany the array.
- FixedPhaseTypeTransition now applies destination bias to all allowed
  destinations in every row. Neural conditioning can therefore change
  terminal tail durations while deterministic chains still enforce minimum
  durations. The fixed baseline uses zero bias and retains its old transition
  matrix.

# Details

## Preparation

Example chr8 preparation, using the paths in
`evaluation/tf_chip/configs/mesc_validation.yaml`:

```bash
cd /home/users/diamant/repos/smf_hmm
PYTHON=/home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python
"$PYTHON" -m evaluation.tf_chip.prepare_dataset \
  --config evaluation/tf_chip/configs/mesc_validation.yaml \
  --cache-root /scratch/users/diamant/mesc_validation/benchmark_cache \
  --chrom chr8 --tf ctcf --tf oct4 --tf klf4 --tf esrrb \
  --window-size 512 --chip-quantile 0.5 \
  --min-valid-fraction 0.9 --min-reads 10 \
  --output-dir /scratch/users/diamant/mesc_validation/prepared/chr8-q50-valid90-reads10
```

Repeat with `--chrom chr16` and a distinct output directory. Existing output
directories are never replaced. The script verifies cache checksums and FASTA
identity, validates nonzero observations against reference CpG/GpC cytosines,
and publishes a complete output directory atomically. Empty selection is a
valid diagnostic result.

The input TSV was produced via `h5_to_calls.py`. That conversion uses stored
region bounds and synthetic `<region>:<row-index>` IDs, with strand `.`; these
are not original alignment spans and BAM molecule IDs. This preparation step
can ensure only that a stored TSV row covers a window. If overlapping source
regions encode the same physical molecule under different synthetic IDs, it
cannot identify or deduplicate those rows. Outputs preserve their input IDs.
This caveat applies to this TSV representation and should be revisited when a
source with original molecule identity and alignment bounds is used.

The output HDF5 has `/windows/<window_id>/dna`, `positions`, `observations`,
`molecule_id`, and `strand`. DNA uses A/C/G/T = 0/1/2/3; positions are sparse
window-relative assayable cytosine offsets; observations at these offsets use
-1 missing, 0 accessible, 1 protected. Dense model observations are made by
filling a `(reads, window_size)` array with -1 and scattering sparse calls at
`positions`. `/window_metadata`, `/selection_summary`, `/thresholds`, and
`/provenance_json` describe selection and provenance. TSV sidecars support
manual inspection.

## Posterior inference

Example neural and fixed-model commands:

```bash
PYTHON=/home/groups/btrippe/diamant/miniforge/envs/hmm/bin/python
INPUT=/scratch/users/diamant/mesc_validation/prepared/chr8-q50-valid90-reads10/dataset.h5
OUT=/scratch/users/diamant/mesc_validation/posteriors

"$PYTHON" -m smf_hmm.cli.infer_posterior \
  "$INPUT" "$OUT/neural.chr8.h5" \
  --checkpoint-dir /scratch/users/diamant/2026-09-08_fixed_phase_type_kernel_9/checkpoints \
  --batch-size 32

"$PYTHON" -m smf_hmm.cli.infer_posterior \
  "$INPUT" "$OUT/geometric.chr8.h5" --baseline geometric --batch-size 32

"$PYTHON" -m smf_hmm.cli.infer_posterior \
  "$INPUT" "$OUT/phase-type.chr8.h5" --baseline phase_type --batch-size 32
```

Inference computes prior terms once per window, then batches reads. The last
read batch is padded with missing dummy observations and trimmed before
writing. For each true read, it computes smoothed microstate marginals, sums
microstates into the model's named biological groups, and saves a probability
for every base, including unobserved positions. The prepared file retains
which positions had observations. The output is published only after success;
existing files are not overwritten and partial runs are not resumed.

The prediction file concatenates rows in `window_ids` order; within a window,
input read order is preserved. `window_offsets[i:i+2]` gives the row interval
for window `i`. `group_names` gives the final probability-axis order. The
prediction provenance records prepared input identity, model specification,
a SHA-256 digest of loaded model arrays, source, versions, and batch size.

## Validation completed

- `tests/test_tf_chip_prepare_dataset.py`: 13 synthetic tests passed, covering
  per-TF quantiles, strict ties, CpG/GpC context at window edges, unsorted TSV
  reads, full-span filtering, validity thresholds, HDF5 metadata/arrays,
  no-overwrite behavior, and failure without publishing partial output.
- Posterior, fixed-HMM, and prior tests: 11 passed for normalized grouped
  probabilities, missing-base predictions, fixed-size batch padding, IDs and
  offsets, empty datasets, and failed-run atomicity.
- The checkpoint at
  `/scratch/users/diamant/2026-09-08_fixed_phase_type_kernel_9/checkpoints`
  loaded on an H100 as a `SimpleCNNPriorArchitecture` with `K=138`. A synthetic
  33-read, 512-bp run produced finite `(33, 512, 3)` output whose largest
  element-sum error was `3.6e-7`. Runtime was about 3.2 seconds including JAX
  compilation; this is a smoke test, not a chromosome throughput estimate.
- The old CPU smoke-test attempt could not restore the Orbax checkpoint
  because it recorded `cuda:0`; rerunning on the available H100 succeeded.

# Related Notes

- [Fixed-HMM posterior-sampling baseline comparison](../2026-09-04/fixed-hmm-posterior-sampling-comparison.md): Records the original baseline models and motivates deterministic posterior summaries.
- [Multi-TF motif caches and state-call visualization plan](../2026-09-16/multi-tf-motif-cache-and-state-call-visualization.md): Documents the shared motif/ChIP caches reused for site selection.
- [CTCF motif and local-control ChIP–nexus signal collection](../2026-09-02/ctcf-motif-chip-nexus-local-controls.md): Describes motif and stranded ChIP cache conventions.

# Open Questions

- How many sites/windows/read pairs survive the default thresholds on chr8 and chr16, and how much time and disk space do preparation and inference use?
- Should the median ChIP filter remain the default after inspecting per-TF selection summaries? It selects a high-signal subset but is not a matched motif-absent control.
- Is 90% validity across assayable CpG/GpC sites and 10 reads per window useful across all four TFs?
- Can a preparation path based on current TSVs be replaced or supplemented with original BAM/SMF data to preserve molecule identity and alignment bounds across overlapping source regions?
- Should grouped posterior marginals and their outputs be compared against established numerical references at longer windows and larger batches?

# Sources

- `smf_hmm/evaluation/tf_chip/prepare_dataset.py` and `tests/test_tf_chip_prepare_dataset.py`.
- `smf_hmm/smf_hmm/posterior_inference.py`, `smf_hmm/cli/infer_posterior.py`, and `tests/test_posterior_inference.py`.
- `smf_hmm/smf_hmm/fixed_models.py`, `smf_hmm/smf_hmm/transition_modules.py`, and fixed-HMM migration tests.
- `smf_hmm/evaluation/tf_chip/configs/mesc_validation.yaml` and existing motif/ChIP cache shards under `/scratch/users/diamant/mesc_validation/benchmark_cache`.
- `smf_hmm/evaluation/datasets/smf/scripts/h5_to_calls.py` and `evaluation/datasets/smf/SPEC.md` (the latter is outdated regarding the assay and should be updated).
- `/scratch/users/diamant/2026-09-08_fixed_phase_type_kernel_9/checkpoints`.
