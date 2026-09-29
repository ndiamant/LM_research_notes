---
title: Chr16 marginal ChIP covariation, measurement width, and window-selection limits
date: 2026-09-28
project: smf-hmm
agent: Codex
status: draft
sources:
  - User discussion and chr16 preparation output, 2026-09-28
  - smf_hmm/evaluation/tf_chip/select_covariation_windows.py
  - smf_hmm/evaluation/tf_chip/analyze_covariation.py
  - smf_hmm/evaluation/tf_chip/plot_posterior_examples.py
  - smf_hmm/evaluation/tf_chip/explore_chip_window_widths.py
  - /scratch/users/diamant/mesc_validation/exploratory_chip_window_widths_v2
tags: [tf-chip, covariation, chip-nexus, posterior-marginals, window-selection]
---

# Summary

The evaluation now uses a shared SMF-window HDF5 for individual TF–ChIP agreement
and motif-pair marginal covariation across neural and fixed HMMs. Today we walked
through the chr16 workflow, investigated ChIP-window overlap, moved toward 51-bp
ChIP measurements, updated example plots, and added a 128-bp motif flank requirement.
The user ran the new preparation: only **1,814** motif-qualified windows remain
from **2,488** coverage-qualified central windows. Coverage already prevents
reaching 5,000 under the current candidate-generation scheme.

The next design discussion is **motif-centered candidates with read collection
independent of predefined SMF regions**. This is not implemented. Current code
still considers exactly one central 512-bp window per stored source region.

This supersedes the single-target-motif evaluation/plotting path documented in
[the September 23 evaluation note](../2026-09-23/posterior-mean-tf-chip-evaluation.md),
while retaining the shared posterior-inference format from
[the inference pipeline note](../2026-09-23/motif-window-posterior-inference.md).

# Key Points

- Target validation chromosome: chr16. Models: old neural run
  `2026-09-08_fixed_phase_type_kernel_9/checkpoints`, updated neural run
  `2026-09-26_fixed_phase_type_kernel_9/checkpoints`, fixed geometric and phase-type
  HMMs. The updated checkpoint subdirectory was assumed in recommended commands.
- User changed YAML `chip.half_window` to 25 and reported rebuilding chr16 caches:
  this means **51 bp**, with half-open bounds `[center-25, center+26)`.
  Existing 101-bp caches are not changed by editing YAML. Downstream interval
  construction uses cache metadata, not the current YAML.
- Analyze **marginals across sites/pairs**, not same-molecule co-binding: no
  ground truth exists for the latter. The model's generic TF state gets its
  TF label from the motif annotation.
- Pair correlations now exclude overlapping motif spans and overlapping ChIP
  intervals. Bins extending beyond the inference-window length and empty bins
  are omitted. Touching half-open intervals are permitted.
- New `--min-motif-flank 128` filters motifs and windows before sampling. With
  512-bp inference the full motif must fit inside the central 256 bp. ChIP
  strength never enters selection. Pair endpoints must both qualify.
- Preparation still has no minimum assayable-site count beyond one. A qualifying
  read has calls at >=90% of reference CpG/GpC cytosines; >=10 such reads are
  required. GCG's central cytosine is counted once.

# Details

## Cache protocol clarified

`benchmark_cache.py` builds one TF/chromosome shard using Biopython `Bio.motifs`.
It loads a MEME motif, applies pseudocount 0.5, derives a PSSM and a score threshold
with nominal background exceedance probability 1e-4, and scans both strands.
This is not an FDR or posterior binding probability. Scores, not per-hit p-values,
are stored.

The order-0 background is A/C/G/T frequencies over the entire supplied FASTA,
excluding ambiguous bases from the denominator but not excluding blacklisted
regions. It is shared via `genome_background.json`. The configured
`mm10-blacklist.v2.bed.gz` excludes sites whose entire ChIP measurement interval
intersects a listed interval; local-file identity against the official download
was not audited in this discussion.

Each shard contains `sites_chip.tsv.gz` and `metadata.json`. Rows hold motif
identity, coordinates, strand, score, and per-replicate positive/absolute-negative
ChIP sums. Metadata records thresholds, background, inputs, versions and checksum.
The cache holds interval sums, not base-resolution profiles. Narrowing a cached
measurement requires rereading the original BigWigs.

## ChIP-window sensitivity report

At the user's request a GPT-6-sol subagent compared 31/51/101-bp windows on fixed
chr8 sites/pairs. Script: `evaluation/tf_chip/explore_chip_window_widths.py`.
Report: `$SCRATCH/mesc_validation/exploratory_chip_window_widths_v2/chip_window_widths.pdf`.
This was run in compute allocation 45539399 (8 CPUs, 125 GB reported at the time).
All **81,294 site-replicate 101-bp sums matched the original cache exactly**.
The six-page PDF was rendered and visually checked.

For identical non-overlapping-motif pairs at 50–100 bp separation:

| TF | Spearman, 101 bp | Spearman, 31 bp |
|---|---:|---:|
| CTCF | .935 | .898 |
| ESRRB | .881 | .753 |
| KLF4 | .890 | .780 |
| OCT4 | .952 | .850 |

All 101-bp intervals overlap in that bin. However, correlations also fall when
narrowing windows on the common subset with disjoint 101-bp intervals. Width
sensitivity cannot be attributed entirely to direct shared-count reuse: signal
capture and spatially correlated binding signal also change. OCT4 has only 22
unique pairs in the 50–100 bin, from 21 windows.

Median per-site fraction of the surrounding 201-bp signal captured:

| TF | 31 bp | 51 bp | 101 bp |
|---|---:|---:|---:|
| CTCF | .188 | .312 | .578 |
| ESRRB | .179 | .289 | .550 |
| KLF4 | .189 | .296 | .545 |
| OCT4 | .189 | .305 | .565 |

The user observed spikes near motif edges for CTCF/OCT4 and more internal peaks
for KLF4/ESRRB. This was a visual interpretation, not a fitted footprint boundary;
the plots did not mark motif boundaries. User selected 51 bp as a practical
working width. No width is established as universally optimal.

## Analysis and example plotting changes

`analyze_covariation.py` averages TF posterior over reads and motif bases,
also saving the sum over motif bases (expected TF-covered bases after read
averaging). Its first two PDF pages show per-TF/model rank scatter plots and
linear occupancy versus log1p pooled ChIP, respectively. Subsequent matrices
compare ChIP and model Spearman correlations on identical pairs. Same-TF pairs
include both orientations; sample counts remain unique unordered genomic pairs.
No independence-based significance tests are reported.

Repeated genomic motifs/pairs use the window maximizing edge margin, then window
ID; a pair always uses both sites from one window. Disjoint ChIP intervals are
now the only pair-analysis scope. The `all` aggregate still includes all eligible
pairs, even if a distance lies outside displayed complete bins. The saved
`matched_pairs.tsv.gz` remains unfiltered for inspection.

`plot_posterior_examples.py` now takes the prepared covariation HDF5 directly,
with multiple motifs per window; old `summarize_posteriors.py` output is no longer
required. It selects up to five non-overlapping example windows per TF using
strong ChIP/motif percentile evidence, requiring >=10 reads observed inside the
target motif. Each model selects **one highest-mean-TF-posterior motif-covered
read**. Shared winners are deduplicated; ties use input order. Each winning read
gets a page with observed calls and every model's posterior tracks, shaded target
motif and selecting-model labels. These are high-scoring examples, not
representative population samples. It now runs without JAX in `smf_clean`.

## 128-bp flank result and limitation

Motif-edge examples motivated requiring enough context to reveal much of a
nucleosome rather than confusing a truncated long footprint with a TF. The
128-bp margin is a conservative working choice, not a guarantee against ambiguity.
`select_covariation_windows.py` annotates all coverage-qualified candidates,
filters motif spans by margin, removes motif-free candidates, then samples.
`candidate_windows.tsv.gz` and provenance preserve the pre-sampling counts.
The default margin is zero; even then the updated workflow requires >=1 motif.

User-reported chr16 result (755 MB TSV; second scan ~1m56s):

- 2,488 deduplicated coverage-qualified central windows.
- 1,814 with >=1 motif satisfying the 128-bp flank requirement (73%).
- 940 with >=2 qualifying motifs; these counts precede disjoint-ChIP filtering.
- All 1,814 selected because fewer than the requested 5,000 were available.

| Pair | Windows | Unique pairs | Motif-disjoint pairs | ChIP-overlapping pairs |
|---|---:|---:|---:|---:|
| CTCF–CTCF | 150 | 217 | 153 | 124 |
| CTCF–ESRRB | 93 | 137 | 127 | 60 |
| CTCF–KLF4 | 290 | 683 | 594 | 359 |
| CTCF–OCT4 | 55 | 67 | 66 | 28 |
| ESRRB–ESRRB | 45 | 59 | 44 | 45 |
| ESRRB–KLF4 | 191 | 328 | 319 | 208 |
| ESRRB–OCT4 | 42 | 50 | 49 | 30 |
| KLF4–KLF4 | 420 | 1275 | 1077 | 974 |
| KLF4–OCT4 | 111 | 180 | 177 | 103 |
| OCT4–OCT4 | 37 | 41 | 36 | 32 |

Do not subtract both overlap columns independently: their excluded sets overlap.
OCT4–OCT4 has at most nine ChIP-disjoint pairs before distance binning. Pair
availability is much sparser than window count alone suggests.

## Proposed next design — not implemented

Generate 512-bp candidates around all cached motif centers, then recheck SMF
coverage at those positions. This could recover motifs near the edge of the
current central-source windows. Pair-midpoint candidates could later recover
longer separations: a motif-centered window plus 128-bp margins only admits
neighboring motif centers roughly within 128 bp, whereas midpoint centering can
fit pairs into the 256-bp interior.

The user objected to assuming future reads are sorted or partitioned into source
regions. Agreed distinction:

- Candidate generation should be independent of SMF storage format.
- A general collector can interval-index candidates and stream reads in arbitrary
  order, evaluating each spanning read's observations within matching windows.
  A second pass can collect reads for the selected candidates. No sorting is
  inherently required.
- Today's TSVs were exported by `h5_to_calls.py`: bounds describe stored regions
  and IDs are synthetic region/row IDs. Physical molecules duplicated across
  overlapping source regions cannot be reliably identified. Selecting one source
  region per candidate remains a data-specific safeguard, not a general API rule.
- Future input with actual read bounds and stable molecule identities could use
  the general collector. An interval index alone cannot recover lost identity.

## Commands and validation provenance

Run from `$HOME/repos/smf_hmm`; CPU preparation/analysis in `smf_clean`, GPU
inference in `hmm`. These are the recommended chr16 reproduction commands;
only the flank128 preparation counts above were explicitly pasted by the user.
Do not assume the new flank128 dataset has already been inferred.

```bash
conda activate smf_clean
RUN_ROOT="$SCRATCH/mesc_validation"
CACHE_ROOT="$RUN_ROOT/benchmark_cache"  # replace if 51-bp caches use another root
for TF in ctcf oct4 klf4 esrrb; do
  python -m evaluation.tf_chip.benchmark_cache \
    --config evaluation/tf_chip/configs/mesc_validation.yaml \
    --cache-root "$CACHE_ROOT" --chrom chr16 --tf "$TF" || break
done
# Cache commands refuse to replace existing shards.
DATA_DIR="$RUN_ROOT/covariation_windows/chr16-512-valid90-reads10-flank128"
python -m evaluation.tf_chip.select_covariation_windows \
  --config evaluation/tf_chip/configs/mesc_validation.yaml \
  --cache-root "$CACHE_ROOT" --chrom chr16 --window-size 512 \
  --min-motif-flank 128 --min-valid-fraction 0.9 --min-reads 10 \
  --max-windows 5000 --seed 123 --output-dir "$DATA_DIR"

conda activate hmm
mkdir -p "$DATA_DIR/posteriors"
python -m smf_hmm.cli.infer_posterior "$DATA_DIR/dataset.h5" \
  "$DATA_DIR/posteriors/neural_old.h5" \
  --checkpoint-dir "$SCRATCH/2026-09-08_fixed_phase_type_kernel_9/checkpoints" --batch-size 32
python -m smf_hmm.cli.infer_posterior "$DATA_DIR/dataset.h5" \
  "$DATA_DIR/posteriors/neural_updated.h5" \
  --checkpoint-dir "$SCRATCH/2026-09-26_fixed_phase_type_kernel_9/checkpoints" --batch-size 32
for BASELINE in geometric phase_type; do
  python -m smf_hmm.cli.infer_posterior "$DATA_DIR/dataset.h5" \
    "$DATA_DIR/posteriors/$BASELINE.h5" --baseline "$BASELINE" --batch-size 32
done

conda activate smf_clean
python -m evaluation.tf_chip.analyze_covariation "$DATA_DIR/dataset.h5" \
  --prediction neural_old="$DATA_DIR/posteriors/neural_old.h5" \
  --prediction neural_updated="$DATA_DIR/posteriors/neural_updated.h5" \
  --prediction geometric="$DATA_DIR/posteriors/geometric.h5" \
  --prediction phase_type="$DATA_DIR/posteriors/phase_type.h5" \
  --output-dir "$DATA_DIR/analysis"
python -m evaluation.tf_chip.plot_posterior_examples "$DATA_DIR/dataset.h5" \
  --run neural_old="$DATA_DIR/posteriors/neural_old.h5" \
  --run neural_updated="$DATA_DIR/posteriors/neural_updated.h5" \
  --run geometric="$DATA_DIR/posteriors/geometric.h5" \
  --run phase_type="$DATA_DIR/posteriors/phase_type.h5" \
  --sites-per-tf 5 --output "$DATA_DIR/posterior_examples.pdf"
```

Changing dataset selection requires new inference; changing only plot settings
or correlation filtering does not. Use new output paths on reruns.
The older no-flank chr16 dataset is under `chr16-512-valid90-reads10`.

Width-report reproduction (use a new output directory):

```bash
python -m evaluation.tf_chip.explore_chip_window_widths \
  --input "$SCRATCH/mesc_validation/covariation_windows/chr8-512-valid90-reads10" \
  --config evaluation/tf_chip/configs/mesc_validation.yaml \
  --output "$SCRATCH/mesc_validation/exploratory_chip_window_widths_rerun"
```

Actual synthetic validation commands used the explicit `smf_clean` interpreter:

```bash
MPLCONFIGDIR=/Users/ndiamant/repos/smf_hmm/.mplconfig \
  /home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python -m pytest -q \
  tests/test_tf_chip_plot_posterior_examples.py tests/test_tf_chip_analyze_covariation.py
# 10 passed
MPLCONFIGDIR=/Users/ndiamant/repos/smf_hmm/.mplconfig \
  /home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python -m pytest -q \
  tests/test_tf_chip_select_covariation_windows.py
# 5 passed, including exact flank boundaries and filtering before sampling
```

The repo-requested macOS-style Matplotlib cache path may fall back to a temporary
cache on Sherlock. Black was unavailable in checked environments; no installation
was performed. Code changes have not been committed by this note task; the code
repo contains unrelated and user-staged work. Only the notes are committed here.

# Related Notes

- [Posterior-mean evaluation and examples](../2026-09-23/posterior-mean-tf-chip-evaluation.md): Previous single-target summary/plot interfaces, now superseded for covariation examples.
- [Motif-window posterior inference](../2026-09-23/motif-window-posterior-inference.md): Shared HDF5 and inference machinery retained here; its ChIP-thresholded selection is distinct from this workflow.
- [Multi-TF motif caches](../2026-09-16/multi-tf-motif-cache-and-state-call-visualization.md): Background on reusable per-TF/chromosome motif/ChIP caches.

# Open Questions

- How much does motif-centered candidate generation recover at the same 128-bp
  margin and coverage thresholds? Count before changing inference workloads.
- What should the general read-collection interface require for molecule identity
  and genomic coordinates, while accommodating current lossy region-based TSVs?
- How should overlapping candidate windows be sampled/deduplicated so increased
  window counts do not masquerade as independent regions?
- Are 128-bp flanks sufficient/necessary for avoiding truncated nucleosome calls?
- Should an explicit minimum CpG/GpC site count be added based on diagnostics?
- Are rare TF-pair counts adequate after ChIP-disjoint filtering and distance bins?

# Sources

- User discussion, pasted chr16 preparation log/table, and code inspected/edited
  in `$HOME/repos/smf_hmm/evaluation/tf_chip/`: `benchmark.py`, `benchmark_cache.py`,
  `configs/mesc_validation.yaml`, `select_covariation_windows.py`,
  `prepare_dataset.py`, `analyze_covariation.py`, `plot_posterior_examples.py`,
  `explore_chip_window_widths.py`, and `workflow/README.md`.
- `$SCRATCH/mesc_validation/exploratory_chip_window_widths_v2/`:
  `chip_window_widths.pdf`, `same_tf_correlations.tsv`, `width_capture.tsv`,
  `cache_audit.tsv`, `site_signals.tsv.gz`, `oriented_profiles.npz`, `provenance.json`.
- [Original ChIP-nexus paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4390430/): Signal reflects exonuclease-stop patterns; motif span alone need not capture it.
- [ENCODE blacklist paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6597582/) and
  [official mm10 blacklist](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz): Background on problematic-region exclusion.
