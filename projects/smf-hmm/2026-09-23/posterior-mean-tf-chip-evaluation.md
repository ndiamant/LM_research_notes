---
title: Posterior-mean TF/ChIP evaluation and example plots
date: 2026-09-23
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /home/users/diamant/repos/smf_hmm/smf_hmm/posterior_inference.py
  - /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/summarize_posteriors.py
  - /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/compare_posteriors.py
  - /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/plot_posterior_examples.py
  - /scratch/users/diamant/mesc_validation/posteriors
  - /scratch/users/diamant/mesc_validation/site_summaries
  - /scratch/users/diamant/mesc_validation/plots/posterior_examples.chr8.pdf
tags:
  - smf-hmm
  - tf-chip
  - posterior-marginals
  - evaluation
  - visualization
  - geometric-hmm
  - neural-prior
---

# Summary

Three scripts now evaluate posterior files from `smf_hmm.cli.infer_posterior`
(see [Motif-window posterior inference pipeline](motif-window-posterior-inference.md)):

- `summarize_posteriors.py` reduces each posterior file to a per-site table.
- `compare_posteriors.py` computes Spearman correlations with ChIP across runs.
- `plot_posterior_examples.py` draws one read plus each method's posterior
  curves at high-evidence sites.

A bug fix was needed first. Posterior inference failed on real data because
GPU TF32 matrix multiplication broke the group sum-to-one check.

The first chr8 comparison (neural checkpoint vs fixed geometric HMM) is mixed.
The neural model is better on CTCF. The geometric model is better on esrrb and
oct4, and it is much better at sites with no CpG/GpC inside the motif.

# Key Points

- **TF32 bug:** `infer_posterior` summed microstate marginals into groups
  with `probs @ membership`. On GPU, JAX's default matmul precision (TF32) gave
  group sums of 1 ± ~3e-4, which failed the `atol=1e-5` check on the first
  batch for both fixed baselines (the neural run evidently completed).
  It now uses `precision=jax.lax.Precision.HIGHEST`, giving errors of about
  5e-7. The per-position `log_softmax` in `HMM` was already correct.
- **Site score:** P(TF) averaged over motif bases, then over reads.
- **Coverage:** a read is "covered" if it has at least one observed position
  inside the motif (no flank). Occupancy is reported for all, covered, and
  uncovered reads.
- **Coverage is effectively a site property in this dataset.** Because
  preparation keeps reads with >=90% of assayable positions observed, almost
  every read at a site with an in-motif CpG/GpC is covered, and none are at
  sites without one. So "uncovered" means sites whose call comes from flanking
  observations plus the prior.
- **Metric:** Spearman against `chip_pooled_sum`, per TF and pooled over TFs.
  A site counts toward a subset only if it has >=10 reads in that subset.
  Pooling mixes TF-specific ChIP scales, so treat it as a rough summary.
- **Point estimates only:** there are no confidence intervals yet.

# Details

## Coverage diagnostic

On 3,000 sampled windows from the chr8 prepared dataset:

| TF | Sites with assayable position in motif | Reads observed in motif | ±10 bp | ±20 bp |
|---|---|---|---|---|
| ctcf | 90% | 89% | 97% | 99% |
| esrrb | 57% | 54% | 87% | 98% |
| klf4 | 63% | 61% | 89% | 96% |
| oct4 | 63% | 60% | 84% | 95% |

Widening the flank quickly makes nearly every read covered, which is why the
definition uses motif bases only.

## Commands

Environment: `hmm` conda env, run from `/home/users/diamant/repos/smf_hmm` on a
Sherlock compute node (job 44690718). Each summary took ~1.5 min and 275 MB RSS.
The example PDF took ~14 s.

```bash
POST=/scratch/users/diamant/mesc_validation/posteriors
SUMMARIES=/scratch/users/diamant/mesc_validation/site_summaries

python -m evaluation.tf_chip.summarize_posteriors \
    "$POST/neural.chr8.h5" "$SUMMARIES/neural.chr8.tsv.gz"
python -m evaluation.tf_chip.summarize_posteriors \
    "$POST/geometric.chr8.h5" "$SUMMARIES/geometric.chr8.tsv.gz"

python -m evaluation.tf_chip.compare_posteriors \
    --run neural="$SUMMARIES/neural.chr8.tsv.gz" \
    --run geometric="$SUMMARIES/geometric.chr8.tsv.gz" \
    --output "$SUMMARIES/metrics.chr8.tsv"

python -m evaluation.tf_chip.plot_posterior_examples \
    --sites "$SUMMARIES/neural.chr8.tsv.gz" \
    --run neural="$POST/neural.chr8.h5" \
    --run geometric="$POST/geometric.chr8.h5" \
    --output /scratch/users/diamant/mesc_validation/plots/posterior_examples.chr8.pdf
```

The posterior files came from
`prepared/chr8-q50-valid90-reads10/dataset.h5`, with the `prepare_dataset`
and `infer_posterior` commands recorded in the
[pipeline note](motif-window-posterior-inference.md). The geometric file was
produced after the TF32 fix. `summarize_posteriors` finds the prepared dataset
through posterior provenance. It errors if that file's size or mtime changed,
or if molecule IDs disagree. No script overwrites existing outputs.

## chr8 results (Spearman with pooled ChIP)

| TF | Subset | Sites | Geometric | Neural |
|---|---|---|---|---|
| ctcf | all | 7,442 | 0.496 | 0.520 |
| ctcf | covered | 6,640 | 0.520 | 0.549 |
| ctcf | uncovered | 1,009 | 0.374 | 0.372 |
| esrrb | all | 3,165 | 0.281 | 0.163 |
| esrrb | covered | 1,797 | 0.194 | 0.097 |
| esrrb | uncovered | 1,508 | 0.354 | 0.198 |
| klf4 | all | 13,189 | 0.109 | 0.109 |
| klf4 | covered | 8,443 | 0.052 | 0.082 |
| klf4 | uncovered | 5,267 | 0.380 | 0.180 |
| oct4 | all | 1,720 | 0.292 | 0.082 |
| oct4 | covered | 1,106 | 0.187 | 0.076 |
| oct4 | uncovered | 669 | 0.464 | 0.147 |
| pooled | all | 25,516 | 0.238 | 0.216 |
| pooled | covered | 17,986 | 0.236 | 0.266 |
| pooled | uncovered | 8,453 | 0.374 | 0.116 |

A site can pass the 10-read minimum in both covered and uncovered subsets.

Interpretation (unverified hypothesis): the geometric prior has no sequence
information, so at uncovered sites its motif P(TF) is driven by flanking
observations. Its high uncovered correlations (0.35–0.46) may reflect local
accessibility, which tracks ChIP within the above-median selection. If so, it
is an "accessibility only" reference that the neural model fails to reach for
esrrb, klf4, and oct4. Covered sites also correlate worse than uncovered ones
for those TFs in both models, which is unexpected if in-motif observations
carry the TF signal.

## Example plots

For each TF, sites are ranked by the minimum of their within-TF ChIP and
motif-score percentiles. Percentiles are computed before the coverage filter.
Only sites with >=10 covered reads are kept, overlapping windows are skipped,
and the top 5 per TF are plotted. Each page shows the covered read with the
most observed positions. Below it is one row per method with open/histone/TF
curves (`smf_hmm.plot.plot_state_prior_distribution`), and the motif is shaded
on every row.

- On the top CTCF site, neural gives a sharp motif-aligned TF call (mean
  P(TF) 0.77). Geometric splits the footprint between histone and TF (0.32).
- On the top oct4 site, the chosen read is protected across the motif, and
  both models call nucleosome. The read choice is method-independent, so it
  does not guarantee a TF-bound molecule.

## Validation

- `tests/test_tf_chip_summarize_posteriors.py`,
  `tests/test_tf_chip_compare_posteriors.py`, and
  `tests/test_tf_chip_plot_posterior_examples.py`: 8 synthetic tests pass.
  They cover exact occupancy values, changed-input and misaligned-run
  rejection, per-subset site counts, ranking and overlap, and page count.
- The 12 existing posterior tests pass after the TF32 fix.
- `black` is not installed in the `hmm` env, so the new files are unformatted.

# Related Notes

- [Motif-window posterior inference pipeline](motif-window-posterior-inference.md): Produces the prepared dataset and posterior files evaluated here.
- [Fixed-HMM posterior-sampling baseline comparison](../2026-09-04/fixed-hmm-posterior-sampling-comparison.md): Older hard-call/sampling comparison via `compare_runs.py`, which these scripts replace for posterior files.
- [Multi-TF motif caches and state-call visualization plan](../2026-09-16/multi-tf-motif-cache-and-state-call-visualization.md): Earlier hard-call example plots and the percentile site-confidence idea reused here.

# Open Questions

- Does local accessibility alone explain the geometric model's uncovered-site correlations? One check is to correlate the raw accessible fraction near each motif with ChIP.
- Why do covered sites correlate worse than uncovered ones for esrrb, klf4, and oct4?
- Do the rankings hold with block-bootstrap intervals, on chr16, and for the phase-type baseline?
- Should example reads be chosen by TF-likeness rather than observation count? That would bias the choice toward whichever method does the ranking.

# Sources

- `smf_hmm/smf_hmm/posterior_inference.py` (TF32 fix) and `tests/test_posterior_inference.py`.
- `smf_hmm/evaluation/tf_chip/summarize_posteriors.py`, `compare_posteriors.py`, `plot_posterior_examples.py`, and their tests.
- `smf_hmm/evaluation/tf_chip/compare_runs.py` (previous comparison design).
- `/scratch/users/diamant/mesc_validation/site_summaries/{neural,geometric}.chr8.tsv.gz` and `metrics.chr8.tsv`.
- `/scratch/users/diamant/mesc_validation/plots/posterior_examples.chr8.pdf`.
