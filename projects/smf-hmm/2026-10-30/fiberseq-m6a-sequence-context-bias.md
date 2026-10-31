---
title: m6A sequence-context bias in GM12878 Fiber-seq — how much context lookup tables need
date: 2026-10-30
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - smf_hmm repo, evaluation/datasets/smf/context_bias/ (README.md, aggregate_sites.py, fit_context_tables.py, plot_context_bias.py, results.json, adj_table_3_3.npz, context_bias.png), branch nd/experiments, uncommitted at time of writing
  - /scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/{active.h5,background.h5,manifest.tsv}
  - /scratch/groups/btrippe/brian/projects/SMF/fiberseq/reference/hg38.fa (+ .fai)
  - /scratch/groups/btrippe/ndiamant/hg38-blacklist.v2.bed.gz
  - /scratch/users/diamant/context_bias/ (sites.npz, fit.log)
  - User discussion, 2026-10-30
tags: [smf-hmm, fiberseq, m6a, gm12878, enzyme-bias, sequence-context, hia5, lookup-table, held-out]
---

# Summary

The user expected strong enzyme bias within about ±5 bp of each adenine. We measured how much
context matters by fitting lookup tables `P(m6A | context)` at increasing widths and scoring
held-out log-likelihood. The data is every called A/T site in the GM12878 active and background
windows: 162M sites and 19.7B calls. We fit two models, a raw table and one adjusted for local
accessibility.

**Most of the bias is within ±3 bp, and essentially all of it within ±4–5 bp.**

- ±1 alone captures about 75% of the best gain; ±3 gets 95%, which is FiberHMM's 3+3
  context.
- The effect is very large. At matched local accessibility, rates run from 0.002 for
  `AA[A]AA` to 0.24 for `CC[A]GC`. In the ±3 table, the 95th-percentile context has 64× the
  odds of m6A of the 5th.
- Sequence context explains more held-out likelihood than a 7–50 bp local accessibility
  estimate does: 0.048 vs 0.026 bits per call.

# Key Points

- **Held-out gain** over a model with no context, in bits per call, on test chromosomes chr8 and
  chr9:

  | context | ±1 | ±2 | ±3 | ±4 | ±5 | ±6 |
  | --- | --- | --- | --- | --- | --- | --- |
  | raw | 0.047 | 0.057 | 0.060 | 0.063 | 0.064 | 0.064 |
  | adjusted | 0.036 | 0.044 | 0.046 | 0.048 | 0.048 | 0.048 |

  The context-free cost is 0.528 bits (raw) and 0.502 bits (adjusted). The adjusted ±6 is
  slightly below ±5, so it starts to overfit.
- **The 3′ side carries more.** One-sided adjusted tables reach 0.029 bits (3′) and 0.022
  bits (5′) at 6 bases, against 0.048 for ±6, so both sides are needed. In the raw grid,
  (kL, kR) = (3, 4) gets 96% and (2, 3) 91%.
- **One A-centred table serves both strands.** Separate 3+3 tables for reference-A sites and
  for reference-T sites (reverse-complemented) correlate at 0.9998, with a mean absolute
  difference of 0.001. Splitting them gains nothing held-out. This extends the chr21 check in
  the data-preparation note genome-wide.
- **The bias is the same in open and closed chromatin.** Active-window and background-window
  tables correlate at 0.996, and splitting them gains 0.0002 bits per call. This points to the
  enzyme and the caller, not chromatin (inferred).
- **What drives it:** A/T runs are depleted (the lowest ±2 contexts are all A/T, `AA[A]AA`
  0.002) and G/C flanks are enriched (`CC[A]GC`, `AG[A]GC` about 0.23). These rates are at the
  median accessibility bin of active windows.

# Details

## Aggregation (`aggregate_sites.py`)

Each called site is collapsed over molecules to (m6A count, call count). Haplotype is ignored,
per the user. The site also gets:

- the A-centred ±6 bp hg38 context, reverse-complemented at T sites so it reads 5′→3′ along
  the methylated strand, and packed as 12 two-bit bases;
- the column m6A and call counts 7–50 bp around the site.

Other details:

- **Coordinates:** h5 position p of window `chrom:start-end` sits at 0-based genomic
  start − 1 + p. Checked: offset 0 puts 100% of positions on A/T, against 57–62% for ±1 and ±2.
- **Overlapping windows:** they share molecules, so any site an earlier window already
  covered is skipped.
- **Blacklist:** zero kept windows overlap it, confirmed by brute force, so the window builder
  had already excluded it. The filter is kept as a guard.
- **Output:** `$SCRATCH/context_bias/sites.npz`, 3.7 GB, built in about 2 min on 8 cores.

## Models (`fit_context_tables.py`)

The split is test chr8 + chr9 (14.0M sites), validation chr10 (8.0M) and train on the rest
(140M).

- **Raw:** a table over kL 5′ and kR 3′ bases, each entry shrunk toward its one-base-shorter
  parent as (m + α·p_parent) / (n + α). Calls at a site are correlated through accessibility,
  so the effective sample size is sites, not calls, and long contexts need the shrinkage. The
  validation set chose α = 30000 calls, the same as 3000 to four digits. The full 7×7 kL × kR
  grid is in `results.json`.
- **Adjusted:** `logit P = f[context] + g[bin]`. The bin is one of 50 quantiles of the 7–50 bp
  column rate, plus a no-data bin, crossed with window class. It is fit by block-Newton on rows
  aggregated by (context, bin), with `f` penalised toward its parent's `f`. Validation chose
  λ = 30, and λ from 3 to 300 gave identical scores at ±4.
- **Splits:** strand (A vs T) and class (active vs background) were tested as an extra key
  digit at ±3, shrunk toward the pooled table.
- **Runtime:** about 65 min inside a 125 GB allocation.

## Caveats

- This is the combined bias of Hia5, the PacBio kinetics-based m6A caller and the `ML >= 230`
  threshold these calls were built at. Lookup tables cannot separate them, and the A/T-run
  depletion may be partly the caller (speculative). A naked-DNA or no-Hia5 control, such as
  Drosophila S2 GSE249742 from the data-preparation note, would be needed.
- The true enzyme effect probably lies between the raw and adjusted estimates. The raw table
  overstates it through chromatin composition. The local accessibility estimate is itself
  depressed in A-rich DNA, because the bias acts on neighbouring sites too, so the adjusted
  `f` probably understates it. Structure finer than 50 bp, such as TF footprints, is not
  adjusted for.
- The FiberHMM checkout is not on this filesystem, so its `hia5_pacbio` 3+3 table was not
  compared.

## Commands

Environment: `PY=/home/groups/btrippe/diamant/miniforge/envs/hmm/bin/python`, run from
`evaluation/datasets/smf/context_bias/` on a Sherlock allocation (sh04-03n07, Slurm job
50179010, 8 CPUs, 125 GB).

```bash
$PY aggregate_sites.py                                   # -> $SCRATCH/context_bias/sites.npz
nohup $PY -u fit_context_tables.py > $SCRATCH/context_bias/fit.log 2>&1 &   # -> results.json, adj_table_3_3.npz
MPLCONFIGDIR=/home/users/diamant/repos/smf_hmm/.mplconfig $PY plot_context_bias.py   # -> context_bias.png
```

A separate check computed raw 3+3 rates for A and T sites straight from `sites.npz` with
`np.bincount` over the middle 6 flank bases. It was not saved as a script, and its result is the
0.9998 correlation above.

# Related Notes

- [Fiber-seq GM12878 data preparation for SMFNet](../../smfnet/2026-10-05/fiberseq-gm12878-data-preparation.md):
  Source of the h5 windows. Its chr21 A-vs-T(rc) check (±3 corr 0.990) is confirmed here
  genome-wide (0.9998), and it names Drosophila S2 as the control for Hia5 bias.
- [Posterior predictive checks of the fixed-blocks SMF HMM](../2026-10-22/posterior-predictive-checks.md):
  Warns that constant emissions let enzyme bias leak into TF states and plans per-k-mer emission
  offsets for Fiber-seq. This note sizes that effect and suggests ±3–4 bp is enough.
- [Flank baseline for profile metrics, logistic regression instead of a k-mer table](../../smfnet/2026-10-05/flank-baseline-logistic-regression.md):
  Chose an additive flank model over a backed-off lookup table for GCH, which is strongly
  strand-asymmetric. Fiber-seq is strand-symmetric, and this shows a backed-off table is
  tractable to ±5 with 162M sites.
- [Cross-read kernel denoising of undermethylated Fiber-seq](fiberseq-cross-read-kernel-denoising.md):
  Same-day experiment whose likelihood kernel assumed a flat P(m6A|acc). The 100-fold
  context spread found here is a likely reason to revisit it.

# Open Questions

- How much of the bias is Hia5 and how much the PacBio caller or the ML threshold? This would
  need naked-DNA or no-Hia5 data, or the same molecules recalled at ML ≥ 125.
- Does the adjusted 3+3 table (`adj_table_3_3.npz`) as per-site emissions improve the HMM's
  held-out site predictions, and does it reduce spurious short TF footprints in A/T-rich DNA?
- Is the bias separable into enzyme accessibility × caller sensitivity? Could a low-rank or
  additive model in place of the full table reach the same held-out likelihood?

# Sources

- smf_hmm `evaluation/datasets/smf/context_bias/` on branch `nd/experiments`, uncommitted:
  `README.md`, `aggregate_sites.py`, `fit_context_tables.py`, `plot_context_bias.py`,
  `results.json`, `adj_table_3_3.npz`, `context_bias.png`.
- `/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/` (`active.h5`,
  `background.h5`, `manifest.tsv`).
- `/scratch/groups/btrippe/brian/projects/SMF/fiberseq/reference/hg38.fa`;
  `/scratch/groups/btrippe/ndiamant/hg38-blacklist.v2.bed.gz`.
- `/scratch/users/diamant/context_bias/sites.npz` and `fit.log`.
- User discussion, 2026-10-30.
