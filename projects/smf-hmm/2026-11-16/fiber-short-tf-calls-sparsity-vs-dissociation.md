---
title: Short TF calls at bound motifs are neither mainly sparse-site artifacts nor detectable partial footprints — except NRF1 (sparse GC-rich motif) and CTCF (real footprints, rare partial ones), bound motifs of the other seven TFs look like accessible background in the reads
date: 2026-11-16
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /scratch/users/diamant/fiber_eval_work/tf_mechanism/short_footprints.py -> short_footprints/{calls.tsv.gz, modelfree.tsv.gz, summary.txt}
  - Viterbi segments of the GELU KL 1000 model on all val windows (viterbi_lengths.py, MODELS=fix_kl1000 STRIDE=1; Slurm 52570388) -> viterbi/viterbi_segments_kl1000_x1_stride1.tsv.gz
  - IDR-bound / unbound benchmark motif instances (analysis_idr_gelu/unique_motif_summaries.tsv.gz), val chr20 + chr8
  - User question, 2026-11-16 (short calls at SPI1/USF2/RFX5/EBF1: partial footprints from dissociation during labelling, or something else?)
tags: [smf-hmm, fiberseq, tf-footprints, tf-duration, dissociation, callable-sites]
---

# Summary

**Question.** With the GELU KL 1000 model, Viterbi TF calls at IDR-bound SPI1, USF2, RFX5 and EBF1
motifs have median lengths of 11-14 bp. Before raising the TF minimum to about 15 bp: are these short
calls
- partial footprints from TFs dissociating during the minutes-long labelling, or
- under-determined lengths at motifs with few callable A/T sites?

**Four tests** (all val windows):
1. Callable-site density.
2. Contradicting m6A in or next to the call.
3. Call position relative to the motif.
4. Model-free protected runs, bound vs unbound. Restricted to accessible reads (protected run < 80 bp)
   to remove the nucleosome confound of unbound motifs.

**Results.**
- **Sparsity is not the general explanation.** Short calls are not concentrated at low-callable motifs;
  SPI1, the most short-heavy (62%), has 22 callable A/T in motif ±10 bp. The exception is **NRF1**:
  - its GC-rich motif leaves about 7 callable sites;
  - its short calls cover the motif centre (89%) with no m6A inside the motif (8%);
  - 97% could extend to 15 bp without an m6A.

  These are real footprints with under-determined length.
- **CTCF and NRF1 show real footprints in the reads.** In accessible reads their bound motifs have long
  protected runs (median 18 / 16 bp between outermost protected calls, vs 6 / 11 unbound) and fewer
  mixed motif spans (0.53 / 0.32 vs 0.77 / 0.46). CTCF's short calls are only 3% of its calls; they are
  off-centre with m6A inside the motif, consistent with partial footprints (dissociation) but rare.
- **For SPI1, USF2, PAX5, EBF1, RFX5, CREM and YY1, bound motifs look like accessible background.**
  - In accessible reads their motif spans are about 80% mixed, with 3-6 bp protected runs, the same as
    at unbound motifs.
  - Their short calls usually have an m6A inside the motif (58-93%, vs 8-55% for long calls) and often
    miss the motif centre (21-61% cover it, vs 52-95%).
  - They match the background protected patches of the
    [open-overdispersion note](fiber-open-state-overdispersion.md), turned into TF by the prior's motif
    bump, rather than partial footprints distinguishable from background.

**Interpretation.**
- Dissociation during labelling may be the biological reason these footprints are weak or absent (short
  residence time), but the reads cannot separate such footprints from background patches. So these short
  calls carry little footprint information either way.
- A TF minimum of about 15 bp should cost little real signal: NRF1's short calls would mostly extend, and
  CTCF's are rare.
- The benchmark AUROC of the weak-footprint TFs may drop where it rode on the prior's motif bump
  (inferred).

# Key Points

- **Model calls at IDR-bound motifs** (GELU KL 1000, Viterbi, all val windows; short = ≤ 14 bp):

  | TF | calls | short | callable A/T, motif ±10 (short / long) | m6A in motif (short / long) | short extendable to 15 bp | covers motif centre (short / long) |
  |---|---|---|---|---|---|---|
  | CTCF | 17,295 | 3% | 14 / 15 | 74% / 11% | 80% | 41% / 97% |
  | NRF1 | 3,061 | 9% | 7 / 7 | 8% / 6% | 97% | 89% / 96% |
  | YY1 | 683 | 11% | 15 / 15 | 74% / 8% | 64% | 21% / 95% |
  | CREM | 572 | 28% | 11 / 9 | 64% / 31% | 84% | 44% / 77% |
  | USF2 | 533 | 30% | 13 / 12 | 58% / 19% | 68% | 61% / 93% |
  | EBF1 | 1,245 | 31% | 15 / 12 | 87% / 15% | 72% | 40% / 91% |
  | RFX5 | 221 | 35% | 13 / 13 | 45% / 35% | 73% | 82% / 92% |
  | PAX5 | 624 | 39% | 16 / 14 | 92% / 55% | 75% | 33% / 52% |
  | SPI1 | 1,801 | 62% | 22 / 22 | 93% / 69% | 60% | 61% / 70% |

- **Model-free, accessible reads** (a protected call in the motif span; protected run < 80 bp):

  | TF | mixed span, bound / unbound | run between outermost protected calls, bound / unbound (median bp) | run between flanking m6A, bound / unbound (median bp) |
  |---|---|---|---|
  | CTCF | 0.53 / 0.77 | 18 / 6 | 23 / 12 |
  | NRF1 | 0.32 / 0.46 | 16 / 11 | 26 / 19 |
  | SPI1 | 0.91 / 0.85 | 3 / 6 | 5 / 9 |
  | USF2 | 0.83 / 0.83 | 3 / 4 | 6 / 7 |
  | PAX5 | 0.86 / 0.82 | 3 / 6 | 9 / 12 |
  | EBF1 | 0.79 / 0.75 | 6 / 8 | 11 / 14 |
  | RFX5 | 0.79 / 0.80 | 5 / 6 | 9 / 11 |
  | CREM | 0.81 / 0.77 | 3 / 5 | 6 / 7 |
  | YY1 | 0.78 / 0.77 | 5 / 7 | 7 / 9 |

# Details

- **Extendable** = some ≥ 15 bp window containing the call has no m6A call in the read.
- **Model-free runs** start from the motif's protected call nearest its centre and extend over
  consecutive protected calls.
- **Unrestricted test 4 is confounded.** Bound vs unbound mixed shares without the < 80 bp restriction
  favour "bound is more mixed" only because unbound motifs mostly sit in nucleosomes (protected runs
  about 105 bp); the restricted version above removes that.

```bash
T=/scratch/users/diamant/fiber_eval_work/tf_mechanism
MODELS=fix_kl1000 STRIDE=1 sbatch $T/viterbi_lengths.sbatch     # 52570388
python $T/short_footprints.py
```

# Related Notes

- [Open-state overdispersion](fiber-open-state-overdispersion.md): The background protected patches.
- [Spacer and duration shaping](fiber-tf-nubs-spacer-and-duration-shaping.md): TF length histograms per
  KL setting.
- [Specificity metrics and debiased kernel](../2026-11-13/fiber-specificity-metrics-debiased-kernel.md):
  Which TFs footprint in this data.

# Open Questions

- Does a TF minimum of 15 bp (with tied TF potentials) keep CTCF and NRF1 sensitivity and remove the
  background nubs? Planned trial.
- Can dissociation be tested directly, e.g. footprint occupancy vs known residence times across more TFs?

# Sources

- Scripts and outputs in the frontmatter.
