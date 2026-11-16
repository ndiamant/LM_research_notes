---
title: Accessible DNA has clustered protection that one open state cannot represent — runs of consecutive protected calls are up to 2.4× more frequent than the trained open emission predicts, outside A/T tracts too; the 8 bp TF nubs are the cheapest patch, and the 40 bp spike without the prior is duration-cap censoring
date: 2026-11-16
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /scratch/users/diamant/fiber_eval_work/tf_mechanism/{open_overdispersion.py -> open_overdispersion.tsv; viterbi_lengths.py, viterbi_lengths.sbatch (Slurm 52534465/52534468), plot_viterbi_lengths.py -> viterbi/, figures/viterbi_lengths*.png}
  - Model gelu_s0 (GELU-head ablation B seed 0): learned context offsets and open protection 0.45
  - User discussion, 2026-11-16 ("something weird is happening for the prior model to put so many TFs at 8 bp; the prior-free model is also clearly wrong with its spike at 40 bp")
tags: [smf-hmm, fiberseq, open-state, overdispersion, tf-duration, model-design]
---

# Summary

**Observation** (Viterbi segment lengths, 64.7k val reads). TF lengths are wrong in both settings:
- **With the AG prior:** 54% exactly 8 bp, the minimum.
- **Without it:** a broad 13-35 bp hump plus 48% at exactly 40 bp, the maximum; other_long also piles
  at its 100 bp maximum.

**1. The 40 bp spike is censoring by hard duration caps** (inferred from the structure).
- TF is ≤ 40 bp, other_long is 40-100 bp, histone is exactly 147 bp.
- Protected stretches just above a cap get the cheapest covering state truncated at its maximum:
  41-60 bp stretches (e.g. two TFs closer than the 21 bp no-chaining minimum) → TF 40 + open;
  100-146 bp stretches → other_long 100.
- Model-free protected-run lengths are smooth across these boundaries (2026-11-13 nucleosome note).

**2. The 8 bp spike patches overdispersion of the open state.**
- **The test.** Within DNA the prior-free model calls open (2.0M called sites, interior open segments ≥
  40 bp, trimmed 10 bp per end), compare runs of consecutive protected calls with independent draws from
  the trained open emission, σ(logit 0.45 + context offset).
- **The mean fits** (0.423 observed vs 0.428), **but protected calls cluster:** runs of ≥ 4 are 1.13×
  more frequent than predicted, ≥ 6 1.47×, ≥ 8 1.84×, ≥ 12 2.37×. This is a lower bound, because Viterbi
  has already moved the strongest runs out of "open".
- **By labelling.** Strongest on lightly labelled reads (3.6× at ≥ 12; protected fraction 0.49 vs 0.43,
  the missing b_m term). Still 2.5× on mid-b_m and 1.7× on heavily labelled reads.
- **Not homopolymers.** Runs mostly inside A/T tracts (≥ 5 bp) are only 1.4× over-represented at
  ≥ 12; runs outside tracts are 3.0×, and almost all runs are outside tracts.
- **What follows.** A run of 3-4 protected calls spans about 8-12 bp. The open state has one protection
  probability, already pinned at its band limit, so the cheapest way to explain a cluster is the shortest
  protected segment the model has: an 8 bp TF. The prior then learns where clusters occur and shortens
  TF durations to the minimum
  ([duration-shaping note](fiber-tf-nubs-spacer-and-duration-shaping.md)).

**What the clusters are is unresolved:** transient or partial protection (small proteins, nucleosome
breathing), m6A deposition correlated along the molecule rather than independent per site, or per-read
labelling variation (part of the light-read excess).

**Design implications (untested):**
- An open sub-state for short partial protection, with no TF meaning, so clusters do not need TF.
- Smooth protected-length distributions instead of hard caps (a soft TF tail past 40 bp, or TF +
  other_long as one block with a learned length distribution).
- Per-read labelling efficiency for the light-read component.

# Key Points

| run of protected calls | ≥ 3 | ≥ 4 | ≥ 5 | ≥ 6 | ≥ 8 | ≥ 10 | ≥ 12 |
|---|---|---|---|---|---|---|---|
| observed / expected, all reads | 1.02 | 1.13 | 1.28 | 1.47 | 1.84 | 2.13 | 2.37 |
| b_m < −0.4 | 1.21 | 1.40 | 1.63 | 1.93 | 2.53 | 3.04 | 3.62 |
| b_m −0.4 to 0.25 | 1.05 | 1.16 | 1.31 | 1.51 | 1.90 | 2.25 | 2.49 |
| b_m > 0.25 | 0.89 | 0.96 | 1.06 | 1.18 | 1.43 | 1.57 | 1.74 |
| runs mostly in A/T tracts (≥ 5 bp) | 0.89 | 0.92 | 0.96 | 1.06 | 1.25 | 1.38 | 1.42 |
| runs outside A/T tracts | 1.04 | 1.17 | 1.36 | 1.60 | 2.12 | 2.57 | 3.02 |

- **Viterbi TF length shape** (interior segments; figures in the duration-shaping note):

  | model | n | exactly 8 bp | exactly 40 bp |
  |---|---|---|---|
  | GELU ×1 | 60,187 | 54% | — |
  | GELU ×0 | 4,969 | about 0% | 48% (rest a 13-35 bp hump) |
  | collapsed (constant prior) | 9,597 | 46% | — |

# Details

```bash
T=/scratch/users/diamant/fiber_eval_work/tf_mechanism
MODELS=gelu_s0 sbatch $T/viterbi_lengths.sbatch; MODELS=collapsed CODE_DIR=$T/smf_hmm_head sbatch $T/viterbi_lengths.sbatch
python $T/plot_viterbi_lengths.py
python $T/open_overdispersion.py
```

- **Expected runs.** The same called sites drawn independently 20 times with each site's predicted
  open protection, so context-driven clustering of poorly labelled sites is already in the expectation.

# Related Notes

- [TF nub mechanism overview](fiber-tf-nub-mechanism-overview.md): All hypotheses and fix trials.
- [Spacer and duration shaping](fiber-tf-nubs-spacer-and-duration-shaping.md): The Viterbi length
  histograms and the 9th-position penalty.
- [Under-methylated reads](fiber-tf-nubs-undermethylated-reads.md): The per-read labelling part.

# Open Questions

- Is the clustering biological (partial protection) or methodological (correlated m6A deposition along
  the molecule)? A test would be the autocorrelation of m6A calls in fully accessible regions (e.g. ATAC
  peak centres), at fixed context.
- Does an open short-patch state remove the 8 bp nubs while keeping motif footprints?

# Sources

- Scripts and outputs in the frontmatter.
