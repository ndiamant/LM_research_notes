---
title: Fixed-HMM posterior-sampling baseline comparison
date: 2026-09-04
project: smf-hmm
agent: Codex
status: complete
sources:
  - /scratch/users/diamant/ctcf_nexus_6_baseline_comparison.pdf
  - /scratch/users/diamant/simple.chr1.ctcf_nexus_6_test
  - /scratch/users/diamant/geometric_hmm.chr1.ctcf_nexus_6_test
  - /scratch/users/diamant/phase_hmm.chr1.ctcf_nexus_6_test
  - /home/users/diamant/repos/smf_hmm/baselines/fixed_hmm
  - /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/compare_runs.py
tags:
  - smf-hmm
  - baseline
  - posterior-sampling
  - geometric-hmm
  - phase-type-hmm
  - ctcf
  - chip-nexus
  - evaluation
---

# Summary

Three chr1 SMF state callers were compared against mESC CTCF ChIP–nexus
replicate 6 using the MA0139.2 motif benchmark: the original linear
interpolation/run-length caller, a fixed three-state geometric HMM with one
posterior path sample per molecule, and a fixed 138-state phase-type HMM with
minimum nucleosome and TF durations.

The **phase-type HMM was strongest**, with Spearman correlation 0.479 across
10,009 adequately covered motifs and 0.554 in the benchmark's
same-track-selected ChIP-supported subset. Its overall correlation exceeded
interpolation by 0.034 and the geometric HMM by 0.054. Both differences
remained above zero in a paired 1 Mb genomic-block bootstrap.

The ranking is encouraging but provisional. This comparison covers one
chromosome, one ChIP replicate, and one posterior-sampling seed. It also uses
the legacy fixed-region input with synthetic molecule identifiers documented
in the [interpolation baseline note](simple-interpolation-ctcf-chip-benchmark.md).

# Key Points

- Overall Spearman correlations were 0.445 for interpolation, 0.425 for the
  geometric HMM, and 0.479 for the phase-type HMM.
- ChIP-supported correlations were 0.473, 0.497, and 0.554, respectively.
- A 500-replicate paired bootstrap over 189 genomic 1 Mb blocks gave a
  phase-type-minus-interpolation difference interval of approximately
  `[0.018, 0.051]` and a phase-type-minus-geometric interval of
  `[0.039, 0.068]`.
- The geometric HMM compressed inferred occupancy: mean occupancy was 0.063,
  compared with 0.111 for interpolation and 0.093 for the phase-type HMM.
- Both posterior-sampling HMMs nearly eliminated sites with exactly zero
  occupancy (0.3–0.4%), versus 8.3% for interpolation.
- Median usable molecule coverage was 153 for both HMMs and 142 for
  interpolation. The difference arises because interpolation represents
  41–99 bp protected runs as ambiguous (`.`), excluding overlapping motifs
  from molecule coverage.
- The phase-type HMM tracked interpolation more closely than the geometric
  HMM in Pearson scale, but all methods differed materially in rank and
  dynamic range.
- The geometric run was incorrectly recorded as method
  `simple_interpolation` in both `run.json` and `metrics.tsv`. The comparison
  report labels it from its run directory; future runs should correct the
  `--method-name` argument.

# Details

## Methods compared

The fixed-HMM implementation is in:

```text
/home/users/diamant/repos/smf_hmm/baselines/fixed_hmm
```

Both HMM callers use fixed parameters and forward-filtering backward-sampling
to draw one latent path per molecule. Random keys are derived from the global
seed and molecule ID, making a molecule's path independent of row ordering and
batch composition. Internal missing positions are inferred by the model;
unobserved sequence outside the first and last informative calls is emitted as
`.`.

The geometric model has accessible, nucleosome, and TF states with geometric
durations. The phase-type model expands nucleosome and TF into deterministic
minimum-length chains followed by geometric tails. Its default nucleosome
duration has minimum 130 bp and mean 150 bp; its TF duration has minimum 7 bp
and mean 10 bp. Both models use stationary probabilities 0.28 accessible,
0.70 nucleosome, and 0.02 TF, with methylation probabilities 0.95 in accessible
DNA and 0.025 in protected DNA.

## Comparable benchmark inputs

The three run manifests used identical benchmark configuration, motif sites,
ChIP values, FASTA, blacklist, and coverage thresholds. The state-call inputs
were:

```text
/scratch/users/diamant/simple.chr1.states.tsv.gz
/scratch/users/diamant/geometric.chr1.states.tsv.gz
/scratch/users/diamant/phase-type.chr1.states.tsv.gz
```

The evaluator used MA0139.2 at motif FPR `1e-4`, a 101 bp ChIP window,
minimum molecule coverage 10, and the paired tracks:

```text
/scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_6_positive.bw
/scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_6_negative.bw
```

All three runs scanned 92,613 chr1 motifs, retained 10,009 adequately covered
sites, and selected 7,331 sites above the same-track ChIP-support threshold
0.02970297.

## Primary results

| Method | Overall rho | ChIP-supported rho | Mean occupancy | Median occupancy | Zero occupancy | Median molecules |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Simple interpolation | 0.445 | 0.473 | 0.111 | 0.052 | 8.3% | 142 |
| Geometric HMM | 0.425 | 0.497 | 0.063 | 0.047 | 0.3% | 153 |
| Phase-type HMM | **0.479** | **0.554** | 0.093 | 0.047 | 0.4% | 153 |

The geometric caller was less correlated overall than interpolation, with a
paired bootstrap difference of -0.020 and 95% interval approximately
`[-0.040, -0.002]`. Nevertheless, its correlation became higher than
interpolation after removing the bottom ChIP-signal quartile. The phase-type
caller was strongest in both subsets.

## ChIP calibration and dynamic range

Mean occupancy in the highest ChIP-signal decile was 0.499 for interpolation,
0.182 for the geometric HMM, and 0.408 for the phase-type HMM. In the ninth
decile, the corresponding values were 0.172, 0.088, and 0.147. Lower deciles
were much closer, generally between 0.04 and 0.06.

Thus, all methods primarily separate the strongest bulk CTCF sites, but the
geometric HMM has substantially less dynamic range. The phase-type duration
constraints restore much of the high-signal occupancy response while
improving rank correlation beyond interpolation.

## Pairwise agreement

| Comparison | Rank rho | Pearson r | Mean absolute occupancy difference |
| --- | ---: | ---: | ---: |
| Interpolation vs geometric | 0.472 | 0.873 | 0.070 |
| Interpolation vs phase-type | 0.557 | 0.923 | 0.045 |
| Geometric vs phase-type | 0.587 | 0.887 | 0.045 |

High Pearson agreement alongside more modest rank agreement reflects shared
large-scale occupancy patterns but meaningful differences among sites with
similar or low occupancy.

## Block-bootstrap sensitivity

The comparison report resampled chr1 in 1 Mb genomic blocks rather than
resampling individual motifs. There were 189 represented blocks and 500
paired replicates. Approximate 95% intervals for overall rho were:

| Method | Point estimate | 95% block-bootstrap interval |
| --- | ---: | ---: |
| Simple interpolation | 0.445 | `[0.427, 0.463]` |
| Geometric HMM | 0.425 | `[0.403, 0.448]` |
| Phase-type HMM | 0.479 | `[0.461, 0.498]` |

The 1 Mb block size is a sensitivity choice, not a fully specified inferential
model. These intervals demonstrate that the phase-type improvement is not
driven by a small number of individual motifs, but they should not be treated
as definitive uncertainty estimates.

## PDF report reproduction

The five-page comparison report is:

```text
/scratch/users/diamant/ctcf_nexus_6_baseline_comparison.pdf
```

It was generated from the SMF repository with:

```bash
cd /home/users/diamant/repos/smf_hmm

MPLCONFIGDIR=/home/users/diamant/repos/smf_hmm/.mplconfig \
/home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python \
  -m evaluation.tf_chip.compare_runs \
  --run 'Simple interpolation=/scratch/users/diamant/simple.chr1.ctcf_nexus_6_test' \
  --run 'Geometric HMM=/scratch/users/diamant/geometric_hmm.chr1.ctcf_nexus_6_test' \
  --run 'Phase-type HMM=/scratch/users/diamant/phase_hmm.chr1.ctcf_nexus_6_test' \
  --output /scratch/users/diamant/ctcf_nexus_6_baseline_comparison.pdf \
  --bootstrap-replicates 500 \
  --block-size 1000000
```

The PDF contains headline metrics, ChIP-decile calibration, occupancy and
coverage distributions, site-level ChIP comparisons, pairwise occupancy
agreement, block-bootstrap intervals, and interpretation caveats.

# Related Notes

- [Simple interpolation baseline against chr1 CTCF ChIP–nexus](simple-interpolation-ctcf-chip-benchmark.md): Establishes the initial baseline result, legacy-region overlap audit, and input-provenance caveats reused here.
- [CTCF motif and local-control ChIP–nexus signal collection](../2026-09-02/ctcf-motif-chip-nexus-local-controls.md): Documents construction and validation of the MA0139.2/mm10 ChIP reference; replicate 6 had the strongest motif-versus-control enrichment.
- [FootprintCharter baseline integration and missingness limitations](../2026-08-18/footprintcharter-baseline-integration.md): Records the earlier external baseline and motivates missing-aware state callers.

# Open Questions

- How much do the HMM benchmark correlations vary across posterior-sampling
  seeds?
- Does the phase-type advantage replicate across ChIP–nexus datasets 1–7 and
  across all chromosomes?
- How should uncertainty account jointly for genomic correlation, overlapping
  legacy regions, and posterior-path sampling?
- Does the phase-type advantage persist on current full-alignment SMF call
  files with BAM QNAME molecule identity?
- Why does posterior sampling nearly eliminate zero-occupancy sites, and would
  posterior motif-event probabilities provide a better evaluator input than a
  single sampled hard path?

# Sources

- `/scratch/users/diamant/ctcf_nexus_6_baseline_comparison.pdf`: Five-page visual comparison report.
- `/scratch/users/diamant/simple.chr1.ctcf_nexus_6_test`: Interpolation benchmark run.
- `/scratch/users/diamant/geometric_hmm.chr1.ctcf_nexus_6_test`: Geometric-HMM benchmark run; recorded method name is incorrect.
- `/scratch/users/diamant/phase_hmm.chr1.ctcf_nexus_6_test`: Phase-type-HMM benchmark run.
- `/home/users/diamant/repos/smf_hmm/baselines/fixed_hmm`: Fixed-HMM model definitions, FFBS implementation, CLIs, and tests.
- `/home/users/diamant/repos/smf_hmm/evaluation/tf_chip/compare_runs.py`: Reproducible PDF comparison generator and block-bootstrap analysis.
- `/scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/references/MA0139.2.meme`: CTCF motif used by the evaluator.
- `/scratch/groups/btrippe/ndiamant/mm10.fa`: Reference genome.
- `/scratch/groups/btrippe/ndiamant/mm10-blacklist.v2.bed.gz`: Excluded genomic regions.
