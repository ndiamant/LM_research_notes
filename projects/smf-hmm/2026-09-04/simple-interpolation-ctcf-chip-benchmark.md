---
title: Simple interpolation baseline against chr1 CTCF ChIP–nexus
date: 2026-09-04
project: smf-hmm
agent: Codex
status: complete
sources:
  - /home/users/diamant/repos/smf_hmm/baselines/simple.py
  - /scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions_TSVs/chr1.tsv
  - /scratch/users/diamant/simple.chr1.states.tsv.gz
  - /scratch/users/diamant/simple.chr1.ctcf_nexus_6_test
  - /scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_6_positive.bw
  - /scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_6_negative.bw
tags:
  - smf-hmm
  - baseline
  - linear-interpolation
  - ctcf
  - chip-nexus
  - mesc
  - evaluation
---

# Summary

An intentionally simple footprint-calling baseline was applied to the legacy
ONT mESC chr1 region TSV and evaluated against mESC CTCF ChIP–nexus replicate
6. In adequately covered MA0139.2 sites, inferred TF occupancy had Spearman
correlation **0.445** with ChIP signal across 10,009 sites. The benchmark's
same-track-selected ChIP-supported subset had correlation **0.473** across
7,331 sites.

This is a promising proof of concept: the simple caller clearly distinguishes
strongly ChIP-supported CTCF sites, and the association remains after a
diagnostic adjustment for motif score and molecule coverage. It is not yet a
formal baseline result because the legacy input consists of fixed region
segments with synthetic region/read identifiers rather than full-alignment,
BAM-QNAME molecules required by the current SMF call specification.

# Key Points

- The baseline linearly interpolates protectedness between informative SMF
  positions, thresholds at 0.5, and classifies protected runs by length.
- Protected runs of at most 40 bp are `T`, runs of at least 100 bp are `N`,
  and 41–99 bp runs are ambiguous (`.`). Accessible sequence is `F`, and
  missing flanks are not extrapolated.
- The full adequately covered set gave Spearman rho 0.445; the ChIP-supported
  subset gave rho 0.473. Reported p-values of `0.0` are floating-point
  underflow at this sample size, not literal zero probabilities.
- Coverage was strong wherever legacy regions overlapped motifs: median 142
  usable molecules per motif and a minimum of 35, despite the configured
  minimum of 10.
- Median inferred occupancy was 5.2%, mean occupancy was 11.1%, and 8.3% of
  evaluated sites had zero inferred TF occupancy.
- Mean occupancy by ascending ChIP quintile was 4.4%, 5.4%, 5.9%, 6.0%, and
  33.5%. The relationship is therefore concentrated at strongly bound sites.
- Spearman rho was 0.174 within the bottom 80% of ChIP signal and 0.652 within
  the top 20%.
- A rank-residual diagnostic controlling jointly for motif score and molecule
  coverage retained correlation 0.367. This was not part of the formal
  evaluator output.
- Of 10,009 evaluated motifs, 208 lay in more than one legacy input region.
  Restricting to the 9,801 single-region motifs gave rho 0.444, so region
  overlap did not materially drive the headline correlation.

# Details

## Baseline definition

The implementation is:

```text
/home/users/diamant/repos/smf_hmm/baselines/simple.py
```

For each input molecule or legacy region/read row, calls are recoded as
accessible (`1 -> 0.0`), protected (`2 -> 1.0`), or missing (`0`). Linear
interpolation fills only positions between the first and last informative
calls. Interpolated protectedness greater than or equal to 0.5 is protected.
Contiguous protected intervals are labeled by their dense base-pair length:

| Length | Output state |
| --- | --- |
| <=40 bp | `T` |
| 41–99 bp | `.` |
| >=100 bp | `N` |

All remaining interpolated positions are `F`; leading and trailing missing
positions remain `.`. The script streams rows, accepts `+`, `-`, and `.`
strands, writes plain or gzip-compressed output, and reports molecule progress
with `tqdm`.

The chr1 state calls were generated from the repository root with:

```bash
cd /home/users/diamant/repos/smf_hmm

python baselines/simple.py \
  /scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions_TSVs/chr1.tsv \
  /scratch/users/diamant/simple.chr1.states.tsv.gz
```

## ChIP benchmark run

The evaluation used the paired positive and negative bigWigs for mESC CTCF
ChIP–nexus replicate 6, the MA0139.2 CTCF motif, mm10, and the mm10 v2
blacklist:

```bash
cd /home/users/diamant/repos/smf_hmm

/home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python \
  -m evaluation.tf_chip.benchmark \
  --state-calls /scratch/users/diamant/simple.chr1.states.tsv.gz \
  --chip-positive-bw /scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_6_positive.bw \
  --chip-negative-bw /scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_6_negative.bw \
  --pwm /scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/references/MA0139.2.meme \
  --fasta /scratch/groups/btrippe/ndiamant/mm10.fa \
  --blacklist /scratch/groups/btrippe/ndiamant/mm10-blacklist.v2.bed.gz \
  --output-dir /scratch/users/diamant/simple.chr1.ctcf_nexus_6_test \
  --method-name simple_interpolation \
  --chip-name mesc_ctcf_nexus_6
```

The run manifest records motif FPR `1e-4`, pseudocount 0.5, score threshold
7.950007, 101 bp ChIP windows, minimum molecule coverage 10, and a bottom
ChIP-signal quantile of 0.25. The resulting same-track support threshold was
0.02970297.

Outputs are:

```text
/scratch/users/diamant/simple.chr1.ctcf_nexus_6_test/metrics.tsv
/scratch/users/diamant/simple.chr1.ctcf_nexus_6_test/sites.tsv.gz
/scratch/users/diamant/simple.chr1.ctcf_nexus_6_test/run.json
```

## Primary results

| Benchmark subset | Sites | Spearman rho |
| --- | ---: | ---: |
| All adequately covered motifs | 10,009 | 0.4453 |
| ChIP-supported, same-track-selected | 7,331 | 0.4725 |

The evaluator scanned 92,613 non-blacklisted MA0139.2 hits on chr1. Exactly
10,009 (10.8%) had any usable molecule coverage, and all exceeded the
10-molecule adequacy threshold. Molecule coverage had median 142, mean 143.2,
minimum 35, 90th percentile 167, and maximum 541.

Occupancy had median 0.0520, mean 0.1105, 90th percentile 0.306, 95th
percentile 0.580, and maximum 0.954. There were 2,603 distinct occupancy
values. No site had occupancy exactly one. Plus- and minus-strand motifs were
balanced: 5,045 plus sites had mean occupancy 0.108, while 4,964 minus sites
had mean occupancy 0.113.

## ChIP-stratified behavior

Adequately covered sites were divided into equal-sized ChIP-signal quintiles
using first-ranked ties for the diagnostic below:

| ChIP quintile | Sites | Median ChIP signal | Mean occupancy | Median occupancy | Fraction with any TF call |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 2,002 | 0.0099 | 0.0444 | 0.0336 | 0.850 |
| 2 | 2,002 | 0.0396 | 0.0536 | 0.0400 | 0.897 |
| 3 | 2,001 | 0.0891 | 0.0593 | 0.0476 | 0.934 |
| 4 | 2,002 | 0.2475 | 0.0599 | 0.0507 | 0.944 |
| 5 | 2,002 | 5.0446 | 0.3354 | 0.2598 | 0.964 |

The baseline therefore separates the highest-signal sites strongly, but the
middle three quintiles have similar occupancy. This suggests weaker dynamic
resolution outside clearly occupied CTCF loci.

## Legacy input caveat

The input contains 1,105,548 rows across 5,461 regions, each exactly 2,114 bp
wide, with a median 201 rows per region. Identifiers have the form
`chr1:<region>:<index>`, rather than BAM QNAMEs, and 368 adjacent region pairs
overlap. It is therefore not possible from this file alone to establish
whether physical molecules are duplicated between overlapping regions.

Only 205 evaluated motifs were completely covered by two input regions and
three by three regions. Their mean nominal molecule coverages were 272 and
402, respectively. The remaining 9,801 single-region motifs retained rho
0.444, essentially the full result of 0.445. The overlap caveat matters for
formal sample-size interpretation, but it does not explain the observed
occupancy–ChIP association.

The current SMF format instead requires one row per full aligned molecule with
the BAM QNAME. A formal benchmark should be repeated on such data and should
include all intended chromosomes.

# Related Notes

- [CTCF motif and local-control ChIP–nexus signal collection](../2026-09-02/ctcf-motif-chip-nexus-local-controls.md): Documents the MA0139.2/mm10 reference construction and shows that replicate 6 had the strongest motif-versus-local-control ratio among the seven ChIP–nexus datasets.
- [FootprintCharter baseline integration and missingness limitations](../2026-08-18/footprintcharter-baseline-integration.md): Records the earlier external baseline, whose complete-case behavior performed poorly on sparse ONT SMF matrices.

# Open Questions

- Does the correlation replicate on current full-alignment SMF call files with
  BAM QNAME molecule identity?
- How stable is the result across all chromosomes and ChIP–nexus replicates
  1–7?
- Would a maximum interpolation-gap rule improve weak-site resolution, or
  merely reduce motif coverage?
- How sensitive is the result to the 40 bp TF and 100 bp nucleosome length
  thresholds?
- Should formal uncertainty use chromosome- or region-level resampling rather
  than treating nearby motif sites as independent?

# Sources

- `/home/users/diamant/repos/smf_hmm/baselines/simple.py`: Simple state-calling implementation.
- `/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions_TSVs/chr1.tsv`: Legacy chr1 SMF region input.
- `/scratch/users/diamant/simple.chr1.states.tsv.gz`: Generated state-call file.
- `/scratch/users/diamant/simple.chr1.ctcf_nexus_6_test/metrics.tsv`: Headline benchmark metrics.
- `/scratch/users/diamant/simple.chr1.ctcf_nexus_6_test/sites.tsv.gz`: Site-level occupancy, coverage, and ChIP signal used for diagnostics.
- `/scratch/users/diamant/simple.chr1.ctcf_nexus_6_test/run.json`: Benchmark configuration, input provenance, and package versions.
- `/scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_6_positive.bw`: Positive-strand ChIP–nexus signal.
- `/scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_6_negative.bw`: Negative-strand ChIP–nexus signal.
- `/scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/references/MA0139.2.meme`: CTCF motif used by the evaluator.
- `/scratch/groups/btrippe/ndiamant/mm10.fa`: Reference genome.
- `/scratch/groups/btrippe/ndiamant/mm10-blacklist.v2.bed.gz`: Excluded genomic regions.
