---
title: CTCF motif and local-control ChIP–nexus signal collection
date: 2026-09-02
project: smf-hmm
agent: Codex
status: complete
sources:
  - https://jaspar.elixir.no/matrix/MA0139.2/
  - https://github.com/Boyle-Lab/Blacklist/releases/tag/v2.0
  - /scratch/groups/btrippe/ndiamant/mm10.fa
  - /scratch/groups/btrippe/ndiamant/mm10-blacklist.v2.bed.gz
  - /scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_*.bw
tags:
  - ctcf
  - chip-nexus
  - motif-scanning
  - matched-controls
  - mesc
  - evaluation
---

# Summary

We designed and ran a genome-wide mm10 scan of the external JASPAR CTCF motif
MA0139.2, followed by extraction of strand-paired CTCF ChIP–nexus signal at
each retained motif and 100 nearby shifted controls. This analysis establishes
whether canonical CTCF motif positions carry greater absolute ChIP–nexus signal
than comparable locations before any single-molecule footprint caller enters
the evaluation.

The primary score is the positive count plus the absolute negative count in a
101 bp window centered on the motif or control. The workflow does not divide
by local-flank signal. It retains individual control values as well as compact
per-motif summaries.

# Key Points

- The scan uses MA0139.2 at FIMO `p <= 1e-4` on mm10 primary chromosomes
  `chr1`–`chr19`, `chrX`, and `chrY`.
- Motif and control signal windows are excluded if they intersect the mm10 v2
  blacklist.
- Each motif receives 100 distinct controls shifted 200–2,000 bp in either
  direction using seed `20260902`.
- A control window may not intersect any retained MA0139.2 motif.
- Signal is collected separately for all seven datasets and both signed
  bigWig strands. CPM normalization uses the positive-track sum plus the
  absolute negative-track sum for the corresponding dataset.
- FIMO returned 1,402,215 sites before filtering; 76,294 blacklist-overlapping
  sites were removed, leaving 1,325,921 motifs and 132,592,100 controls.
- Mean motif CPM exceeded mean matched-control CPM in every dataset by
  2.99-fold to 8.59-fold. This is descriptive rather than an inferential test.
- Bulk ChIP–nexus provides a locus-level reference; it cannot validate an
  individual single-molecule binding event.

# Details

## Inputs and checksums

```text
b7866a8a5c7392f25c3ab818375c143feb15c1a6fcaf649fadf5a33c797c1a59  MA0139.2.meme
9a9559cbdf9e6775f0dfbd1f42f17d6b01a8fe050833a66d6aeb5bb7f3d09a74  mm10.fa
febafb843c6df492f9a9fc418f8796762ee899d9864330fb509ae2d38ddc0b46  mm10-blacklist.v2.bed.gz
```

The 14 bigWigs form seven pairs named
`mesc_ctcf_nexus_{1..7}_{positive,negative}.bw`. Their chromosome sizes agree
with the supplied mm10 FASTA. Positive tracks contain nonnegative signal and
negative tracks contain nonpositive signal.

## Durable workflow

Workflow files are in the SMF repository:

```text
/home/users/diamant/repos/smf_hmm/evaluation/datasets/ctcf_mESC/README.md
/home/users/diamant/repos/smf_hmm/evaluation/datasets/ctcf_mESC/run_motif_chip_signal.sbatch
/home/users/diamant/repos/smf_hmm/evaluation/datasets/ctcf_mESC/scripts/write_primary_fasta.py
/home/users/diamant/repos/smf_hmm/evaluation/datasets/ctcf_mESC/scripts/collect_motif_chip_signal.py
/home/users/diamant/repos/smf_hmm/evaluation/datasets/ctcf_mESC/scripts/validate_motif_chip_signal.py
```

Large inputs and outputs are rooted at:

```text
/scratch/groups/btrippe/ndiamant/ctcf_chip_nexus
```

The downloaded motif and isolated MEME environment were created with:

```bash
root=/scratch/groups/btrippe/ndiamant/ctcf_chip_nexus
mkdir -p "$root"/{references,scan,results,logs,envs}

curl -fsSL \
  'https://jaspar.elixir.no/api/v1/matrix/MA0139.2.meme' \
  -o "$root/references/MA0139.2.meme"

mamba create -y \
  -p "$root/envs/meme-5.5.7" \
  --override-channels \
  -c conda-forge \
  -c bioconda \
  meme=5.5.7

conda list --explicit \
  -p "$root/envs/meme-5.5.7" \
  > "$root/references/meme-5.5.7-linux-64-explicit.txt"
```

The complete scan and extraction command is preserved in
`run_motif_chip_signal.sbatch`. It was submitted from the SMF repository with:

```bash
cd /home/users/diamant/repos/smf_hmm
sbatch evaluation/datasets/ctcf_mESC/run_motif_chip_signal.sbatch
```

The first submission, job `41775314`, failed before computation because
`conda activate smf_clean` referenced unset variable `ADDR2LINE` while Bash
`nounset` was active. The launcher was changed from `set -euo pipefail` to
enable `set -u` only after Conda activation. Recovery job `41775457` was then
submitted with the same `sbatch` command and ran on the `stat` partition.

Job `41775457` was cancelled after FIMO reported that the
`--max-stored-scores 1000000` limit had been reached and hits with
`p >= 4.2e-05` were being discarded. The cap was increased to 5,000,000 and
the launcher was changed to reuse completed FASTA/background files. Final job
`41779679` ran on `stat` and completed successfully:

```text
elapsed: 01:35:21
total CPU: 01:34:25
maximum RSS: 11,785,436 KB
exit code: 0:0
node: sh04-15n16
```

## Design and output schema

FIMO coordinates are converted from one-based inclusive coordinates to
zero-based half-open BED coordinates. The motif center is
`(bed_start + bed_end) // 2`, and the signal interval is
`[center - 50, center + 51)`.

The HDF5 result contains:

```text
dataset_names
motif_center
control_center
control_shift
library_total
motif_positive_counts
motif_negative_abs_counts
control_positive_counts
control_negative_abs_counts
```

The first dimension of every signal array follows `ctcf_motif_sites.tsv.gz`.
Control arrays have a second dimension of length 100, and signal arrays have a
final dimension of length seven. Missing bigWig values are interpreted as zero.

## Results and validation

Final outputs are:

```text
/scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/results/ctcf_motif_sites.tsv.gz
/scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/results/ctcf_motif_local_controls.h5
/scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/results/ctcf_motif_chip_summary.tsv.gz
/scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/results/validation.json
/scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/results/SHA256SUMS
```

Validation reported:

```text
FIMO sites before blacklist: 1,402,215
Blacklist-overlapping sites removed: 76,294
Chromosome-boundary sites removed: 0
Retained motifs: 1,325,921
Controls per motif: 100
Total controls: 132,592,100
Absolute shift range: 200–2,000 bp
All controls unique within motif: true
All signal finite: true
All signal nonnegative after taking abs(negative): true
```

An independent exhaustive coordinate audit was run after the Slurm job:

```bash
cd /home/users/diamant/repos/smf_hmm
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate smf_clean

python evaluation/datasets/ctcf_mESC/scripts/validate_motif_chip_signal.py \
  --sites /scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/results/ctcf_motif_sites.tsv.gz \
  --values /scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/results/ctcf_motif_local_controls.h5 \
  --blacklist /scratch/groups/btrippe/ndiamant/mm10-blacklist.v2.bed.gz \
  --chrom-sizes /scratch/groups/btrippe/ndiamant/ctcf_chip_nexus/references/mm10.primary.chrom.sizes \
  --half-window 50 \
  --min-shift 200 \
  --max-shift 2000
```

It validated all 1,325,921 motifs and 132,592,100 controls. In particular,
control centers reproduce the stored shifts, complete signal windows remain
inside chromosome bounds, and no control signal window intersects a blacklist
interval or retained motif. Both gzip files passed `gzip -t`, all four result
checksums passed `sha256sum -c`, and final-job stderr contained no FIMO
score-cap warning.

Output checksums are:

```text
a1dadec530b911ff68dbdd680af2aa1285bb3c0a9ba7f2e4a7e2a3cde836cecd  ctcf_motif_local_controls.h5
dc6c3818e9ad24cd4671b11b0a6e4a1ac0b9144d62230b004b527a367c39dd69  ctcf_motif_sites.tsv.gz
31266ef417c13ad62f4733083a63a93a9cace3eea6c742f5b216284f74ffb6e9  ctcf_motif_chip_summary.tsv.gz
81bc3635d23aeeb822af7e2764a5b5eb73c3bbfbec4253fa0f77246accee6b74  validation.json
```

### Descriptive motif-versus-control check

The full `p <= 1e-4` scan is deliberately permissive. Nevertheless, mean
motif CPM exceeded the mean of the 100 matched-control CPMs in all datasets:

| Dataset | Mean motif CPM | Mean control CPM | Ratio | Mean control percentile |
| --- | ---: | ---: | ---: | ---: |
| 1 | 0.2153 | 0.03743 | 5.75 | 0.619 |
| 2 | 0.1062 | 0.03556 | 2.99 | 0.580 |
| 3 | 0.1437 | 0.03745 | 3.84 | 0.611 |
| 4 | 0.2357 | 0.04058 | 5.81 | 0.652 |
| 5 | 0.3338 | 0.03997 | 8.35 | 0.649 |
| 6 | 0.3305 | 0.03845 | 8.59 | 0.672 |
| 7 | 0.3269 | 0.03916 | 8.35 | 0.653 |

ChIP support increased monotonically with motif strength:

| FIMO p-value cutoff | Retained sites | Median combined control percentile | Fraction at or above 0.95 |
| --- | ---: | ---: | ---: |
| `1e-4` | 1,325,921 | 0.640 | 0.188 |
| `5e-5` | 753,437 | 0.710 | 0.279 |
| `1e-5` | 215,929 | 0.935 | 0.487 |
| `1e-6` | 37,595 | 1.000 | 0.829 |
| `1e-7` | 10,726 | 1.000 | 0.970 |

There were 11,741 retained sites with FIMO `q <= 0.05`; their median combined
control percentile was 1.0. These summaries demonstrate the expected signal,
but formal uncertainty should use region- or chromosome-level resampling rather
than treating nearby motif sites as independent.

# Related Notes

- [CTCF motif experiment summary](../../hiddenfoot/2026-07-21/ctcf-experiment-summary.md): Shows that motif selection and thresholds previously changed which CTCF-family states HiddenFoot could express.
- [FootprintCharter baseline integration](../2026-08-18/footprintcharter-baseline-integration.md): Documents one of the state-calling baselines that will eventually be evaluated against the CTCF reference.

# Open Questions

- Should the final SMF benchmark be restricted to MA0139.2, or include the
  extended CTCF motifs MA1929.2 and MA1930.2 as a sensitivity analysis?
- Which SMF region set should be applied to the genome-wide motif table for the
  first state-caller comparison?
- Should `chrY` be retained once the sex and provenance of the mESC data are
  confirmed?

# Sources

- [JASPAR MA0139.2](https://jaspar.elixir.no/matrix/MA0139.2/)
- [Blacklist v2 release](https://github.com/Boyle-Lab/Blacklist/releases/tag/v2.0)
- `/scratch/groups/btrippe/ndiamant/mm10.fa`
- `/scratch/groups/btrippe/ndiamant/mm10-blacklist.v2.bed.gz`
- `/scratch/groups/btrippe/ndiamant/julia_lab/bw/mesc_ctcf_nexus_*.bw`
- `/home/users/diamant/repos/smf_hmm/evaluation/evaluation_directions.md`
