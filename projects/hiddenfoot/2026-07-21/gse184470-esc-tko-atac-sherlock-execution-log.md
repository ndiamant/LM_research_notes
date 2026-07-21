---
title: GSE184470 ESC TKO ATAC Sherlock execution log
date: 2026-07-21
project: hiddenfoot
agent: Codex
status: draft
sources:
  - /home/users/diamant/repos/HiddenFoot/GSE184470_ESC_TKO_ATAC_Sherlock_agent_plan.md
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/scripts/common.sh
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/logs/submitted_jobs.tsv
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/preflight.txt
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_1.raw_read_pairs.tsv
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_2.raw_read_pairs.tsv
tags:
  - hiddenfoot
  - atac-seq
  - gse184470
  - sherlock
  - chrombpnet
  - workflow-log
---

# Summary

Started execution of the `GSE184470` ESC DNMT-TKO ATAC-seq processing plan on Sherlock. The workflow root is `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC`, with all scripts, logs, references, intermediates, and outputs kept under that scratch directory.

As of the latest check on 2026-07-21, preflight, reference preparation, SRA download, and trimming had completed successfully. Both replicate alignment/filtering jobs were running on the `akundaje` partition, and merge, track generation, peak calling, and finalization were still dependency-held.

# Key Points

- Created the scratch layout: `scripts`, `logs`, `ref`, `sra`, `fastq`, `trimmed`, `bam`, `peaks`, `tracks`, `qc`, `tmp`, and `checksums`.
- Verified required Sherlock modules were available: SRA Toolkit 3.0.7, pigz 2.4, cutadapt 5.1, bowtie2 2.5.4, samtools 1.16.1, bedtools 2.30.0, UCSC `bedGraphToBigWig`/`bigWigInfo`, MACS2 2.2.9.1, FastQC, and MultiQC.
- MACS3 was not available as a Sherlock module, so the workflow uses MACS2, consistent with the plan's fallback.
- Added caveat-driven workflow improvements before execution: explicit Slurm log paths, SRA home/cache under scratch, mtDNA pre-filter counts, FRiP on the reproducible peak set, and chromosome-edge filtering for centered 2,114 bp training windows.
- Submitted jobs first on `normal`, then moved active/pending workflow jobs to `stat,akundaje,hns` after Nate reported access to those freer queues.
- Patched the saved `.sbatch` scripts to use `#SBATCH --partition=stat,akundaje,hns` for future reruns.
- Patched `scripts/common.sh` to use `module purge` rather than `module --force purge`; a full force purge removed Sherlock's sticky base toolchain and broke clean loading of py312/py39 biology modules.

# Details

## Submitted Jobs

The workflow was submitted as a dependency chain:

| role | job id |
| --- | --- |
| preflight | `35001640` |
| reference | `35001642` |
| download | `35001645` |
| trim | `35001651` |
| align | `35001653` |
| merge | `35001716` |
| tracks | `35001802` |
| peaks | `35001804` |
| finalize | `35001806` |

The dependency order was:

1. Preflight.
2. Reference and SRA download after preflight.
3. Trim after download.
4. Align/filter after trim and reference.
5. Merge after both align tasks.
6. Tracks and peaks after merge.
7. Finalize after tracks and peaks.

The partition switch was applied in place with `scontrol update JobId=<job> Partition=stat,akundaje,hns` for jobs `35001642`, `35001645`, `35001651`, `35001653`, `35001716`, `35001802`, `35001804`, and `35001806`. After the switch, reference and both SRA download array tasks started immediately on `akundaje`.

## Completed Work

Preflight completed with exit code `0:0` and wrote `qc/preflight.txt`. It recorded scratch storage availability and module paths/versions. Scratch had about 100 TB available at workflow start.

Reference preparation completed with exit code `0:0`. The job downloaded UCSC `mm10.fa.gz` and `mm10.chrom.sizes`, built the Bowtie2 index, generated `mm10.primary.chrom.sizes`, downloaded and sorted `mm10-blacklist.v2.bed`, and wrote reference checksums.

SRA download completed for both runs:

| sample | run | raw read pairs |
| --- | --- | ---: |
| `ESC_TKO_1` | `SRR15974966` | 36,126,571 |
| `ESC_TKO_2` | `SRR15974967` | 40,641,484 |

Trim/FastQC completed for both runs and wrote cutadapt reports plus raw and trimmed FastQC outputs under `qc/`.

## Current State

At the latest Slurm check, both alignment/filtering array tasks were running on `akundaje`:

| job | state | partition | note |
| --- | --- | --- | --- |
| `35001653_0` | running | `akundaje` | `ESC_TKO_1` align/filter |
| `35001653_1` | running | `akundaje` | `ESC_TKO_2` align/filter |
| `35001716` | pending | `stat,akundaje,hns` | waits on align |
| `35001802` | pending | `stat,akundaje,hns` | waits on merge |
| `35001804` | pending | `stat,akundaje,hns` | waits on merge |
| `35001806` | pending | `stat,akundaje,hns` | waits on tracks and peaks |

## Script Set

The scratch workflow scripts are:

- `scripts/common.sh`: shared root, sample arrays, module loading helpers, and validation helpers.
- `scripts/00_preflight.sbatch`: storage, modules, tool paths, and version logging.
- `scripts/01_reference.sbatch`: mm10 reference/blacklist setup and validation.
- `scripts/02_download_fastq.sbatch`: `prefetch`, `fasterq-dump`, compression, checksums, and raw-pair counts.
- `scripts/03_trim_qc.sbatch`: cutadapt trimming plus optional FastQC.
- `scripts/04_align_filter.sbatch`: bowtie2 alignment, samtools duplicate removal, MAPQ/flag/chrM filtering, fragment-length files, BAM checksums.
- `scripts/bam_to_tn5_bigwig.sh`: unnormalized base-resolution Tn5 insertion BED, bedGraph, and BigWig generation with count validation.
- `scripts/05_merge.sbatch`: pooled filtered BAM generation.
- `scripts/06_tracks.sbatch`: replicate and pooled Tn5 BigWig array.
- `scripts/07_peaks.sbatch`: MACS2 peak calling, reproducible peak set, blacklist and chromosome-edge filtering, FRiP.
- `scripts/08_finalize.sbatch`: final validation, manifest, README, optional MultiQC.

# Related Notes

- [mESC data preparation and MAP decoder addition](../2026-07-20/mesc-data-preparation-and-map-decoder.md): Documents the mESC HiddenFoot data context that motivated keeping mm10 and accessibility-derived inputs consistent.
- [CTCF motif experiment summary](ctcf-experiment-summary.md): Records the active HiddenFoot mESC/CTCF experiment context surrounding this ATAC processing work.

# Open Questions

- Alignment/filtering was still running when this note was written; final BAM, Tn5 BigWig, peak, FRiP, and manifest outputs still need to be checked after dependent jobs finish.
- The final report should still validate filtered fragment counts, duplicate rates, fragment-length periodicity, peak counts, reproducible peak counts, and total Tn5 insertions in each BigWig.
- TSS enrichment QC was not implemented in this pass because no trusted annotation path was selected yet.

# Sources

- `/home/users/diamant/repos/HiddenFoot/GSE184470_ESC_TKO_ATAC_Sherlock_agent_plan.md`: Original execution plan reviewed and implemented.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/scripts/`: Generated Slurm workflow scripts.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/logs/submitted_jobs.tsv`: Submitted job ID mapping.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/preflight.txt`: Tool/module and storage preflight output.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_1.raw_read_pairs.tsv`: Raw read-pair count for replicate 1.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_2.raw_read_pairs.tsv`: Raw read-pair count for replicate 2.
