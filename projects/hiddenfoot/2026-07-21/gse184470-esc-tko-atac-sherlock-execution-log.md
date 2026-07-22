---
title: GSE184470 ESC TKO ATAC Sherlock execution log
date: 2026-07-21
project: hiddenfoot
agent: Codex
status: complete
sources:
  - /home/users/diamant/repos/HiddenFoot/GSE184470_ESC_TKO_ATAC_Sherlock_agent_plan.md
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/scripts/common.sh
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/logs/submitted_jobs.tsv
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/README.md
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/manifest.tsv
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/preflight.txt
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_1.raw_read_pairs.tsv
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_2.raw_read_pairs.tsv
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_1.filtered_counts.tsv
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_2.filtered_counts.tsv
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/tracks/ESC_TKO_pooled.tn5_counts.tsv
  - /scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_reproducible.peak_count.txt
tags:
  - hiddenfoot
  - atac-seq
  - gse184470
  - sherlock
  - chrombpnet
  - workflow-log
---

# Summary

Executed the `GSE184470` ESC DNMT-TKO ATAC-seq processing plan on Sherlock. The workflow root is `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC`, with all scripts, logs, references, intermediates, and outputs kept under that scratch directory.

All Slurm jobs completed successfully with exit code `0:0`. The final deliverables include filtered replicate BAMs, a pooled filtered BAM, unnormalized base-resolution Tn5 insertion BigWigs for each replicate and the pool, relaxed MACS2 peaks, a reproducible peak set, a manifest, README, checksums, FastQC/MultiQC reports, and QC summaries. The pooled training targets are `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/bam/ESC_TKO_pooled.filtered.bam`, `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/tracks/ESC_TKO_pooled.tn5_counts.bw`, and `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/peaks/ESC_TKO_reproducible.narrowPeak`.

# Key Points

- Created the scratch layout: `scripts`, `logs`, `ref`, `sra`, `fastq`, `trimmed`, `bam`, `peaks`, `tracks`, `qc`, `tmp`, and `checksums`.
- Verified required Sherlock modules were available: SRA Toolkit 3.0.7, pigz 2.4, cutadapt 5.1, bowtie2 2.5.4, samtools 1.16.1, bedtools 2.30.0, UCSC `bedGraphToBigWig`/`bigWigInfo`, MACS2 2.2.9.1, FastQC, and MultiQC.
- MACS3 was not available as a Sherlock module, so the workflow uses MACS2, consistent with the plan's fallback.
- Added caveat-driven workflow improvements before execution: explicit Slurm log paths, SRA home/cache under scratch, mtDNA pre-filter counts, FRiP on the reproducible peak set, and chromosome-edge filtering for centered 2,114 bp training windows.
- Submitted jobs first on `normal`, then moved active/pending workflow jobs to `stat,akundaje,hns` after Nate reported access to those freer queues.
- Patched the saved `.sbatch` scripts to use `#SBATCH --partition=stat,akundaje,hns` for future reruns.
- Patched `scripts/common.sh` to use `module purge` rather than `module --force purge`; a full force purge removed Sherlock's sticky base toolchain and broke clean loading of py312/py39 biology modules.
- Final BAMs passed `samtools quickcheck`; final narrowPeak files had 10 columns and nonnegative summit offsets.
- Tn5 insertion validation passed exactly: insertion BED counts, bedGraph sums, and saved TSV counts matched for each replicate and the pooled track.
- Final reproducible peak set contains 90,463 peaks after requiring overlap with both replicate peak sets, removing 1,057 bp padded blacklist overlaps, and filtering chromosome-edge windows for 2,114 bp centered model inputs.

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

## Final Outcome

All workflow jobs completed with exit code `0:0`:

| job | state | partition | note |
| --- | --- | --- | --- |
| `35001640` | completed | `normal` | preflight |
| `35001642` | completed | `akundaje` | reference preparation |
| `35001645_0`, `35001645_1` | completed | `akundaje` | SRA download |
| `35001651_0`, `35001651_1` | completed | `akundaje` | trimming and FastQC |
| `35001653_0`, `35001653_1` | completed | `akundaje` | alignment, deduplication, filtering |
| `35001716` | completed | `akundaje` | pooled BAM |
| `35001802_0`, `35001802_1`, `35001802_2` | completed | `akundaje` | Tn5 tracks |
| `35001804` | completed | `akundaje` | MACS2 peaks and reproducible peak set |
| `35001806` | completed | `akundaje` | manifest, README, MultiQC |

Key final metrics:

| metric | `ESC_TKO_1` | `ESC_TKO_2` | pooled |
| --- | ---: | ---: | ---: |
| raw read pairs | 36,126,571 | 40,641,484 | n/a |
| filtered alignments | 48,904,510 | 53,393,970 | 102,298,480 |
| filtered fragments | 24,452,255 | 26,696,985 | 51,149,240 |
| duplicate rate | 18.58% | 20.89% | n/a |
| Tn5 insertions | 48,904,510 | 53,393,970 | 102,298,480 |
| relaxed peaks | 109,808 | 115,028 | 127,207 |
| FRiP on reproducible peaks | 0.4977 | 0.4856 | 0.4914 |

The reproducible peak set contains 90,463 peaks.

Fragment-length summaries looked consistent between replicates and ATAC-like: short-fragment modes were around 44-54 bp, about 31-32% of filtered fragments were under 100 bp, and both replicate libraries retained mono- and di-nucleosome-range fragments. Replicate means were 222.4 bp and 224.0 bp.

Primary deliverable paths:

- Pooled BAM: `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/bam/ESC_TKO_pooled.filtered.bam`
- Pooled Tn5 BigWig: `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/tracks/ESC_TKO_pooled.tn5_counts.bw`
- Reproducible peak set: `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/peaks/ESC_TKO_reproducible.narrowPeak`
- MultiQC report: `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/multiqc_report.html`
- Manifest: `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/manifest.tsv`
- README: `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/README.md`

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

- TSS enrichment QC was not implemented in this pass because no trusted annotation path was selected.
- MACS2 was used because MACS3 was unavailable as a Sherlock module. This is expected from the plan's fallback, but it should be remembered when comparing against MACS3-based pipelines.
- The preflight job stderr contains module/version-probe warnings from mixing py312 tools with MACS2/py39 in one process. The downstream workflow scripts avoided this by using cleaner module-loading paths, and all downstream jobs completed successfully.

# Sources

- `/home/users/diamant/repos/HiddenFoot/GSE184470_ESC_TKO_ATAC_Sherlock_agent_plan.md`: Original execution plan reviewed and implemented.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/scripts/`: Generated Slurm workflow scripts.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/logs/submitted_jobs.tsv`: Submitted job ID mapping.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/README.md`: Final run summary, output paths, counts, FRiP values, and Slurm job IDs.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/manifest.tsv`: Checksummed manifest for final files and QC artifacts.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/preflight.txt`: Tool/module and storage preflight output.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_1.raw_read_pairs.tsv`: Raw read-pair count for replicate 1.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/ESC_TKO_2.raw_read_pairs.tsv`: Raw read-pair count for replicate 2.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/*.filtered_counts.tsv`: Filtered read and fragment counts.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/tracks/*.tn5_counts.tsv`: Tn5 insertion count validation summaries.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/*.peak_count.txt`: Peak counts for replicate, pooled, and reproducible peak sets.
- `/scratch/users/diamant/smf_data/GSE184470_ESC_TKO_ATAC/qc/*.frip_on_reproducible.tsv`: FRiP values against the final reproducible peak set.
