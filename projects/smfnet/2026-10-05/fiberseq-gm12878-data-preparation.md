---
title: Fiber-seq GM12878 data preparation for SMFNet, design choices and pipeline run
date: 2026-10-05
project: smfnet
agent: Claude Code
status: draft
sources:
  - SMFNet repo, data_preparation/fiberseq/ (README.md, select_regions.py, fiberbam_to_calls.py, fiberbam_to_calls.sbatch, calls_io.py, calls_to_h5.py), branch nd/more_assays, uncommitted at time of writing
  - SMFNet repo, src/smf_net/data.py and src/smf_net/evaluation/flank_baselines.py (fiberseq_m6a assay), same branch, uncommitted
  - smf_hmm repo, evaluation/datasets/archive/fiberseq_survey/FIBERSEQ_DATA_REPORT.md
  - smf_hmm repo, evaluation/datasets/smf/ (README.md, SPEC.md, samples.tsv, scripts/fiberbam_to_calls.py, scripts/make_regions.py)
  - smf_hmm repo, 20260420_create_h5_from_mBED.py (mouse ONT H5 builder)
  - /scratch/users/diamant/gch_rc_analysis/ (GCH strand symmetry scripts and results)
  - /scratch/users/diamant/fiberseq_strand/ (analysis scripts from this session)
  - User discussion, 2026-10-04 to 2026-10-05
tags: [smfnet, fiberseq, m6a, gm12878, data-preparation, strand-symmetry, rc-augmentation, region-selection, gc-matching]
---

# Summary

We built training data for a new SMFNet assay type, `fiberseq_m6a`, from the
released GM12878 FIRE Fiber-seq BAM. The output is two H5 files laid out like
the mouse ONT files: 84,125 active and 83,343 GC-matched background 2114 bp
windows, about 130 molecules each. The pipeline lives in
`data_preparation/fiberseq/` in the SMFNet repo and starts from the BAM, so it
is reproducible for publication. `data.py` gained the assay type, with RC
augmentation on by default. Training has not been run yet; the recommended
command is at the end of Details.

The main findings behind the design:

- Fiber-seq m6A is RC-symmetric to within noise, unlike our GCH data. That
  asymmetry is the one documented in the
  [flank baseline note](flank-baseline-logistic-regression.md).
- Mirroring the mouse bait design with the ENCODE cCRE registry gives a
  background with less composition shortcut than GC-matched random tiles.

# Key Points

- **Sample: GM12878** (Vollger et al. 2025, FIRE). Chosen for depth (about
  110–130 molecules per 2.1 kb window), clean strand symmetry, haplotype tags,
  and the most ENCODE orthogonal data. The alternatives:
  - COLO829BL is deeper but has about half the m6A rate.
  - Jurkat has about half the depth and a read-strand calling asymmetry.
  - Drosophila S2 (GSE249742) has naked-DNA and no-Hia5 controls. It is the
    natural second dataset for calibrating Hia5 sequence bias.
- **RC augmentation is valid for Fiber-seq.** Each HiFi molecule reports both
  strands: m6A at reference A positions for the top strand, at reference T
  positions for the bottom. On chr21, rates at A in context k match rates at T
  in rc(k):
  - ±1 bp context: mean |d| 0.0004, against a split-half noise floor of 0.0001.
  - ±3 bp context: corr 0.990, against a noise ceiling of 0.998.
  - Our GCH data, by contrast: mean |d| 0.034 against 0.0018 noise at ±1 bp,
    and corr 0.095 against 0.60 at ±3 bp.
  - Jurkat shows a read-strand effect (0.190 vs 0.181). It is RC-equivariant,
    so pooling strands remains valid.
- **RC for an A/T mark is a flip with no shift**, unlike the CpG/GpC case,
  where a call moves to the other cytosine of its dinucleotide.
- **Region design mirrors the mouse ONT data.** The mouse data uses the
  Sönmezer 2021 bait regions split by ATAC overlap. The human equivalent of
  the bait library is the full ENCODE cCRE registry (1.06M elements).
  - **Active:** windows overlapping an IDR ATAC peak.
  - **Background:** windows overlapping no ATAC peak, relaxed calls included.
    These are mostly `Low-DNase` cCREs, i.e. regulatory-looking DNA that is
    closed in GM12878.
- **Matching: GC in 1% bins, background downsampled 1:1** (a toggle, since the
  mouse set is not GC-matched).
  - Composition-only AUC (3-mers) separating active from background was 0.74
    for the cCRE split, against 0.77 for ATAC peaks versus GC-matched random
    tiles.
  - GC matching costs about 1–4% of windows and lowers the AUC to 0.71.
  - CpG o/e is the main remaining signal (AUC 0.67). We handle it in
    evaluation (the flank baseline), not by stricter matching, which would cost
    about 40% of peaks.
- **Extract from the BAM, not the existing `smf_hmm` TSVs.** The old TSVs only
  covered molecules near GM12878-active cCREs; 34% of ATAC peaks had essentially
  no reads. They were also built without the SNP fix below.
- **Converter fixes relative to the `smf_hmm` `fiberbam_to_calls.py`:**
  - A position is assayable only if the read base matches the reference A/T.
    The old version made "protected" calls at heterozygous SNPs. On two chr21
    windows the fix changed 3.2e-4 of bases, 71% of them at 29 recurrent
    mostly single-haplotype sites.
  - Adds the `HP` haplotype column.
  - Also excludes duplicate and QC-fail reads.
  - Handles explicit-mode `MM` segments.
  - Writes atomically.
- **Intermediate format: zstd-compressed TSV** (`-12 --long=27`), read and
  written by streaming through the `zstd` command-line tool. We kept the TSV
  step, rather than going straight from BAM to H5, for three reasons: it
  overlaps the HMM format, it keeps windowing out of the BAM stream, and it
  allows re-cropping in minutes. Long mode compresses 9–15× against about 6×
  for gzip. samtools stays a module dependency (`ml biology samtools/1.16.1`):
  pysam wheels need glibc ≥ 2.24 and Sherlock has 2.17.
- **Splits: ChromBPNet fold 0.**
  - Train: chr2, 4, 5, 7, 9–19, 21, 22, X, Y — 120,130 windows.
  - Val: chr8, chr20 — 12,304 windows.
  - Test: chr1, chr3, chr6 — 35,034 windows.

# Details

## Pipeline and run

All outputs are under `$GROUP_SCRATCH/ndiamant/fiberseq/GM12878/` (that is,
`/scratch/groups/btrippe/ndiamant/...`). The input BAM is
`$OAK/datasets/SMF_fiberseq/raw_bams/GM12878/GM12878.fire.bam` (326 GB). The
reference is `/scratch/groups/btrippe/brian/projects/SMF/fiberseq/reference/hg38.fa`.
The ENCODE inputs are listed in `data_preparation/fiberseq/README.md`. The hg38
blacklist v2 was downloaded to `$GROUP_SCRATCH/ndiamant/hg38-blacklist.v2.bed.gz`.

**Step 1, regions** (about 80 s, run on the interactive node):

```bash
F=/scratch/groups/btrippe/brian/projects/SMF/fiberseq
python data_preparation/fiberseq/select_regions.py \
    --ccre $F/encode/ccre/ENCFF733BFV_GM12878_cCRE_v3.bed.gz \
    --active-peaks $F/encode/GM12878_peaks/ENCFF614SMH.bed.gz \
    --exclude-peaks $F/encode/GM12878_peaks/ENCFF158ORB.bed.gz \
    --blacklist $GROUP_SCRATCH/ndiamant/hg38-blacklist.v2.bed.gz \
    --fasta $F/reference/hg38.fa \
    --out $GROUP_SCRATCH/ndiamant/fiberseq/GM12878/candidates.bed
```

1,062,844 cCREs → 1,058,400 after blacklist, N and edge filters → 583,549
after thinning (centres at least 1057 bp apart) → 84,241 active and 455,183
background candidates, with 44,125 dropped as ambiguous (relaxed peak only).

Window names follow the mouse ONT convention: `[a, a+2114)` is named
`chr:{a+1}-{a+2115}`, which `load_dna` reads correctly with the default
`dna_sub_one_from_end_idx=True`. This was verified on the ONT files and on 200
random GM12878 windows.

**Step 2, BAM → calls** (Slurm; submitted from `data_preparation/fiberseq/`):

```bash
W=$GROUP_SCRATCH/ndiamant/fiberseq/GM12878
BAM=$OAK/datasets/SMF_fiberseq/raw_bams/GM12878/GM12878.fire.bam REGIONS=$W/candidates.bed OUT=$W/calls \
  sbatch --export=ALL --array=20 -o $W/calls/logs/%A_%a.log fiberbam_to_calls.sbatch        # job 46646815, chr21 test
BAM=... REGIONS=... OUT=... \
  sbatch --export=ALL --array=0-19,21-22 -o $W/calls/logs/%A_%a.log fiberbam_to_calls.sbatch  # job 46652102
```

- **Resources:** 8 CPUs, 32 GB, 12 h, on partitions
  `owners,normal,akundaje,stat,hns,btrippe`. `btrippe` is a single GPU node and
  can never take these jobs, because Sherlock adds `NO_GPU` to CPU-only jobs;
  it can be dropped from the list.
- **Run:** every task completed on its first attempt. Tasks took 5–31 min;
  peak memory was about 12 GB.
- **Output:** 23 `<chrom>.tsv.zst` files, 37 GB, about 21.0M molecules.
- **chr21 check against the old TSV:**
  - Molecules present in both files agree on m6A rate (0.1044 new vs 0.1043 old).
  - Molecules new to this run have a rate of 0.079, as expected for closed
    background windows.
  - 883 old molecules (2%) are absent, presumably because they fell outside
    every new candidate window. This was not checked individually.
  - 91% of molecules carry a haplotype tag, split evenly.

**Step 3, calls → H5** (run directly in the interactive allocation, 16.5 min
on 8 CPUs, 1.6 GB peak):

```bash
python data_preparation/fiberseq/calls_to_h5.py --calls-dir $W/calls --candidates $W/candidates.bed \
    --fasta /scratch/groups/btrippe/brian/projects/SMF/fiberseq/reference/hg38.fa \
    --out-dir $W/h5 --workers 8
```

- **Storage rules:**
  - A molecule is stored in a window if it overlaps it by at least 500 bp.
  - A window is kept if at least 32 molecules span all of it.
  - Background is GC-matched only after the depth filter, so the 1:1 ratio
    holds for the windows actually written.
- **Groups per window:** `smf_mat` (int8: 1 accessible / 0 protected / −1 no
  call), `smf_pos`, `read_id`, `strand`, `haplotype`.
- **A/T check:** every kept window is loaded back through `load_dna` and
  checked with `get_valid_mask(..., "fiberseq_m6a")`. All windows passed.
- **Output:** `h5/active.h5` (84,125 windows, 2.9 GB), `h5/background.h5`
  (83,343 windows, 2.8 GB; 99.1% of a 1:1 match), and `h5/manifest.tsv`.
- **chr21 sanity check:**
  - The two classes have similar depth (median 132 vs 129 stored molecules,
    about 1,100 A/T columns per window) and similar GC (0.490 vs 0.486).
  - Active windows have the higher m6A rate (0.137 vs 0.110).
  - Loading with the assay filter drops 0 columns.

## Code changes (SMFNet, branch `nd/more_assays`, uncommitted)

- **`src/smf_net/data.py`:**
  - `fiberseq_m6a` assay with target rule `at`, RC on by default.
  - `get_valid_mask` returns every A and T.
  - `reverse_complement_dna_and_smf(dna, smf, assay_type)` flips the calls
    without a shift for `at`; `get_item` passes the assay type through.
- **`src/smf_net/evaluation/flank_baselines.py`:**
  - Separate `A` and `T` site models for `at`.
  - Unknown target rules now raise instead of falling through to the CpG/GpC
    branch.
- **Tests:**
  - New: `tests/test_fiberbam_to_calls.py` (synthetic BAM, hand-derived rows),
    `tests/test_select_regions.py`, `tests/test_calls_to_h5.py`.
  - Updated: `test_data.py`, `test_model.py` and
    `test_evaluation_flank_baselines.py` cover the new assay.
  - The full suite passed (359 passed, 3 skipped) before the H5 tests and the
    model-test case were added; those were then run on their own and pass.
- **Bug fixed in passing:** `calls_to_h5.py` first raised `SystemExit` inside a
  `multiprocessing.Pool` worker, which hung the parent, because Pool only
  passes `Exception` subclasses back. It now raises `ValueError`.

## Recommended training command (not yet run)

```bash
W=$GROUP_SCRATCH/ndiamant/fiberseq/GM12878/h5
python scripts/train.py \
  hydra.run.dir=/scratch/users/diamant/smf_models/GM12878_fiberseq \
  wandb_logger.save_dir=/scratch/users/diamant/wandb \
  wandb_logger.name=GM12878_fiberseq \
  assay_type=fiberseq_m6a \
  "data.cfg.h5_path=[$W/active.h5,$W/background.h5]" \
  data.cfg.fasta_path=/scratch/groups/btrippe/brian/projects/SMF/fiberseq/reference/hg38.fa \
  "data.cfg.train_chrm=[chr2,chr4,chr5,chr7,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr21,chr22,chrX,chrY]" \
  "data.cfg.val_chrm=[chr8,chr20]" \
  "data.cfg.test_chrm=[chr1,chr3,chr6]"
```

All windows are held in memory: about 24 GB of calls, about 17 GB of it in the
training split. A batch job should request at least 64 GB. Load time at this
size has not been measured. For a smoke test, add
`data.cfg.max_regions_per_h5_per_split=200 +trainer.max_steps=20` (the `+` may
need dropping) and use a separate `hydra.run.dir`.

# Related Notes

- [Flank baseline for profile metrics](flank-baseline-logistic-regression.md):
  - It documents the GCH strand asymmetry that made RC augmentation invalid for
    GCH. This note shows Fiber-seq lacks that asymmetry.
  - The flank baseline is also the evaluation-side answer to Hia5's strong
    context preference (0.011 at AAA to 0.22 at GAC) and to the leftover CpG
    o/e shortcut. This work added its A/T site types.
- [Positional covariation baselines](../2026-09-29/positional-covariation-baselines-and-noise-matched-corr-of-corr.md):
  the other train-set baselines a Fiber-seq model will be scored against.

# Open Questions

- **Preserve the outputs.** `$GROUP_SCRATCH` is purged after 90 days. The H5
  files, `candidates.bed`, the manifest and possibly the 37 GB of calls should
  be copied to `$OAK`.
- **Training is not run yet.** Load time and memory at 167k windows are
  untested.
- **Background is GC-matched genome-wide, not per split.** Test has 18.7k
  active against 16.4k background windows, val the reverse. This only matters
  if test-set active-vs-background comparisons become a headline number.
- **ML ≥ 230 threshold.** It was inherited from the `smf_hmm` pipeline. Methods
  should justify or cite it; it has not been checked against FIRE's own
  threshold.
- **The 883 old chr21 molecules missing from the new extraction** were not
  checked individually.
- **Second dataset.** Drosophila S2, with its naked-DNA and no-Hia5 controls,
  would calibrate Hia5 sequence bias directly. Other human samples (COLO829,
  Jurkat, THP1, A549) can be run through the same pipeline from their BAMs.
- **Remaining CpG o/e shortcut** (AUC 0.67): watch region-to-region metrics
  against the composition and flank baselines.

# Sources

- SMFNet: `data_preparation/fiberseq/README.md` (full step-by-step
  documentation, inputs table, strand-symmetry table) and the scripts beside it.
- smf_hmm: `evaluation/datasets/archive/fiberseq_survey/FIBERSEQ_DATA_REPORT.md`
  (dataset survey).
- smf_hmm: `evaluation/datasets/smf/README.md`, `SPEC.md`, `samples.tsv`
  (existing TSVs and their GM12878-anchored 10 kb windows).
- smf_hmm: `20260420_create_h5_from_mBED.py` (mouse ONT H5 layout).
- Analysis scripts in `/scratch/users/diamant/fiberseq_strand/`:
  `fiberseq_symmetry.py`, `depth.py`, `chrombpnet_regions.py`, `peak_depth.py`,
  `matching_diag.py`, `ccre_split.py`. GCH results are in
  `/scratch/users/diamant/gch_rc_analysis/sym_*.txt`.
- ENCODE: SCREEN v3 cCREs `ENCFF733BFV`; GM12878 ATAC `ENCFF614SMH` (IDR) and
  `ENCFF158ORB` (pseudoreplicated); hg38 blacklist v2,
  <https://github.com/Boyle-Lab/Blacklist>.
- Vollger et al. 2025, *A haplotype-resolved view of human gene regulation*
  (FIRE), `s3://stergachis-manuscript-data/2024/Vollger_et_al/FIRE/`.
