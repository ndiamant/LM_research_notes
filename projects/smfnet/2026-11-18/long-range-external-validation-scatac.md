---
title: "External validation of 16 kb SMFNet long-range coupling: scATAC cell-level co-accessibility of the CRE pairs (24k GM12878 cells), ±250 bp element windows, a raw-diffusion 50-epoch run, and a haplotype-imbalance pilot. All models tie on scATAC (Spearman 0.08-0.09, about 70% of Fiber-seq's own agreement); the no-b_m models' Fiber-seq lead does not transfer"
date: 2026-11-18
project: smfnet
agent: Claude Code
status: draft
sources:
  - SMFNet branch nd_scatac (dev remote), commit 9968206: src/smf_net/evaluation/scatac.py, scripts/scatac_pileup.py, scripts/scatac_pairs.py, tests/test_scatac.py, pixi environment scatac (bioconda AMULET 1.1)
  - 10x Genomics public human (GM12878) / mouse barnyard scATAC runs (cf.10xgenomics.com), CC BY 4.0
  - GIAB HG001 (NA12878) v4.2.1 GRCh38 benchmark variants (ftp-trace.ncbi.nlm.nih.gov)
  - Outputs under $SCRATCH/smfnet_figwig/long_range/ (scatac/, cre_pairs_val_hw250/, scores_val_hw250_*, *_raw16_e50*) and $SCRATCH/smfnet_figwig/haplotype/
  - Jobs 52590039 (model pair scores at ±250 bp), 52625990 (raw 50-epoch training), 52859311 (evaluation of the raw run)
  - User discussion, 2026-11-16 to 2026-11-18
tags: [smfnet, fiberseq, gm12878, evaluation, long-range, co-accessibility, scatac, amulet, haplotype, imprinting, labelling-efficiency, raw-diffusion]
---

# Summary

This note adds checks of the 16 kb models' long-range coupling from outside Fiber-seq. It builds on the
[long-range evaluation suite](../2026-11-16/long-range-evaluation-suite.md).

**Main findings:**

- **scATAC works as an external check of the CRE pair test.**
  - **Data:** 24,045 GM12878 cells from seven public 10x human/mouse runs.
  - **Statistic:** for each val pair, the correlation across cells of the two elements' Tn5 insertion counts, after regressing out depth within each run.
  - **Agreement with Fiber-seq:** it agrees with the Fiber-seq partial r (Spearman 0.121 at ±250 bp, about 0.5-0.7 after correcting for both measures' reliabilities). It does not agree with the raw r, which carries labelling.
- **±250 bp element windows double the reliability of both measures.** Fiber-seq split-half reliability goes from 0.10 to 0.20 and scATAC's from 0.13 to 0.27 (ceilings 0.45 and 0.52). They should be the suite's default.
- **Models against scATAC.** Every 16 kb model reaches Spearman 0.081-0.091 against scATAC at matched ±250 bp windows. That is about 70% of the independent Fiber-seq data's own agreement with scATAC (0.121). No model differs from another by more than ±0.01, with intervals of about ±0.03.
- **The no-b_m lead is specific to Fiber-seq.**
  - **On Fiber-seq:** the two models without b_m (the pilot, and the new raw-diffusion jitter + short-batch + background run) lead on the Fiber-seq partial r: 0.149-0.150 against 0.084-0.103 for the additive models. The raw run leads already at epoch 24.
  - **On scATAC:** neither leads.
  - **Interpretation (inferred, not directly tested):** a model without b_m has to infer each molecule's labelling from its calls, and the linear labelling adjustment in the partial r leaves some labelling in. That would explain why the [suite's open question](../2026-11-16/long-range-evaluation-suite.md#open-questions) about the pilot's lead does not reflect better long-range biology.
- **The raw-diffusion run trained for 50 epochs** (no emission head, no b_m, same data recipe as the additive background run):
  - **Val loss:** 0.15668, equal to the additive run at 25 epochs (0.15665).
  - **Short range:** +0.96% against the 2 kb additive model, against +0.26% for the additive run and +0.84% for the original raw 2 kb model.
  - **Labelling spread:** 0.64 of the observed, as the pilot.
  - **Position:** no position bias.
- **Haplotype pilot (chr20 and chrX).**
  - **The test works:** the haplotype-label null is calibrated, chrX's inactive X shows in 36% of elements, and the imprinted GNAS control appears.
  - **Autosomal imbalance is rare:** 0.5% of chr20 elements at q < 0.05, and few have no heterozygous variant within 1 kb.
  - **The user chose scATAC first;** the genome-wide haplotype run was not done.

# Key Points

- **What each assay measures.**
  - **Estimand:** with a cell-level part `c_cell` (shared cell state) and a molecule-level part `c_mol`:
    - Fiber-seq measures `c_cell + c_mol`.
    - Unphased diploid scATAC measures `4·c_cell + 2·c_mol`, which gives the cell-level part twice the weight. The extra term comes from copies on the two homologs.
  - **No single-molecule signal at this range:** scATAC fragments are under about 1 kb, so for pairs 2 kb or more apart no fragment links the two elements. Phasing scATAC is impractical: about 8% of fragments carry a heterozygous SNP, so about 0.6% of co-detections are phased.
  - **Hi-C and Micro-C measure contact, not co-accessibility:**
    - **Few loop calls at this range:** loop calls barely exist below 20-40 kb.
    - **Small expected effect even with true loops:** a loop present 10% of the time that raises both elements' openness by 0.3 gives ρ ≈ 0.03.
    - **Use:** only as a stratifier (insulation, short-range dots).
  - **Data availability:** 4DN has no GM12878 Micro-C. ENCODE has GM12878 intact Hi-C with MNase (ENCSR916MFV, 1.87B contacts).
- **scATAC data and cell calling.**
  - **GRCh38 runs:** `10k_hgmm_ATACv2_nextgem_Chromium_Controller`, `10k_hgmm_ATACv2_nextgem_Chromium_X`, `10k_hgmm_ATACv1p1_nextgem_Chromium_X` and `atac_hgmm_10k_nextgem`. After dropping mixed-species barcodes: 18,585 human-only cells, a median of about 23k fragments each.
  - **hg19 runs:** `atac_v1_hgmm_10k`, `_1k` and `_500` add 5,460 cells, counted at element centres lifted to hg19. 116 of 84,241 elements do not lift; their pairs are skipped.
  - **AMULET multiplets:** it flags 563-734 cells per run (14%), about 3× the roughly 5% human-human doublets implied by the mixed-species rate. Flagged cells have 1.7-2.0× the median depth.
    - Removing them slightly lowers every metric, which fits a doublet being the sum of two independent cells: its expected correlation is unchanged and only its depth changes. All cells are used.
  - **Capture:** about 0.06 insertions per element per cell at ±100 bp (0.12 at ±250 bp), about a third of the napkin estimate. scATAC is therefore no more powerful per pair than Fiber-seq.
- **scATAC pair statistics** (val, 3,502 pairs from the ±100 bp Fiber-seq table, all human-only cells):

  | configuration | cells | split-half reliability | ceiling | p < 0.01 and r > 0 | Spearman with Fiber-seq partial r |
  | --- | --- | --- | --- | --- | --- |
  | ±100 bp, 4 runs | 18,585 | 0.127 | 0.36 | 6.3% | 0.078 |
  | ±100 bp, 4 runs, AMULET multiplets removed | 15,980 | 0.121 | 0.35 | 5.7% | 0.063 |
  | ±100 bp, 7 runs | 24,045 | 0.152 | 0.39 | 7.2% | 0.088 |
  | ±250 bp, 4 runs | 18,585 | 0.206 | 0.45 | 10.3% | 0.091 |
  | ±250 bp, 7 runs | 24,045 | 0.266 | 0.52 | 11.9% | 0.092 |

  - **Calibration:** p < 0.01 with r < 0 is 0.1% of pairs, against 6.3% with r > 0.
  - **Raw r:** the Spearman with the Fiber-seq raw r is 0.008 (n.s.).
- **Fiber-seq pair table at ±250 bp** (`cre_pairs_val_hw250`):
  - 3,378 val pairs; molecules must now cover each element ±250 bp.
  - 152 positives (4.5%).
  - Split-half reliability 0.199, ceiling 0.45.
  - Spearman with scATAC (±250 bp, 7 runs, 3,366 pairs): 0.121 [0.087, 0.152].
- **Models at ±250 bp throughout** (Spearman with 95% pair-bootstrap intervals; noise masking):

  | model | vs scATAC r | vs Fiber-seq partial r | AUROC / AUPRC on Fiber-seq positives |
  | --- | --- | --- | --- |
  | raw, 50 epochs (epoch 49) | +0.091 [+0.055, +0.128] | +0.150 [+0.117, +0.183] | 0.703 / 0.157 |
  | raw, 50-epoch schedule, epoch 24 | +0.082 [+0.048, +0.117] | +0.150 [+0.117, +0.186] | 0.678 / 0.154 |
  | 16 kb pilot | +0.081 [+0.046, +0.116] | +0.150 [+0.117, +0.186] | 0.700 / 0.143 |
  | 16 kb additive | +0.084 [+0.050, +0.121] | +0.086 [+0.052, +0.122] | 0.644 / 0.128 |
  | + jitter + short batches | +0.089 [+0.054, +0.122] | +0.088 [+0.053, +0.120] | 0.666 / 0.126 |
  | + 10% background | +0.085 [+0.051, +0.120] | +0.103 [+0.068, +0.135] | 0.661 / 0.134 |
  | observed Fiber-seq partial r | +0.121 [+0.087, +0.152] | | |

  - **Paired differences on scATAC** (each model minus the raw 50-epoch model): −0.009 to −0.002, with intervals within about ±0.03.
  - **Beyond observed Fiber-seq:** the partial Spearman with scATAC given the Fiber-seq partial r is 0.064-0.079. Part of what the models share with scATAC is not explained by the noisy observed Fiber-seq.
  - **Ranking scATAC positives** (404 pairs with p < 0.01 and r > 0): AUROC 0.56-0.57 for every model, 0.60 for the Fiber-seq partial r.
- **Window mismatch matters.**
  - **Mismatched windows:** with model scores at ±100 bp against scATAC at ±250 bp, model agreement was lower (0.031-0.057) even though scATAC was more reliable.
  - **Matched windows:** with ±100 bp models against ±100 bp scATAC over 7 runs, it is 0.044-0.072.
  - **Conclusion:** model scores, Fiber-seq and scATAC should all use the same window.
- **Fiber-seq pair metrics, both windows:**

  | model | ±100 bp: Spearman / given raw r / AUROC / AUPRC | ±250 bp: same |
  | --- | --- | --- |
  | 16 kb pilot | 0.138 / 0.114 / 0.635 / 0.108 | 0.149 / 0.113 / 0.700 / 0.143 |
  | 16 kb additive | 0.071 / 0.080 / 0.643 / 0.099 | 0.084 / 0.062 / 0.644 / 0.128 |
  | + jitter + short batches | 0.082 / 0.065 / 0.629 / 0.090 | 0.087 / 0.070 / 0.666 / 0.126 |
  | + 10% background | 0.066 / 0.077 / 0.620 / 0.094 | 0.102 / 0.071 / 0.661 / 0.134 |
  | raw, 50 epochs | 0.133 / 0.106 / 0.649 / 0.080 | 0.149 / 0.104 / 0.703 / 0.157 |
  | raw, epoch 24 | 0.120 / 0.092 / 0.641 / 0.101 | 0.148 / 0.111 / 0.678 / 0.154 |

- **The raw 50-epoch run** (job 52625990, 22.1 h on one H100):
  - **Val loss:**

    | epoch | 4 | 9 | 14 | 19 | 24 | 29 | 34 | 39 | 44 | 49 |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | val/loss | 0.15879 | 0.15801 | 0.15772 | 0.15753 | 0.15735 | 0.15719 | 0.15702 | 0.15685 | 0.15673 | 0.15668 |

    - Still falling slightly at the end.
    - **Comparisons:** the additive background run reached 0.15665 at its epoch 24 and the pilot 0.15686.
    - **8 kb short-window val loss:** 0.16037, against the additive run's 0.16024.
  - **Short range** (central 1 kb of the 12,304 2 kb val regions; loss per call compared with the 2 kb additive model; intervals about ±0.02-0.05%):

    | model | all sites | CRE sites | other sites |
    | --- | --- | --- | --- |
    | 16 kb pilot | +4.16% | +3.86% | +4.23% |
    | + 10% background (additive) | +0.26% | +0.35% | +0.24% |
    | raw, epoch 49 | +0.96% | +0.87% | +0.99% |
    | raw, epoch 24 | +1.39% | +1.35% | +1.40% |

    The original raw 2 kb model scored +0.84% in the [suite note](../2026-11-16/long-range-evaluation-suite.md).
  - **Labelling spread:**
    - **Spread:** generated reads vary 0.64 as much as observed reads (0.68 at epoch 24), and r(read mean, b_m) = −0.01. Both are as for the pilot (0.63, −0.01).
    - **Raw co-accessibility at 0.5-1 kb:** 0.087, against observed 0.162. The raw model undershoots, while additive models with matched b_m overshoot (0.195).
    - **Adjusted for the read mean:** about 0.
  - **Position:** the largest relative loss change across offsets up to ±4 kb is 0.014% (pilot 0.07%, additive background run 0.11%).
- **The standard val sample metrics have not been run** (`eval_16kb_val.sh`: 1,000 val regions × 512 generated reads, then `metrics.py`). Of the 16 kb models only the pilot has them. Its sampling took about 35 min on an H100, and the metrics step about 8 h.

# Details

## scATAC method

`smf_net.evaluation.scatac` (commit 9968206) does three things:

1. **Cell calling** (`read_cells`): keeps barcodes called cells of one species only, from 10x `singlecell.csv`. Covariates are passing fragments (linear and log), fraction in peaks, and fraction at TSSs.
2. **Pileup** (`element_cut_counts`): streams the fragment file once and counts both Tn5 insertions of each fragment (`start`, `end - 1`; 10x fragments are already shifted) within `half_width` of each active cCRE centre (84,241 elements genome-wide). This takes about 3 min per run.
   - **Lifted runs:** `remap_columns` returns counts made at lifted centres to the hg38 element columns and marks unlifted ones as missing.
3. **Pair statistics** (`pair_correlations`):
   - **Depth within each run:** depth covariates are regressed out of each element's counts within each run, so chemistry and depth differences don't create co-accessibility.
   - **Pooling:** residuals are pooled over runs.
   - **Statistic:** Pearson r per pair, with a Fisher-z p-value (dof = cells − 2 − runs × (covariates + 1)). r is NaN when an element has fewer than 20 insertions in total.
   - **Split halves:** r is also computed on two random halves of the cells.
   - **Model comparison:** `model_agreement` adds the model comparison and its bootstrap.

**Napkin power math versus what was measured.** scATAC's observed correlation is roughly the true correlation times √(d_A·d_B·(1−p_A)(1−p_B)), where d is the insertions per open copy and p the open probability, and its standard error is 1/√N_cells.
- **Predicted:** with d ≈ 0.2, about 5k-10k cells would match the current Fiber-seq per pair.
- **Measured:** d was about 0.06-0.07 at ±100 bp. 24k cells give a reliability (0.15) only slightly above Fiber-seq's (0.10), which agrees with the napkin math once the measured d is used.

## Why the no-b_m lead looks like labelling (inferred)

The pattern across the six models:

| models | Fiber-seq partial r | scATAC |
| --- | --- | --- |
| pilot, raw run (no b_m) | 0.15 | 0.08-0.09 |
| additive runs (with b_m) | 0.08-0.10 | 0.08-0.09 |

- **Not training length:** the raw run's lead is present at epoch 24 of its 50-epoch schedule.
- **Not the data recipe:** the raw and additive background runs share data and training settings.
- **Proposed mechanism (inferred):** a model without b_m treats the partner element's calls as evidence about the molecule's labelling. The partial r removes labelling only linearly (level and b_m), so pairs whose co-variation is partly non-linear labelling get high partial r. A no-b_m model ranks those pairs highly; scATAC has no m6A labelling and does not.
- **Not yet tested:** a non-linear labelling adjustment of the partial r, or scoring an additive model at a fixed b_m.

## scATAC datasets for pooling

The 10x mixing runs pooled here are one source among several. From a subagent inventory on 2026-11-16, ranked by a rough 10x-cell equivalent. Per-cell power scales with about the square of in-peak insertions; the 10x runs have about 13k in-peak fragments per cell. These estimates come from the reported depths, not measurements.

| dataset | GM12878 cells | estimated equivalent | effort |
| --- | --- | --- | --- |
| dscATAC (Lareau 2019, GSE123581) | ~5.9k at 55-79k fragments | ~25-30k | High: no fragment file; realign and redo barcode merging |
| METATAC (Wu 2022, PRJNA789047) | 1.1k at 142k fragments | ~25-45k | High: raw reads only, plate barcodes |
| SHARE-seq rep3 (Ma 2020, GSE140203) | ~67k barcodes at ~8k fragments | ~10k | Moderate: hg19 fragments; call cells ourselves; barcode collisions |
| ENCODE4 ENCSR680NPV + ENCSR371ZNU | ~5.5k + ~23k | ~4k + ~3k | Moderate: hg38 fragments; our own cell calling; high multiplet flags |
| BABEL multiome (Wu 2021, GSE166797) | ~3.5-7k (unverified) | ~2-5k | Low to moderate |

- **Code needed for non-10x sets:** cell calling and depth covariates from fragments alone.
- **Skip:** the perturbed sets (Spear-ATAC, Perturb-ATAC) and the very sparse ones.

## Haplotype pilot

**Question:** are there enough strong haplotype-resolved accessibility differences in the Fiber-seq data to test the model's long-range sequence effects? The proposed test predicts B's haplotype difference from haplotype-specific DNA, with variants only near A.

**Method** (`scripts/haplotype_imbalance.py`):
- **Data:** molecules already carry an HP haplotype tag in `<chrom>.molecules.npz`. 89% of chr20 molecules are phased; the mean span is 17.1 kb.
- **Test:** for each active cCRE, the element's accessible fraction (±100 bp) is regressed by OLS on a haplotype indicator plus b_m, or plus b_m and the labelling level.
- **Null:** the same test with haplotype labels shuffled within the element.
- **Variants:** distance to the nearest GIAB HG001 heterozygous variant, using only elements with ≥ 90% confident coverage over ±1 kb.

**Results** (chr20, a val chromosome, and chrX; about 50 molecules per haplotype):
- **Calibration:** the shuffled-label p-values are uniform (KS p = 0.11), and 0 null elements reach q < 0.05.
- **chrX:** 670 of 1,844 elements at q < 0.05; about half are estimated to be imbalanced.
- **chr20:** 12 of 2,197 elements at q < 0.05; 3-5% are estimated to be imbalanced.
  - **Effect sizes:** the hits differ by a median of 0.11 in accessible fraction; the per-element detection limit is about 0.076.
  - **Imprinted GNAS element** (chr20:58,850,945): 0.15 against 0.54 accessible.
  - **Nearest heterozygous variant among the hits:**
    - 2 have one inside ±100 bp.
    - 4 have none closer than 1.3-2.0 kb; only the 16 kb model could see those.
    - 1 is pericentromeric with no confident calls, likely an artefact.
    - 0 of 599 elements whose nearest variant is over 2 kb away (near-identical haplotypes, fewer phased reads).
- **Reliability** of the per-element haplotype difference: 0.23 overall, 0.13 for elements with no variant within 1 kb.
- **Extrapolation:** about 400 strong imbalances among 82,219 active autosomal elements, perhaps 100-150 without a variant within 1 kb. Enough to pool over, not for per-element claims.
- **Caveats:**
  - Imprinting and the inactive X are parent-of-origin effects that the DNA-only model cannot predict, so they must be excluded.
  - GIAB lacks structural variants.
  - HP labels are consistent only within a phase block, so directional tests need read-level assignment from `GM12878.fire.bam`.

## Commands

The environment is as in the [suite note](../2026-11-16/long-range-evaluation-suite.md#commands): pixi with `--frozen --manifest-path ~/repos/SMFNet_unet_figwig/pyproject.toml`, `PYTHONPATH=<checkout>/src`, and `CONDA_OVERRIDE_CUDA=12.0` on CPU nodes. CPU steps ran in the interactive allocations 52445748 and 52719102.

```bash
X=$SCRATCH/scatac/GM12878/10x; L=$SCRATCH/smfnet_figwig/long_range; S=$L/scatac
B=https://cf.10xgenomics.com/samples/cell-atac
# downloads: <version>/<run>/<run>_{fragments.tsv.gz,fragments.tsv.gz.tbi,singlecell.csv,summary.csv}
#   2.1.0: 10k_hgmm_ATACv2_nextgem_Chromium_Controller, 10k_hgmm_ATACv2_nextgem_Chromium_X,
#          10k_hgmm_ATACv1p1_nextgem_Chromium_X; 2.0.0: atac_hgmm_10k_nextgem;
#   1.2.0 (hg19): atac_v1_hgmm_10k, atac_v1_hgmm_1k, atac_v1_hgmm_500
# AMULET (pixi -e scatac): barcode CSV with is__cell_barcode = human-only cells; chromosome list
# and repeat filter (AMULET v1.0 release RestrictionRepeatLists, hg38) prefixed with GRCh38_
pixi run -e scatac amulet $X/<run>_fragments.tsv.gz amulet/<run>_human_cells.csv \
  amulet/autosomes_GRCh38_prefixed.txt amulet/repeats_GRCh38_prefixed.bed amulet/<run>/
# element centres hg38 -> hg19 (ml biology ucsc-utils/v489; UCSC hg38ToHg19.over.chain.gz), kept if
# same chromosome, written as liftover_active_hg38_to_hg19.tsv (chrom, centre, lifted_chrom, lifted_centre)
liftOver centres_hg38.bed hg38ToHg19.over.chain.gz centres_hg19.bed centres_unlifted.txt
# pileups (GRCh38 runs; hg19 runs add species=hg19 contig_prefix=hg19_ liftover=$X/liftover_active_hg38_to_hg19.tsv)
pixi run python scripts/scatac_pileup.py hydra.run.dir=$S/pileup_hw250/<run> half_width=250 \
  fragments=$X/<run>_fragments.tsv.gz singlecell=$X/<run>_singlecell.csv \
  elements_bed=/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/candidates.bed
# Fiber-seq pair table at +/- 250 bp
pixi run python scripts/build_cre_pairs.py hydra.run.dir=$L/cre_pairs_val_hw250 figwig_dir=$F/figwig \
  efficiency_dir=$F/efficiency elements_bed=.../candidates.bed "chromosomes=[chr8,chr20]" pairs.half_width=250
# model pair scores at +/- 250 bp: scripts/score_pairs_hw250.sbatch (job 52590039), eval_raw16_e50.sbatch (job 52859311)
# scATAC pair statistics and model agreement
pixi run python scripts/scatac_pairs.py hydra.run.dir=$S/pairs_fs250_pileup_hw250_7runs_raw50 \
  pairs_path=$L/cre_pairs_val_hw250/cre_pairs.tsv "pileups=[$S/pileup_hw250/<run>,...]" \
  "model_scores=[raw16_e50:$L/scores_val_hw250_raw16_e50,pilot16:$L/scores_val_hw250_pilot16,...]"
# haplotype pilot
pixi run python haplotype_imbalance.py pilot_chr20_chrX.tsv chr20 chrX; pixi run python summarize.py pilot_chr20_chrX.tsv
```

**Scripts and runs:**
- **Scripts:** copied to [long-range-external-validation-scatac/scripts/](long-range-external-validation-scatac/scripts/) (the haplotype pilot, the two evaluation `.sbatch` files, and the raw run's training `.sh`/`.sbatch`).
- **Raw run:**
  - checkpoint: `/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_16kb_flank1024_jitter_mixed_bg10_e50_alphagenome_unet_attn_muon_ema_k9_dilated/epoch49-step182550.ckpt`, periodic checkpoints every 5 epochs in `periodic/`;
  - code: trained from the detached checkout `~/repos/SMFNet_scatac_runs` at 52518fb, whose training code is identical to 759ba00, the additive runs'.
- **Purge risk:** the per-pair tables (`scatac_pairs.tsv`, 1.2 MB; `cre_pairs_val_hw250/cre_pairs.tsv`, 0.9 MB) stay on `$SCRATCH`. The commands above regenerate them from the public 10x files and the Fiber-seq inputs.

# Related Notes

- [Long-range evaluation suite](../2026-11-16/long-range-evaluation-suite.md): The pair test, ceilings and models this note validates externally. Its open question on the pilot's Spearman lead is answered here (inferred: labelling, not biology).
- [Additive emission head](../2026-11-15/fiberseq-additive-emission-head.md): b_m and the additive models whose Fiber-seq scores trail the no-b_m models.
- [16 kb pilot and labelling efficiency](../2026-11-12/fiberseq-16kb-pilot-labelling-efficiency.md): The pilot and the origin of the CRE-pair analysis.

# Open Questions

- **Non-linear labelling:** does a non-linear labelling adjustment of the Fiber-seq partial r (e.g. b_m and level splines, or matching on them) remove the no-b_m models' lead? Does an additive model scored at fixed b_m gain it?
- **Making ±250 bp the default:** should the suite switch, and should the window change the checkpoint-scan and ceiling numbers?
- **More scATAC:** pooling SHARE-seq, ENCODE and BABEL (fragment files exist) might add about 20k cell equivalents; dscATAC or METATAC (reprocessing reads) maybe 25-45k each. Reliability at ±250 bp could reach about 0.4.
- **Cell-level versus molecule-level:** about half the true pair-to-pair variation is shared between scATAC and Fiber-seq. How much of the rest is cell-level co-variation, which scATAC weights twice?
- **Standard val metrics:** these have not been run for the 16 kb additive runs or the raw 50-epoch run (sampling about 35 min, metrics about 8 h per model).
- **Haplotype run:** the genome-wide haplotype-imbalance run is on hold. It would need imprinted regions excluded, ±250 bp windows, structural-variant calls and read-level haplotype assignment.
- **Test set:** test chromosomes stay unused; the user does not want to spend the test set. Larger evaluation sets should come from more val data or more external data.

# Sources

- **SMFNet code:** `nd_scatac` at 9968206 (dev remote): `src/smf_net/evaluation/scatac.py`, `scripts/scatac_pileup.py`, `scripts/scatac_pairs.py`, `tests/test_scatac.py`, `pyproject.toml` (scatac pixi environment).
- **10x Genomics public datasets** (cf.10xgenomics.com/samples/cell-atac/...), CC BY 4.0.
- **AMULET:** Thibodeau et al. 2021, Genome Biology, <https://doi.org/10.1186/s13059-021-02469-x>; <https://github.com/UcarLab/AMULET>; bioconda `amulet` 1.1.
- **GIAB HG001 v4.2.1:** <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/>.
- **UCSC liftOver chain:** <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz>.
- **Outputs:** `$SCRATCH/smfnet_figwig/long_range/` and `$SCRATCH/smfnet_figwig/haplotype/`; downloads in `$SCRATCH/scatac/GM12878/10x/` and `$SCRATCH/fiberseq/GM12878/variants/`.
