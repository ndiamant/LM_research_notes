---
title: FootprintCharter baseline integration and missingness limitations
date: 2026-08-18
project: smf-hmm
agent: Codex
status: complete
sources:
  - https://github.com/ndiamant/smf_hmm/tree/79c44d4/baselines/footprint_charter
  - https://bioconductor.org/packages/3.22/bioc/vignettes/SingleMoleculeFootprinting/inst/doc/FootprintCharter.html
  - https://git.embl.org/grp-krebs/nf-smfont
  - https://nanoporetech.github.io/modkit/filtering.html
  - /scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5
tags:
  - smf-hmm
  - footprintcharter
  - nanopore
  - baseline
  - missing-data
  - sherlock
---

# Summary

A reproducible FootprintCharter 2.4.0 baseline was installed in an Apptainer sandbox on Stanford Sherlock and connected directly to the `smf_hmm` H5 format. The integration successfully reads H5 matrices, runs FootprintCharter, and writes its native footprint tables and three plotting-utility outputs.

The resulting state calls are not currently a useful baseline. FootprintCharter's strict complete-case filtering leaves only 30-32 molecules in the three most compatible regions found in a scan of 77,030 H5 regions. It reduces the requested 16 partitions to 3-4, and the resulting calls include implausibly long nucleosome segments and unstable short TF/accessibility segments. The user's visual assessment is that the state calls are quite poor; low retained molecule count is a leading explanation, although sparse jointly observed positions and interpolation likely also contribute.

The principal next direction is an opt-in missing-aware FootprintCharter wrapper: retain partially observed reads and positions, compute pairwise NaN-aware Euclidean distances with overlap normalization and a minimum-overlap requirement, then reuse FootprintCharter's clustering output structure, footprint detection, and plotting functions.

# Key Points

- Sherlock's available R modules were not a clean match for the required Bioconductor release, so the baseline uses `bioconductor/bioconductor_docker:RELEASE_3_22` through Apptainer.
- The installed stack pins `SingleMoleculeFootprinting` 2.4.0, `plyranges` 1.30.0, `stringfish` 0.16.0, and `qs` 0.27.3; `rhdf5` was added for direct H5 reading.
- Apptainer 1.5.2 completed package installation and container tests but crashed in `mksquashfs` with exit 139. Building an uncompressed sandbox instead of a SIF avoids that final compression failure.
- H5 calls map to FootprintCharter sparse values as `-1 -> structural zero/NA`, `0 accessible -> 2`, and `1 protected -> 1`.
- `rhdf5` reverses the dimensions of these h5py-authored matrices, so the R adapter transposes them back to reads by positions.
- The H5 coordinates use a 1-based start and exclusive end, matching `load_dna(..., sub_one_from_end=True)`. A stored offset maps to genomic coordinate `start + offset`.
- All assayed columns in the inspected examples are CpG, GpC, or GCG sites. DNMT-TKO mESCs allow both CpG and GpC methylation to report accessibility.
- FootprintCharter first requires every retained molecule to have every assayed site in the padded 80 bp ROI, then requires every retained position to be observed in every retained molecule across the supplied matrix.
- The compatibility scan found 107 regions meeting the default criteria: at least 30 complete ROI reads, at least 250 bp of complete-case span, no retained-site gap over 60 bp, and unique read IDs.
- Even the three highest-ranked valid regions retain only 30-32 molecules. FootprintCharter reduced `k=16` to `k=3`, `4`, and `3`, respectively.
- Across 1,000 random H5 regions, 32.2% of matrix cells were missing. Approximately 24.9 percentage points were leading/trailing read-boundary gaps and 7.4 points were internal gaps. The median within-span internal missing rate was 8.5%.
- nf-smfONT's current default `llr=2` converts modification probabilities between approximately 0.119 and 0.881 to an ambiguous `x`, subsequently represented as `-1`. This likely explains a substantial part of the internal missingness, while most total missing cells arise because reads do not span the full 2,114 bp matrix region.
- A pasted shell continuation prompt (`>`) accidentally truncated the R runner, region list, and H5. The repository files were restored, and the H5 was restored from its verified `.zst` copy and made read-only.

# Details

## Durable implementation

The baseline implementation is committed in `smf_hmm` at commit `79c44d4` under `baselines/footprint_charter/`:

- `FootprintCharter.def`: reproducible Bioconductor image and pinned package stack.
- `build_container.sh`: Sherlock compute-allocation build into a sandbox directory.
- `run_footprint_charter.R`: direct H5 adapter, caller, result writer, and plotting entry point.
- `find_compatible_regions.py`: complete-case compatibility scanner.
- `example_regions.txt`: the three highest-ranked valid regions.
- `README.md`: build, execution, and scan commands.

The container sandbox is:

```text
/home/groups/btrippe/diamant/containers/footprint-charter_bioc-3.22.sandbox
```

The input data and reference are:

```text
/scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5
/scratch/users/diamant/smf_data/mm10.fa
```

The FASTA was used to validate coordinate and CpG/GpC context conventions but is not required by the current motif-unannotated FootprintCharter run.

## Installation history

The initial package build exposed three compatibility problems:

1. `qs` had been removed from current CRAN; archived `qs` 0.27.3 required compatible archived `stringfish` 0.16.0.
2. `SingleMoleculeFootprinting` 2.4.0 imports `filter` from `plyranges`, but `plyranges` 1.30.1 removed that export; version 1.30.0 works.
3. After successful installation and `%test`, Sherlock's bundled `mksquashfs` crashed while creating a SIF. A directory-format sandbox works with `apptainer exec`, `shell`, and `run` without this compression step.

The definition verifies archived source checksums and tests the pinned versions plus the exported `FootprintCharter` function.

## FootprintCharter filtering behavior

FootprintCharter converts sparse structural zeros to `NA`, stored 1 to protected/unmethylated 0, and stored 2 to accessible/methylated 1. Its internal `filter.dense.matrix` then:

1. Expands the 80 bp region of interest by 15 bp on each side.
2. Retains only reads with no missing calls at any represented site in that 110 bp interval.
3. Across the full supplied matrix, retains only cytosines with no missing call in any retained read.
4. Again removes reads with any missing value among the remaining columns.
5. Calculates a 40 bp rolling average and fills at most 20 consecutive empty rolling positions.

A retained-site gap over approximately 60 bp therefore leaves `NA` values in the smoothed matrix. `parallelDist::parDist` propagates those NAs and `cluster::pam` fails because it requires a finite dissimilarity matrix.

For the original HMM comparison regions:

| region | complete reads | complete sites | complete span | maximum gap | outcome |
| --- | ---: | ---: | ---: | ---: | --- |
| `chr16:96818819-96820933` | 45 | 11 | 198 bp | 115 bp | fails smoothing |
| `chr8:15107871-15109985` | 97 | 8 | 51 bp | 19 bp | technically clusterable but too short for a nucleosome call |
| `chr16:20764490-20766604` | 39 | 21 | 609 bp | 164 bp | fails smoothing |

Thus raw read count alone is not sufficient; jointly observed span and gap structure determine compatibility.

## Compatibility scan and selected examples

The scanner evaluated all 77,030 region matrices and excluded regions with duplicate read UUIDs. The default ranking first maximizes complete-case span, then retained reads and retained sites. The selected examples are:

| region | total reads | complete reads | ROI sites | complete sites | span | maximum gap | resulting partitions |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `chr8:87810742-87812856` | 188 | 30 | 20 | 30 | 431 bp | 51 bp | 3 (`15/7/8`) |
| `chr9:102834133-102836247` | 197 | 32 | 28 | 40 | 427 bp | 59 bp | 4 (`5/6/16/5`) |
| `chr8:119393978-119396092` | 188 | 30 | 30 | 42 | 393 bp | 54 bp | 3 (`13/7/10`) |

All three complete end-to-end and produce `footprints.tsv`, `results.rds`, `footprints.pdf`, `partitioned_molecules.pdf`, and `state_calls.pdf` under `baselines/footprint_charter/results/`. Those result files are currently untracked in `smf_hmm`.

The calls should not yet be treated as a quantitative baseline. Examples include inferred nucleosome intervals hundreds of bases long (including 701 and 735 bp segments) and TF calls supported by only one or a few retained cytosines. These widths are technically allowed by the configured thresholds, but the visual output is poor and partition sizes are small.

## Missingness provenance

The H5 contains the union of CpG-only, GpC-only, and GCG sites. In the three initially inspected regions, every matrix column belonged to at least one of these contexts.

An analysis of 1,000 randomly selected region matrices separated missing calls into leading/trailing gaps versus gaps between the first and last observed call of each read:

| quantity | fraction of all cells |
| --- | ---: |
| total missing | 32.2% |
| approximate read-boundary missing | 24.9% |
| internal missing | 7.4% |

Internal gaps were about 22.9% of all missing cells, and the median read had 8.5% internal missingness within its observed span. This classification is approximate and likely underestimates ambiguous calls near read boundaries.

The checked-out nf-smfONT source provides a likely mechanism for internal gaps. `BAM_to_mBED.py` converts BAM modification probabilities to log-likelihood ratios and `SMFONT_utils.py` emits `m` or `u` only when `abs(LLR) > 2`; intermediate values become `x`. Downstream scripts map `x` to `-1`. The converter also uses `exclude_zero_calls=True`, skipping positions without an emitted modification entry.

Even modest internal missingness is catastrophic under a complete-case-across-reads rule. If each call were observed independently with probability 0.92, a position would be complete across 45 reads with probability only `0.92^45`, approximately 2.3%.

## Reproduction commands

Build from a Sherlock compute allocation:

```bash
sh_dev -c 4
cd /home/users/diamant/repos/smf_hmm
./baselines/footprint_charter/build_container.sh
```

If using the sandbox built before `rhdf5` was added to the definition:

```bash
CONTAINER="$GROUP_HOME/$USER/containers/footprint-charter_bioc-3.22.sandbox"
apptainer exec --writable "$CONTAINER" Rscript -e \
    'BiocManager::install("rhdf5", ask = FALSE, update = FALSE, Ncpus = 4L)'
```

Scan compatible regions:

```bash
cd /home/users/diamant/repos/smf_hmm
/home/groups/btrippe/diamant/miniforge/envs/hmm/bin/python \
    baselines/footprint_charter/find_compatible_regions.py \
    /scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5 \
    --output baselines/footprint_charter/compatible_regions.tsv
```

Run the selected examples. The leading `>` characters shown by Bash as continuation prompts must not be copied because Bash interprets them as truncating output redirections.

```bash
CONTAINER="$GROUP_HOME/$USER/containers/footprint-charter_bioc-3.22.sandbox"

apptainer exec "$CONTAINER" Rscript \
    baselines/footprint_charter/run_footprint_charter.R \
    /scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5 \
    baselines/footprint_charter/example_regions.txt \
    baselines/footprint_charter/results
```

The accidentally truncated H5 was recovered from the adjacent tested archive and then protected against accidental writes:

```bash
zstd -d \
    /scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5.zst \
    -o /scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.recovered.h5

mv /scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.recovered.h5 \
   /scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5

chmod a-w /scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5
```

# Related Notes

- No earlier `smf-hmm` research note was found. The [project overview](../README.md) indexes this initial baseline note.

# Open Questions

- Can a missing-aware FootprintCharter wrapper produce stable partitions on the original HMM comparison regions?
- What minimum observed fraction per read and minimum pairwise overlap preserve useful coverage without clustering molecules by read boundaries?
- Should NaN-aware Euclidean distance be calculated only over the central 80 bp ROI or over the surrounding 500 bp calling window?
- How closely do the H5 `-1` values correspond to intermediate BAM modification probabilities versus absent MM/ML entries?
- Are state calls stable under read subsampling or bootstrap resampling? With only 30-32 retained molecules, instability is expected but has not yet been quantified.
- Is FootprintCharter still a meaningful baseline for this ONT representation if its complete-case assumption must be substantially changed?

# Sources

- [FootprintCharter baseline implementation in `smf_hmm`](https://github.com/ndiamant/smf_hmm/tree/79c44d4/baselines/footprint_charter): Container definition, build script, H5 adapter, compatible-region scanner, example list, and usage documentation.
- [FootprintCharter Bioconductor vignette](https://bioconductor.org/packages/3.22/bioc/vignettes/SingleMoleculeFootprinting/inst/doc/FootprintCharter.html): Package inputs, default parameters, result objects, and plotting utilities.
- [SingleMoleculeFootprinting Bioconductor package](https://bioconductor.org/packages/3.22/bioc/html/SingleMoleculeFootprinting.html): Versioned package source and documentation.
- [nf-smfONT](https://git.embl.org/grp-krebs/nf-smfont): BAM modification extraction, default LLR threshold, motif selection, and read-level filtering pipeline.
- [Modkit modification-call filtering](https://nanoporetech.github.io/modkit/filtering.html): Interpretation of pass/fail thresholds for modification probabilities.
- `/scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5`: ONT mESC accessibility H5 used for dimension, missingness, context, and compatibility diagnostics.
- `/scratch/users/diamant/smf_data/mm10.fa`: Reference used to validate genomic offsets and CpG/GpC contexts.
- `/home/groups/btrippe/diamant/containers/footprint-charter_bioc-3.22.sandbox`: Tested Sherlock Apptainer sandbox.
