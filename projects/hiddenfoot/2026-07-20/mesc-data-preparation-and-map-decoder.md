---
title: mESC data preparation and MAP decoder addition
date: 2026-07-20
project: hiddenfoot
agent: Codex
status: draft
sources:
  - /home/users/diamant/repos/smf_diffusion/exploration/hiddenfoot_data.ipynb
  - /home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_DNA.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_SMF.txt
  - /home/users/diamant/repos/HiddenFoot/source/hiddenfoot.cpp
  - /home/users/diamant/repos/HiddenFoot/results/chr16_50073256_50075370.map_state_labels.tsv
  - /home/users/diamant/repos/HiddenFoot/results/chr16_50073256_50075370.map_states.tsv
  - /home/users/diamant/repos/HiddenFoot/results/chr16_50073256_50075370.map_starts.tsv
tags:
  - hiddenfoot
  - mesc
  - smf
  - viterbi
  - map-decoder
  - visualization
---

# Summary

The mESC HiddenFoot test data for `chr16:50073256-50075370` was prepared in `/home/users/diamant/repos/smf_diffusion/exploration/hiddenfoot_data.ipynb`. HiddenFoot was then extended with a Viterbi-style MAP decoder that writes compact integer state-call matrices plus a state-label lookup table for visualization.

# Key Points

- The source notebook exports `/home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_DNA.txt` and `/home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_SMF.txt`.
- The notebook uses `SparseSMFDataModule` with mm10 FASTA and ONT mESC accessibility-region HDF5 inputs, selects validation region `chr16:50073256-50075370`, densifies the sparse SMF matrix, converts missing values from `-1` to `2`, and writes HiddenFoot-compatible whitespace-delimited files.
- HiddenFoot was extended in `source/hiddenfoot.cpp` with `get_map_state_calls`, called after `get_binding` whenever `run` or `bnd` requests binding outputs.
- The MAP decoder writes `*.map_states.tsv`, `*.map_starts.tsv`, and `*.map_state_labels.tsv`, avoiding the raw motif-name filename issue in per-motif posterior files.
- The first chr16 MAP output has 128 molecule rows and 2114 base-position columns, matching the prepared region length.

# Details

The data-preparation notebook configures the sparse mESC data module with:

```text
fasta_path="/scratch/users/diamant/smf_data/mm10.fa"
h5_path="/scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5"
reads_per_region=512
h5_region_length=2114
output_len=2048
region = "chr16:50073256-50075370"
```

The notebook then selects the region from `data.val_regions`, obtains `data.val_dna`, `data.val_mats`, and `data.val_pos`, converts DNA to a string, and writes the HiddenFoot files. Missing SMF entries are encoded as `2` before export, which is consistent with the current HiddenFoot parser because only `0` and `1` contribute to methylation counts.

The MAP decoder is a Viterbi-style dynamic program over the candidate state elements already selected by HiddenFoot. For each molecule, it maximizes the globally consistent non-overlapping path score:

```text
previous_score + log(probM(state, molecule)) + log(probB(state))
```

It uses the final transformed parameter vector available to the binding stage. In the chr16 test run, the decoder was run using:

```bash
./binary/hiddenfoot_omp bnd \
  -s 'data/chr16:50073256-50075370_DNA.txt' \
  -m 'data/chr16:50073256-50075370_SMF.txt' \
  -w data/pwms_clusters_jaspar.txt \
  -p results/chr16_50073256_50075370.params_MCMC.txt \
  -o results/chr16_50073256_50075370 \
  -numthreads 4
```

The new output contract is:

- `*.map_state_labels.tsv`: `state_id`, original `pwm_id`, footprint `length`, `state_type`, and motif name.
- `*.map_states.tsv`: one row per molecule, one integer state ID per base position. Occupied intervals are filled across every covered base.
- `*.map_starts.tsv`: one row per molecule, one integer state ID only at decoded event starts; `0` indicates no non-free start at that base.

One implementation detail to remember: `0` is also the `free` state ID in `map_state_labels.tsv`. In `map_starts.tsv`, zeros are used as the no-start background, so a decoded free-state start is not represented as a separate event.

# Related Notes

- [Compile and binary verification](compile-and-binary-verification.md): Establishes the GCC/OpenMP compile commands used before adding and running the MAP decoder.

# Open Questions

- Should MAP decoding be controlled by a command-line flag, or should it remain a default side effect of `run` and `bnd`?
- Should `map_starts.tsv` use `-1` for no-start so that a start of state `0` can be represented explicitly?
- Should future visualization code use region-relative offsets or emit genomic coordinates directly in a long-form BED-like MAP file?

# Sources

- `/home/users/diamant/repos/smf_diffusion/exploration/hiddenfoot_data.ipynb`: Notebook that prepared the mESC chr16 DNA and SMF matrix files for HiddenFoot.
- `/home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_DNA.txt`: Exported sequence input for the chr16 region.
- `/home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_SMF.txt`: Exported HiddenFoot-format SMF matrix for the chr16 region.
- `/home/users/diamant/repos/HiddenFoot/source/hiddenfoot.cpp`: Adds `get_map_state_calls` and invokes it after binding-profile generation.
- `/home/users/diamant/repos/HiddenFoot/results/chr16_50073256_50075370.map_state_labels.tsv`: Integer state-label lookup output from the chr16 run.
- `/home/users/diamant/repos/HiddenFoot/results/chr16_50073256_50075370.map_states.tsv`: Per-molecule base-level MAP state ID matrix from the chr16 run.
- `/home/users/diamant/repos/HiddenFoot/results/chr16_50073256_50075370.map_starts.tsv`: Per-molecule MAP start-event ID matrix from the chr16 run.
