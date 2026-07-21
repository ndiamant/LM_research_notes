---
title: CTCF motif experiment summary
date: 2026-07-21
project: hiddenfoot
agent: Codex
status: draft
sources:
  - /home/users/diamant/repos/HiddenFoot/source/hiddenfoot.cpp
  - /home/users/diamant/repos/HiddenFoot/data/pwms_clusters_jaspar.txt
  - /home/users/diamant/repos/HiddenFoot/data/pwms_ctcf_jaspar.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_DNA.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_SMF.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr1:6405368-6407482_DNA.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr1:6405368-6407482_SMF.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr16:29559618-29561732_DNA.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr16:29559618-29561732_SMF.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr16_29559618_29561732_default_plus_ctcf.states.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr8:71708782-71710896_DNA.txt
  - /home/users/diamant/repos/HiddenFoot/data/chr8:71708782-71710896_SMF.txt
  - /home/users/diamant/repos/HiddenFoot/results_chr8_consensus_all_tfs/chr8_71708782_71710896_all_tfs.map_state_labels.tsv
  - /home/users/diamant/repos/HiddenFoot/results_chr8_consensus_all_tfs/chr8_71708782_71710896_all_tfs.map_states.tsv
  - /home/users/diamant/repos/HiddenFoot/results_chr8_consensus_all_tfs/chr8_71708782_71710896_all_tfs.map_starts.tsv
tags:
  - hiddenfoot
  - ctcf
  - mesc
  - map-decoder
  - motif-scanning
  - experiments
---

# Summary

We ran several HiddenFoot fits on mESC single-molecule footprinting inputs to test whether CTCF-family motifs would be selected and called as binding states. Default all-TF runs on the first chr16 locus and the matched chr1 locus did not select CTCF-family states. A CTCF-only diagnostic run could force the state space to `free`, `nucleosome`, and `CTCF`, but that result is not directly comparable to all-TF fits because it used only one motif class and an extremely permissive motif cutoff.

For the chr16 region containing the sequence `CGCCCCCTAGTGG`, the default motif scanner still did not include CTCF-family states because the matching PWM scores were below HiddenFoot's default `-cutoffwms 0.001`. A forced all-TF run that added two CTCF-family candidate states produced sparse hard MAP calls: `INSM1_CTCFL` had 5 MAP starts and `CTCF_CTCF_CTCF` had 0.

For the chr8 region `chr8:71708782-71710896`, which contains the exact shorter `INSM1_CTCFL` consensus `TAGCAGGGGGCGC`, HiddenFoot selected `INSM1_CTCFL` at the default cutoff. Even in this consensus motif region, CTCF-family hard MAP calls remained sparse: 26 start calls across 232 molecules, spanning 338 bases. Motif presence was sufficient to include a candidate state, but posterior/MAP binding calls remained conservative.

# Key Points

- The MAP decoder addition writes integer state calls per molecule and position in `*.map_states.tsv`, start-only calls in `*.map_starts.tsv`, and integer-to-state-name metadata in `*.map_state_labels.tsv`.
- The default all-TF chr16 and chr1 runs did not include CTCF-family states because no CTCF-family motif passed the scanner threshold in those regions.
- The chr16 sequence `CGCCCCCTAGTGG` looked like a CTCF consensus-like site by substring, but the relevant PWM scores were below the default cutoff: `CTCF_CTCF_CTCF` scored `6.64156e-08` and `INSM1_CTCFL` scored `0.000416358`.
- Lowering the all-TF cutoff to `0.0001` on the chr16 CTCF-like region selected more motifs but produced `NaN` likelihoods, so that fit was not usable.
- Forcing CTCF-family candidates into the chr16 CTCF-like region produced marginal posterior signal for some molecules, but very few hard MAP state calls.
- The chr8 exact `INSM1_CTCFL` consensus was selected by the default scanner with score `0.5`, but the fitted model still produced only 26 `INSM1_CTCFL` MAP starts across 232 molecules.

# Details

| Run | Output prefix | State selection | CTCF-family result |
| --- | --- | --- | --- |
| chr16:50073256-50075370 all TFs | `results/chr16_50073256_50075370` | 57 selected states, 5001 elements | No CTCF-family state selected |
| chr16:50073256-50075370 CTCF-only | `results_ctcf/chr16_50073256_50075370_ctcf_only` | `free`, `nucleosome`, `CTCF`; 7062 elements | Diagnostic only; 819 CTCF MAP starts under CTCF-only state space |
| chr1:6405368-6407482 matched all TFs | `results_chr1_matched_all_tfs/chr1_6405368_6407482_all_tfs` | 67 selected states, 5037 elements | No CTCF-family state selected |
| chr16:29559618-29561732 default all TFs | `results_chr16_ctcf_locus_all_tfs/chr16_29559618_29561732_all_tfs` | 61 selected states, 5008 elements | No CTCF-family state selected |
| chr16:29559618-29561732 forced CTCF | `results_chr16_ctcf_locus_all_tfs_forced_ctcf/chr16_29559618_29561732_all_tfs_forced_ctcf` | 63 selected states, 5010 elements | `INSM1_CTCFL`: 5 MAP starts; `CTCF_CTCF_CTCF`: 0 MAP starts |
| chr8:71708782-71710896 consensus all TFs | `results_chr8_consensus_all_tfs/chr8_71708782_71710896_all_tfs` | 62 selected states, 5060 elements | `INSM1_CTCFL` selected; 26 MAP starts and 338 MAP bases |

There is also an accidental mixed-input run using chr1 DNA with chr16 SMF reads at `results_chr1_all_tfs/chr1_6405368_6407482_chr16_smf_all_tfs`. That run should be treated as an input mismatch and ignored for biological interpretation.

The main all-TF command pattern used for the matched runs was:

```bash
./binary/hiddenfoot_omp run \
  -s 'data/<region>_DNA.txt' \
  -m 'data/<region>_SMF.txt' \
  -w data/pwms_clusters_jaspar.txt \
  -o <output-prefix> \
  -numthreads 4
```

The CTCF-only diagnostic command used:

```bash
./binary/hiddenfoot_omp run \
  -s 'data/chr16:50073256-50075370_DNA.txt' \
  -m 'data/chr16:50073256-50075370_SMF.txt' \
  -w data/pwms_ctcf_jaspar.txt \
  -o results_ctcf/chr16_50073256_50075370_ctcf_only \
  -numthreads 4 \
  -cutoffwms 1e-300
```

For the chr16 CTCF-like region, the exact substring `CGCCCCCTAGTGG` was not enough to trigger default CTCF-family inclusion. Under the PWMs used by `data/pwms_clusters_jaspar.txt`, the best CTCF-family candidates were below the default threshold:

- `CTCF_CTCF_CTCF`, offset 1422, length 35, score `6.64156e-08`
- `INSM1_CTCFL`, offset 1424, length 13, score `0.000416358`

The shorter `INSM1_CTCFL` consensus substring to search for is `TAGCAGGGGGCGC`, or its reverse complement `GCGCCCCCTGCTA`. The chr8 region contains the forward consensus at offset 1470 and was selected by the default scanner:

```text
1470 248 13 0.5 INSM1_CTCFL
```

However, selected candidate states do not guarantee hard MAP binding calls. In the chr8 consensus run, `INSM1_CTCFL` was state ID 24 in `*.map_state_labels.tsv`, and `*.map_starts.tsv` contained 26 start calls for that state across 232 molecules. The strongest marginal `INSM1_CTCFL` posterior was `0.919607` at molecule 156, offset 790, not at the exact consensus offset 1470. This suggests that the emission evidence from the SMF reads and competition with other states can dominate the sequence prior when producing hard MAP paths.

Some all-TF runs failed to write a few per-motif `hfbinding` or `hfprofile` files because raw JASPAR cluster names exceeded filesystem filename length limits. The compact MAP outputs and core state/parameter/probability files were still written.

# Related Notes

- [Compile and binary verification](../2026-07-20/compile-and-binary-verification.md): Documents compiler overrides used for HiddenFoot builds.
- [mESC data preparation and MAP decoder addition](../2026-07-20/mesc-data-preparation-and-map-decoder.md): Records data export provenance and the MAP decoder output contract.

# Open Questions

- Should future CTCF analyses use a curated CTCF PWM and an empirically chosen threshold instead of the broad clustered JASPAR motif set?
- Why are hard MAP calls sparse at the exact chr8 `INSM1_CTCFL` consensus site despite motif selection?
- Should visualizations compare marginal posterior binding probabilities against hard MAP paths, rather than using MAP calls alone?
- Are local SMF methylation patterns in the consensus-containing molecules consistent with CTCF occupancy or better explained by nucleosome/free states?

# Sources

- `/home/users/diamant/repos/HiddenFoot/source/hiddenfoot.cpp`
- `/home/users/diamant/repos/HiddenFoot/data/pwms_clusters_jaspar.txt`
- `/home/users/diamant/repos/HiddenFoot/data/pwms_ctcf_jaspar.txt`
- `/home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_DNA.txt`
- `/home/users/diamant/repos/HiddenFoot/data/chr16:50073256-50075370_SMF.txt`
- `/home/users/diamant/repos/HiddenFoot/data/chr1:6405368-6407482_DNA.txt`
- `/home/users/diamant/repos/HiddenFoot/data/chr1:6405368-6407482_SMF.txt`
- `/home/users/diamant/repos/HiddenFoot/data/chr16:29559618-29561732_DNA.txt`
- `/home/users/diamant/repos/HiddenFoot/data/chr16:29559618-29561732_SMF.txt`
- `/home/users/diamant/repos/HiddenFoot/data/chr16_29559618_29561732_default_plus_ctcf.states.txt`
- `/home/users/diamant/repos/HiddenFoot/data/chr8:71708782-71710896_DNA.txt`
- `/home/users/diamant/repos/HiddenFoot/data/chr8:71708782-71710896_SMF.txt`
- `/home/users/diamant/repos/HiddenFoot/results_chr8_consensus_all_tfs/chr8_71708782_71710896_all_tfs.map_state_labels.tsv`
- `/home/users/diamant/repos/HiddenFoot/results_chr8_consensus_all_tfs/chr8_71708782_71710896_all_tfs.map_states.tsv`
- `/home/users/diamant/repos/HiddenFoot/results_chr8_consensus_all_tfs/chr8_71708782_71710896_all_tfs.map_starts.tsv`
