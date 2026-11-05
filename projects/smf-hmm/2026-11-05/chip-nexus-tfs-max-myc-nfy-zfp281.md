---
title: More ChIP-nexus TFs (MAX::MYC, NF-Y, ZFP281) — none footprints like CTCF; NF-Y has a small narrow footprint
date: 2026-11-05
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - smf_hmm repo, branch nd/experiments, commit 5c657b9 (evaluation/tf_chip/configs/mesc_validation_footprint.yaml, references MA0059.2, MA0060.4, MA1630.3)
  - Caches /scratch/users/diamant/mesc_validation/benchmark_cache/mesc_tf_chip_validation_footprint_v1/{max_myc,nfy,zfp281}/{chr8,chr16}
  - Datasets and analyses /scratch/users/diamant/mesc_validation/covariation_windows/footprint_{chr8,chr16}-512-valid90-reads10 (analysis/, analysis/motif_call_diagnostics/)
  - Pipeline script /scratch/users/diamant/mesc_validation/footprint_tfs_2026-11-05.sh; log /scratch/users/diamant/mesc_validation/logs/footprint_tfs_2026-11-05.log
  - JASPAR API (https://jaspar.elixir.no/api/v1/matrix/<id>.meme), accessed 2026-11-05
  - User request, 2026-11-05
tags: [smf-hmm, tf-chip, benchmark, chip-nexus, footprint, max-myc, nfy, zfp281, prior-occupancy-kl]
---

# Summary

The user asked to evaluate other TFs that footprint and that have mESC ChIP-nexus: REST, RFX,
THAP11, ZFP524 and MAX::MYC.
- **Only MAX::MYC had data:** `mesc_cmyc_nexus`, 2 replicates, in
  `/scratch/groups/btrippe/ndiamant/julia_lab/bw`. No ChIP-nexus exists there, or elsewhere under
  the group's data directories, for REST, RFX, THAP11 or ZFP524.
- **Two substitutes were added**, other available ChIP-nexus TFs that may footprint:
  - NF-Y (NFYB ChIP-nexus, 5 replicates, with the NFYA CCAAT motif), a known strong footprinter;
  - ZFP281 (3 replicates).

The full benchmark ran on new chr8 and chr16 covariation datasets built for these TFs:
caches, window selection, 9 models' posteriors, controlled and protected-read metrics, and
meta-profiles.

**None of the three shows a CTCF-like binding footprint in the SMF reads:**
- **NF-Y:** a narrow footprint of about ±5 bp at the CCAAT box. Bound-site excess over context
  is +0.047 / +0.021 (chr8 / chr16), against about 0 at unbound sites; CTCF's is about +0.4.
  Bound NF-Y sites are also strongly accessible.
- **MAX::MYC:** no binding footprint. The CpG-containing E-box shows excess protection that is
  largest at unbound sites.
- **ZFP281:** no footprint.

On ChIP:
- **NF-Y** is the only TF where the HMMs clearly beat accessibility. The partial Spearman given
  accessibility is 0.37–0.48 on chr8, and also on protected reads. But the prior alone gives
  about the same, so this is the CNN prior's sequence model, not read evidence.
- **MAX::MYC:** about 0.
- **ZFP281:** like OCT4, its control-model signal disappears with MAP-rate KL targets.

CTCF remains the only TF here where footprint modelling clearly pays off on the ChIP benchmark.

# Key Points

- **Motif sites** on chr8 (genome / in the 5,000 selected windows):
  - MAX::MYC: 30,629 / 949;
  - NF-Y: 21,698 / 645;
  - ZFP281: 255,928 / 10,194. Its G-rich motif is very common at the cache's scan threshold
    (FPR 1e-4).
- chr16 has 2,488 windows, with 484 / 296 / 5,053 motif instances for MAX::MYC / NF-Y / ZFP281.
- **ChIP agreement** (Spearman with pooled ChIP-nexus, chr8 / chr16; controls as defined in
  the [controlled-metrics note](../2026-11-04/controlled-chip-metrics-and-context-offsets.md) and
  the [protected-read note](protected-read-chip-metric.md)):

| | MAX::MYC | NF-Y | ZFP281 |
|---|---|---|---|
| sites (all / ≥3 protected reads), chr8 | 949 / 842 | 645 / 207 | 10194 / 5030 |
| **raw** | | | |
| accessibility | 0.25 / 0.21 | 0.70 / 0.66 | 0.33 / 0.29 |
| motif score | 0.04 / 0.22 | 0.25 / 0.18 | 0.17 / 0.18 |
| phase_type | 0.06 / 0.11 | 0.64 / 0.54 | 0.28 / 0.29 |
| neural_global_rerun (original KL target) | 0.09 / 0.15 | 0.75 / 0.62 | 0.31 / 0.33 |
| neural_global_map_targets | 0.04 / 0.16 | 0.59 / 0.41 | 0.10 / 0.14 |
| single_tf | 0.04 / 0.14 | 0.64 / 0.46 | 0.24 / 0.27 |
| **partial given accessibility** | | | |
| motif score | 0.04 / 0.21 | 0.23 / 0.17 | 0.17 / 0.18 |
| neural_global_rerun | −0.02 / 0.06 | 0.48 / 0.28 | 0.16 / 0.23 |
| neural_global, prior only | −0.00 / 0.12 | 0.41 / 0.24 | 0.13 / 0.17 |
| neural_global_map_targets | −0.01 / 0.12 | 0.37 / 0.17 | 0.02 / 0.09 |
| **protected reads, partial** | | | |
| footprint fraction (reads only) | 0.00 / 0.07 | 0.08 / 0.10 | 0.03 / 0.03 |
| phase_type | −0.01 / 0.09 | 0.31 / 0.10 | 0.06 / 0.06 |
| neural_global_rerun | −0.02 / 0.10 | 0.47 / 0.16 | 0.04 / 0.08 |
| neural_global_map_targets | 0.00 / 0.10 | 0.45 / 0.16 | 0.03 / 0.08 |
| neural_global, prior only | −0.02 / 0.14 | 0.43 / 0.14 | 0.08 / 0.14 |

  Prior-only rows use the original `neural_global` checkpoint (2026-09-29), whose posteriors
  the rerun matches within about 0.03 on the main benchmark.
- **Footprint over context**: protection rise from 20–40 bp out to the motif center (±7 bp),
  minus the rise a ±3 bp context table predicts. All reads, by ChIP quartile Q1 → Q4.

| | chr8 Q1 / Q2 / Q3 / Q4 | chr16 Q1 / Q2 / Q3 / Q4 |
|---|---|---|
| MAX::MYC | +0.056 / +0.050 / +0.037 / +0.028 | +0.034 / +0.056 / +0.019 / −0.000 |
| NF-Y | −0.010 / −0.006 / +0.001 / **+0.047** | +0.014 / −0.011 / +0.023 / **+0.021** |
| ZFP281 | +0.033 / +0.021 / +0.011 / +0.014 | +0.026 / +0.005 / +0.004 / +0.030 |

  For comparison, CTCF is +0.42 / +0.46 at Q4 and OCT4 +0.02 / +0.01.

# Details

**Interpretation (inferred):**
- **NF-Y:** bound sites are in very accessible DNA (P(protected) about 0.25, against 0.65 at
  unbound sites in the meta-profile), with a sharp narrow peak at the motif. The HMMs' advantage
  over accessibility comes from the CNN prior: prior-only equals posterior on protected reads
  (0.43 vs 0.45–0.47 on chr8). It recognises CCAAT boxes in open DNA. The footprint itself is too
  small and narrow, with few calls under it, to drive per-read calls.
- **MAX::MYC:** the E-box CACGTG contains a CpG. Excess protection there falls with ChIP, so it
  most likely reflects CpG-specific calling or endogenous CpG methylation rather than binding.
  Unverified: whether CpG calls in this assay measure accessibility only.
- **ZFP281:** behaves like ESRRB and OCT4. The control's partial correlation (0.16–0.23)
  vanishes with MAP-rate KL targets (0.02–0.09), so it came from the weak-TF floor of the
  [occupancy-target note](../2026-11-04/controlled-chip-metrics-and-context-offsets.md).
- **Overall:** in this mESC SMF data, ChIP-nexus intensity for these TFs mostly tracks
  accessibility and motif strength rather than per-molecule protection. CTCF is the exception.

**Choices and caveats:**
- REST, RFX, THAP11 and ZFP524 could not be evaluated: there is no ChIP-nexus for them in
  `julia_lab/bw` or elsewhere under `/scratch/groups/btrippe/ndiamant`. Public ChIP-seq (e.g.
  ENCODE REST in mESC) was not searched or downloaded, pending the user's decision.
- NF-Y pairs NFYB ChIP with the NFYA motif `MA0060.4`, since NF-Y binds CCAAT as a heterotrimer.
- `MA0059.2` (MAX::MYC) is built from only 21 sites.
- The new datasets use the main datasets' settings: 512 bp windows, valid fraction ≥ 0.9, ≥ 10
  reads, ≤ 5,000 windows, seed 123, no motif-flank constraint. Windows were selected for these
  TFs, so they differ from the main datasets.
- NF-Y's protected-read samples are small (207 / 106 sites), so its chr16 numbers are noisy.

# Commands

Everything ran in Slurm job 51132065 (sh04-06n09, 1 L40S, 4 CPUs, 24 GB): about 75 minutes
end to end, from 09:54 to 11:09.

```bash
cd /home/users/diamant/repos/smf_hmm
bash /scratch/users/diamant/mesc_validation/footprint_tfs_2026-11-05.sh   # caches, windows, 9 posteriors, analysis, diagnostics
```

The script runs, per TF and chromosome:
- `evaluation.tf_chip.benchmark_cache --config evaluation/tf_chip/configs/mesc_validation_footprint.yaml`;
- `evaluation.tf_chip.select_covariation_windows` with the main settings.

Then, per chromosome:
- `smf_hmm.cli.infer_posterior` for:
  - the geometric and phase_type baselines;
  - `neural_global` and its `--prior-only` reference;
  - `neural_global_rerun`;
  - `neural_global_map_targets` and its `--prior-only` reference;
  - `single_tf`;
  - `single_tf_147_merged_cgc`, with `--assay-type cg_gc_accessibility_merged_cgc`;
- `evaluation.tf_chip.analyze_covariation`;
- `evaluation.tf_chip.motif_call_diagnostics` with `--model neural_global_rerun --model neural_global_map_targets`.

The tables above came from short read-outs of `individual_correlations.tsv` and the diagnostics'
`results.json`.

# Related Notes

- [Protected-read ChIP metric](protected-read-chip-metric.md): The protected-read variant used here, and the CTCF and OCT4 reference numbers.
- [Controlled ChIP metrics and occupancy targets](../2026-11-04/controlled-chip-metrics-and-context-offsets.md): The weak-TF floor and the MAP-target model that explain the ZFP281 pattern.
- [ChIP evaluation: single-TF vs multi-TF](../2026-11-03/chip-eval-single-tf-vs-multi-tf.md): The footprint-over-context measure and the CTCF/OCT4 aggregate footprints compared here.

# Open Questions

- Should public ChIP-seq be used for REST and the other requested footprinting TFs, e.g. ENCODE
  mESC REST? That needs the user's go-ahead, and a decision on mixing ChIP-seq with ChIP-nexus.
- Is the NF-Y footprint learnable with more receptive field, or is it too narrow for the calls
  available under it?
- Do CpG calls in this assay carry endogenous methylation? That would explain the MAX::MYC
  E-box pattern.

# Sources

- Config and motifs: smf_hmm `5c657b9` (`evaluation/tf_chip/configs/mesc_validation_footprint.yaml`, `evaluation/tf_chip/references/{MA0059.2,MA0060.4,MA1630.3}.meme`, `references/README.md`)
- ChIP-nexus bigWigs: `/scratch/groups/btrippe/ndiamant/julia_lab/bw/{mesc_cmyc_nexus,mesc_nfyb_nexus,mesc_zfp281_nexus}_*.bw`
- Outputs: `/scratch/users/diamant/mesc_validation/covariation_windows/footprint_{chr8,chr16}-512-valid90-reads10/`
- JASPAR: https://jaspar.elixir.no/api/v1/matrix/MA0059.2.meme, https://jaspar.elixir.no/api/v1/matrix/MA0060.4.meme, https://jaspar.elixir.no/api/v1/matrix/MA1630.3.meme
