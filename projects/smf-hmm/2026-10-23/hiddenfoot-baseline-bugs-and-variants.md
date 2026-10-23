---
title: HiddenFoot baseline: bugs, fixed reruns, and favorable-setting variants
date: 2026-10-23
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - smf_hmm/third_party/hiddenfoot/VENDORED.md
  - smf_hmm/third_party/hiddenfoot/source/hiddenfoot.cpp
  - smf_hmm/evaluation/tf_chip/workflow/HIDDENFOOT.md
  - smf_hmm/evaluation/tf_chip/workflow/run_hiddenfoot.py
  - smf_hmm/evaluation/tf_chip/workflow/hiddenfoot_variants.py
  - smf_hmm/evaluation/tf_chip/analyze_covariation.py
  - $SCRATCH/mesc_validation/covariation_windows/{chr8,chr16}-512-valid90-reads10/analysis_hiddenfoot_fixed
  - $SCRATCH/mesc_validation/hiddenfoot_variants/chr16_subset500/analysis
  - https://github.com/MolinaLab-IE/HiddenFoot (commit d1499e8)
tags: [hiddenfoot, baseline, tf-chip, benchmark, smf, reproducibility]
---

# Summary

HiddenFoot (official repo, MIT) was added as an external baseline to the
posterior-based TF/ChIP benchmark on mESC validation chromosomes chr8 and chr16.
It is far below the HMMs, and is negatively correlated with ChIP for KLF4 and
OCT4.

Investigating the gap uncovered three upstream bugs:

- RNG seeding from `clock()`.
- An integer-division SGD batch count. It skipped SGD entirely on windows with
  fewer than 64 reads, which is 55% of chr16 windows.
- An uninitialized gradient row.

All three are fixed in the vendored copy. Rerunning with the fixed build
improved CTCF only modestly (chr8 0.24 → 0.30, chr16 0.19 → 0.25) and left
KLF4/OCT4 negative.

In the most favorable setting tried (fit and evaluate on the same 500 chr16
windows), three changes did not close the gap:

- MCMC-averaged posteriors.
- Consensus parameters shared across windows.
- A four-motif "oracle" PWM set.

The best HiddenFoot variant per TF reached CTCF 0.44, ESRRB 0.22, KLF4 −0.05
and OCT4 0.16. The HMM neural_global model reached 0.71, 0.31, 0.15 and 0.57.

Diagnosis, inferred from state occupancies at motifs: with the full 258-motif
set, HiddenFoot's TF states absorb nucleosome-scale protection. At motifs its
mean P(TF) is 0.30 and P(nucleosome) is 0.24, against 0.04 and 0.53 for the
phase-type HMM. Its "TF" signal therefore tracks closed chromatin and
anticorrelates with ChIP. Restricting to four motifs restores the
nucleosome/accessibility pattern, but then almost no TF is called (0.03), and
for CTCF that little TF anticorrelates with ChIP.

Conclusion: report the fixed default HiddenFoot as the baseline, and treat the
variants as a supplementary robustness check.

# Key Points

- **Benchmark population**: covariation windows (512 bp, ≥90% valid, ≥10
  reads). Spearman ρ between pooled ChIP and read-averaged motif occupancy,
  over unique motif sites. P(TF) sums all TF groups. See
  [the chr16 covariation note](../2026-09-28/chr16-covariation-chip-width-and-window-selection.md).
- **Upstream bugs** (line references are to vendored `source/hiddenfoot.cpp`):
  1. **Seeding.** `default_random_engine generator((unsigned) clock())` at five
     sites. `clock()` has 10 ms resolution and sums CPU time over threads.
     Runs were therefore irreproducible even single-threaded (up to 0.15 in
     P(TF)), and with 4 threads differed by up to 0.42. Fixed with one shared
     generator seeded with 1.
  2. **`numbatches = ceil(nummol/batchsize)`** divides integers. With fewer than
     64 reads there are 0 SGD batches, so MCMC starts at the initial
     parameters. With 64–127 reads one random batch is used per epoch. Fixed to
     integer ceiling. This changes P(TF) by up to 0.9 on affected windows. It is
     the most consequential bug.
  3. **`forward_data_grad`/`forward_conf_grad`** never zero row 0 of the
     gradient arrays (valgrind). The effect is small (≤0.02).

  Two bounds fixes were already applied before these: the DNA-string loop and
  the msites search.
- **Even after the fixes, results depend on the seed.** Seeds 1 and 2 differ by
  up to 0.34 in P(TF). The default posterior uses only the final MCMC sample of
  an unconverged chain.
- **Pooled joint fit is infeasible as-is.** HiddenFoot's forward recursion
  multiplies unscaled probabilities. Concatenating 5 windows (~4 kb) gives NaN
  after the first SGD update; 2–3 windows (≤2.1 kb) are fine. A true pooled fit
  would need a rescaled or log-space forward-backward. "Pooled" below is
  therefore a consensus-parameter heuristic, not a joint fit.
- **The input encoding is correct.** At initialization, HiddenFoot's `qU`
  (P(emit 0 | unbound)) is 0.05, so 0 means protected. The exporter maps
  dataset 1 → 0 and 0 → 1. Group indices are correct too: state 0 is free,
  state 1 is nucleosome, states ≥2 are PWMs.

# Details

## Rerun with fixed build (full chromosomes)

Summed-TF Spearman ρ, HiddenFoot before → after fixes. Other models shown for
context (chr8 / chr16):

| TF | HF chr8 | HF chr16 | neural_global | phase_type | geometric |
|---|---|---|---|---|---|
| CTCF | 0.24 → 0.30 | 0.19 → 0.25 | 0.67 / 0.69 | 0.64 / 0.64 | 0.60 / 0.63 |
| ESRRB | 0.08 → 0.08 | 0.08 → 0.05 | 0.25 / 0.30 | 0.25 / 0.30 | 0.25 / 0.30 |
| KLF4 | −0.23 → −0.18 | −0.20 → −0.17 | 0.24 / 0.21 | 0.24 / 0.20 | 0.11 / 0.08 |
| OCT4 | −0.06 → −0.05 | −0.05 → −0.06 | 0.50 / 0.57 | 0.32 / 0.38 | 0.26 / 0.30 |

Site counts: chr8 has 3,780 / 1,674 / 6,655 / 810 for CTCF / ESRRB / KLF4 /
OCT4; chr16 has 1,754 / 859 / 3,329 / 520.

## Favorable-setting variants (chr16 subset, 500 random windows, seed 0)

The subset has 30,708 reads and 333 / 185 / 669 / 112 sites for CTCF / ESRRB /
KLF4 / OCT4. Uncertainty is about ±0.1 for ESRRB and OCT4.

The variants:

- **HF default**: per-window fit, last MCMC sample.
- **MCMC mean**: average of the marginals over 20 samples from the second half
  of the chain.
- **Consensus**: per motif, and for β, qU and qB, the median over the subset's
  windows of each window's mean log-parameter. That single vector is used for
  every window.
- **Oracle**: per-window fit with only the CTCF, POU (incl. Pou5f1::Sox2),
  KLF/SP/EGR and nuclear-receptor (incl. ESRRB) JASPAR clusters, then MCMC
  mean.

| TF | phase_type | neural_global | HF default | HF MCMC mean | HF consensus | HF oracle |
|---|---|---|---|---|---|---|
| CTCF | 0.68 | **0.71** | 0.34 | 0.31 | 0.44 | −0.19 |
| ESRRB | **0.36** | 0.31 | −0.07 | −0.08 | −0.04 | 0.22 |
| KLF4 | 0.13 | **0.15** | −0.11 | −0.12 | −0.12 | −0.05 |
| OCT4 | 0.39 | **0.57** | −0.20 | −0.17 | −0.11 | 0.16 |

Mean state probability at motif bases (all TFs pooled):

| | P(open) | P(nucleosome) | P(TF) |
|---|---|---|---|
| phase_type | 0.43 | 0.53 | 0.04 |
| HF default | 0.46 | 0.24 | 0.30 |
| HF consensus | 0.46 | 0.30 | 0.25 |
| HF oracle | 0.40 | 0.57 | 0.03 |

Spearman ρ of P(open) and P(nucleosome) at motifs with ChIP:

- **OCT4**: phase_type gives +0.46 and −0.46. HF default gives +0.22 and +0.22.
  HF oracle gives +0.45 and −0.46, matching the HMM.
- **CTCF**: phase_type gives +0.25 and −0.57. HF default gives −0.11 and −0.17.
  HF oracle gives +0.41 and −0.43.

This supports the absorption diagnosis. The full motif set lets overlapping PWM
hits explain protected stretches. That displaces the nucleosome state, and TF
becomes a closed-chromatin signal. This is an interpretation of occupancy
patterns, not a direct proof.

**Inferred, not tested:** the oracle's poor CTCF performance may relate to
HiddenFoot's motif-selection cutoff. The earlier
[CTCF motif experiment](../../hiddenfoot/2026-07-21/ctcf-experiment-summary.md)
found CTCF-family states often not selected at the default `-cutoffwms 0.001`,
and found sparse CTCF calls even where selected.

**Data scale** likely matters too. HiddenFoot's example amplicon has 3,715
molecules over 276 bp, while benchmark windows have about 60–90 reads over
512 bp. One window (`w0001179`) selected 34 motifs, i.e. 1,569 candidate
elements and about 37 parameters, from 69 reads.

## Provenance and commands

smf_hmm commits (local on `main` at the time of writing; not yet pushed):

| Commit | Change |
|---|---|
| `62bb436` | Vendor HiddenFoot d1499e8 |
| `8a8552d` | Bounds fixes and `-grouped` output |
| `4d121c7` | `run_hiddenfoot` workflow |
| `9717a8e` | Seed and gradient-row fixes |
| `d9b4e18` | `numbatches` fix |
| `b46ec1f` | `HIDDENFOOT.md` output provenance |
| `2553a70` | `hiddenfoot_variants.py` |

Earlier output schema mismatch: pre-fix `hiddenfoot.h5` files lacked
`schema_version` and `prepared` provenance. They were patched in place after
verifying the dataset SHA-256 and read IDs. `run_hiddenfoot.py` now writes both.

Fixed-build reruns (Slurm, partitions `hns,stat,akundaje`, 4 CPUs, 8 GB). The
binary was `$SCRATCH/mesc_validation/hiddenfoot_build/hiddenfoot_d9b4e18`, built
with `g++ -O3 -DUSE_OMP -fopenmp`, GCC 14.2.0. Scripts are in
`$SCRATCH/mesc_validation/hiddenfoot_build/` (`run_hiddenfoot.sbatch`,
`analyze.sbatch`).

| Job | Step | Runtime |
|---|---|---|
| 49215217 | chr16 HiddenFoot | 34 min |
| 49215220 | chr8 HiddenFoot | 1 h 38 min |
| 49216545 | chr16 analysis | `afterok` dependency |
| 49216556 | chr8 analysis | `afterok` dependency |

```bash
B=$SCRATCH/mesc_validation/hiddenfoot_build; D=$SCRATCH/mesc_validation/covariation_windows
sbatch --time=03:00:00 --job-name=hiddenfoot_chr16 --output=$B/logs/%x-%j.out $B/run_hiddenfoot.sbatch $D/chr16-512-valid90-reads10
sbatch --time=08:00:00 --job-name=hiddenfoot_chr8 --output=$B/logs/%x-%j.out $B/run_hiddenfoot.sbatch $D/chr8-512-valid90-reads10
sbatch --time=01:00:00 --dependency=afterok:<hiddenfoot job> --output=$B/logs/%x-%j.out $B/analyze.sbatch <root>
```

Outputs:

- `posteriors/hiddenfoot.h5` and `analysis_hiddenfoot_fixed/` under each
  chromosome root.
- The pre-fix chr8 file is kept as `posteriors/hiddenfoot_prefix_build.h5`.
- In the shared advisor folder `$GROUP_SCRATCH/smf_hmm_tf_chip_chr16_2026-10-01`,
  the fixed chr16 files replaced the old ones, and the old ones went to
  `superseded/`.

The variant experiment ran interactively in allocation 49031037 (about 30 min).
Its outputs are in `$SCRATCH/mesc_validation/hiddenfoot_variants/chr16_subset500/`
(`posteriors/`, `analysis/`). Commands are in `HIDDENFOOT.md`, section
"Diagnostic variants", using `python -m evaluation.tf_chip.workflow.hiddenfoot_variants`
with subcommands subset, select-pwms, shard, consensus, pooled and merge. The
original exploratory script and log are in
`$SCRATCH/mesc_validation/hiddenfoot_variants/` (`hf_variants.py`,
`logs/pipeline.log`). The committed version was checked against them: identical
oracle PWM file and consensus, and shard posteriors within about 1e-7.

# Related Notes

- [Chr16 covariation, ChIP width, and window selection](../2026-09-28/chr16-covariation-chip-width-and-window-selection.md):
  defines the covariation-window population and analysis used here.
- [Posterior-mean TF/ChIP evaluation](../2026-09-23/posterior-mean-tf-chip-evaluation.md):
  the posterior-file evaluation conventions HiddenFoot was adapted to.
- [Motif-window posterior inference pipeline](../2026-09-23/motif-window-posterior-inference.md):
  the shared posterior schema HiddenFoot outputs must match.
- [HiddenFoot CTCF motif experiment](../../hiddenfoot/2026-07-21/ctcf-experiment-summary.md):
  earlier evidence that CTCF-family states are rarely selected or called. This
  may explain the oracle's CTCF result.
- [HiddenFoot compile and binary verification](../../hiddenfoot/2026-07-20/compile-and-binary-verification.md):
  the macOS makefile issue that led to the Linux build command used here.
- [HiddenFoot mESC data prep and MAP decoder](../../hiddenfoot/2026-07-20/mesc-data-preparation-and-map-decoder.md):
  origin of the local MAP-decoder edits in the user's HiddenFoot checkout. They
  are not part of the vendored copy.

# Open Questions

- Would a true pooled fit, after adding rescaling to HiddenFoot's
  forward-backward, behave differently from the consensus heuristic?
- How much of the oracle's CTCF failure is due to motif-selection cutoffs
  (`-cutoffwms`) versus footprint emission modeling?
- Would fitting on training chromosomes (the HMMs' data source) change
  anything? The favorable chr16 → chr16 setting already failed, so this seems
  unlikely to rescue HiddenFoot.
- Should HiddenFoot's multi-sample MCMC averaging be the default reported
  variant? It made no measurable difference here.

# Sources

- HiddenFoot upstream: https://github.com/MolinaLab-IE/HiddenFoot (commit
  `d1499e8`; MIT, Nacho Molina). Vendored at `smf_hmm/third_party/hiddenfoot/`.
- `smf_hmm/third_party/hiddenfoot/VENDORED.md`: list of local fixes.
- `smf_hmm/evaluation/tf_chip/workflow/HIDDENFOOT.md`,
  `run_hiddenfoot.py`, `hiddenfoot_variants.py`: workflow and variants.
- Analysis tables: `individual_correlations.tsv` and
  `unique_motif_summaries.tsv.gz` under the output roots listed above.
- Sub-agent bug investigation: valgrind 3.14 (`module load devel
  valgrind/3.14.0`), ASan/UBSan builds, seed-printing builds. The scratch
  artifacts in `/lscratch/diamant/49031037/hf_debug/` are node-local and
  temporary.
