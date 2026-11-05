---
title: Accessibility-controlled ChIP metrics; the OCT4 weak-TF floor comes from the occupancy KL target, not context bias
date: 2026-11-04
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /scratch/users/diamant/mesc_validation/covariation_windows/{chr8,chr16}-512-valid90-reads10/analysis_2026-11-04, analysis_2026-11-04_context3 and analysis_2026-11-04_map_targets (individual_correlations.tsv, unique_motif_summaries.tsv.gz, motif_call_diagnostics/)
  - Runs /scratch/users/diamant/2026-11-03_fixed_blocks_global_kl3_{context3,rerun} (wandb 3727hyo4, 19wq03el); /scratch/users/diamant/2026-11-04_fixed_blocks_global_kl3_map_targets (wandb 8cz8h7q6)
  - smf_hmm repo, branch nd/experiments, commits d160bfc, 8626d7f, 4436579, 53b9552, 61f2e83
  - User discussion, 2026-11-03 and 2026-11-04
tags: [smf-hmm, tf-chip, benchmark, metrics, accessibility, prior-only, context-bias, emissions, oct4, ctcf, klf4, prior-occupancy-kl, weak-tf]
---

# Summary

This note follows the [ChIP evaluation note](../2026-11-03/chip-eval-single-tf-vs-multi-tf.md).
There, the multi-TF model's OCT4 advantage was traced to a weak-TF "floor": about 0.4
P(TF) placed on fully accessible reads at OCT4 motifs. The hypothesis was that a
sequence-context bias at those motifs causes it.

**New metrics.** We added accessibility-controlled ChIP metrics and prior-only references to
the benchmark. Across the four TFs, the benchmark mostly measures the **sequence prior plus
accessibility**, and the reads add little at the level of site rankings:
- **CTCF:** the prior alone reaches 0.66–0.67, against 0.67–0.69 for the posterior. The result
  survives accessibility control (0.63–0.65).
- **KLF4:** beyond accessibility, every HMM is about 0.
- **OCT4:** controlled for accessibility, `neural_global` (0.33–0.40) matches the external
  motif score (0.28–0.41). The reads add about 0.05–0.07 over its prior.

**Context-offset experiment.** We retrained `neural_global` with a learned ±3 bp context
logit-offset table, plus a control retrain with no table.
- Validation NLL improves slightly: 12.337 vs 12.356 per read.
- ChIP agreement changes within noise.
- **The OCT4 weak-TF floor does not drop** (0.212 vs 0.208 mean P(TF) on all-accessible motif
  reads).
- The learned table has no special offset at OCT4 motif calls.

So the context-bias explanation of the floor is not supported.

**MAP-target experiment: the occupancy KL target is the cause.** We retrained `neural_global`
with prior-occupancy KL targets set to the control model's MAP group occupancy. Total TF went
from 0.06 to 0.020; everything else was unchanged.
- **The floor disappears.** The Q4 contribution of all-accessible reads to OCT4 site-mean P(TF)
  falls from 0.063 / 0.092 (chr8 / chr16) to 0.002 / 0.002.
- **Footprints are kept.** Reads protected at the motif keep their P(TF), for OCT4 and for CTCF.
- **Likelihood is slightly better:** validation NLL 12.351 vs 12.356 per read.
- **The weak TF groups leave their floor:** `TF_short_weak` 0.51 and `TF_long_weak` 0.74
  P(protected), previously pinned at 0.2.
- **ChIP agreement:**
  - CTCF is unchanged.
  - OCT4 falls to single-TF's level: raw 0.48 / 0.57 → 0.22 / 0.23; partial given
    accessibility 0.31 / 0.41 → 0.15 / 0.20.
  - ESRRB falls to about 0 (0.28 / 0.30 → 0.00 / 0.01).
  - KLF4 falls from 0.24 / 0.20 to 0.15 / 0.11.

The multi-TF model's OCT4/ESRRB/KLF4 advantage was therefore produced by a KL target about 3.5
times the TF rate the reads support. That target pushed the sequence prior to spread weak-TF
mass over motifs in accessible DNA.

# Key Points

- **Clarification of the earlier OCT4 reading** (user correction, 2026-11-03).
  - At the example site `chr8:24571331-24571346:-`, reads with a footprint (protected at the
    motif, accessible flanks, e.g. `aPPaaaP`, `aaPaaaa`) get `neural_global` motif P(TF)
    0.67–0.82. So the model does detect real OCT4 footprints in individual reads.
  - What drives its site-level ChIP correlation is the 0.42 floor (mostly `TF_short_weak`) on
    the ~46 fully accessible reads at that site.
  - Encoding throughout: 1 = protected, 0 = accessible.
- **Controlled correlations** (chr8 / chr16; Spearman with pooled ChIP-nexus):

| | CTCF | ESRRB | KLF4 | OCT4 |
|---|---|---|---|---|
| **raw** | | | | |
| accessibility (±25 bp of motif) | 0.29 / 0.28 | 0.26 / 0.29 | **0.47 / 0.41** | 0.41 / 0.46 |
| motif score (external PWM) | 0.54 / 0.52 | 0.25 / 0.28 | 0.22 / 0.20 | 0.29 / 0.41 |
| neural_global | 0.67 / 0.69 | 0.25 / 0.30 | 0.24 / 0.21 | 0.50 / 0.57 |
| neural_global, prior only | 0.66 / 0.67 | 0.10 / 0.12 | 0.14 / 0.11 | 0.39 / 0.43 |
| single_tf | 0.64 / 0.67 | 0.27 / 0.30 | 0.13 / 0.11 | 0.34 / 0.41 |
| **partial, given accessibility** | | | | |
| motif score | 0.53 / 0.52 | 0.26 / 0.28 | 0.17 / 0.17 | 0.28 / 0.41 |
| neural_global | 0.63 / 0.65 | 0.13 / 0.17 | 0.02 / 0.04 | 0.33 / 0.40 |
| neural_global, prior only | 0.63 / 0.65 | 0.02 / 0.04 | 0.01 / 0.02 | 0.28 / 0.33 |
| single_tf | 0.62 / 0.65 | 0.16 / 0.18 | −0.01 / 0.01 | 0.14 / 0.22 |

  - Stratifying by accessibility quintile agrees with the partial correlation.
  - Restricting occupancy to reads with a call in the motif barely changes `neural_global`.
- **Context-offset experiment** (`2026-11-03_fixed_blocks_global_kl3_context3` vs
  `_rerun`; chr8 shown, chr16 the same pattern):

| | control | ±3 context table |
|---|---|---|
| validation NLL per read | 12.356 | 12.337 |
| raw Spearman, CTCF / ESRRB / KLF4 / OCT4 | 0.66 / 0.28 / 0.24 / 0.48 | 0.66 / 0.29 / 0.25 / 0.50 |
| partial Spearman given accessibility | 0.63 / 0.16 / 0.02 / 0.31 | 0.63 / 0.17 / 0.02 / 0.33 |
| OCT4 per-read P(TF) on all-accessible motif reads | 0.208 | 0.212 |
| learned offset at OCT4 motif calls vs all calls in those windows | — | −0.26 vs −0.29 |

  - The rerun reproduces the original `neural_global` within about 0.03 on every
    correlation, so code changes since 2026-09-29 do not move the results.
- **Occupancy of the control model** (1,000 val regions):
  - **MAP:** TF_short_strong 0.007, TF_short_weak 0.0001, TF_long_strong 0.0127, TF_long_weak
    0.0005 (total TF 0.020). Open 0.383, other_long 0.102, histone 0.494.
  - **Posterior marginals and prior:** each TF group about 0.017–0.019 (total TF about 0.072),
    matching the KL targets of 0.015 each.
  - So the weak TF groups are diffuse posterior mass that is almost never on the MAP path.

# Details

**New metrics** (`evaluation/tf_chip/analyze_covariation.py`, commit `53b9552`):
- Per site:
  - `accessibility`: the fraction of accessible calls within ±25 bp of the motif, over all
    reads;
  - `n_reads_covered`;
  - `<model>_covered_mean`: mean P(TF) over reads with a call in the motif.
- Per TF and model: `spearman_rho`, `spearman_rho_covered`, `partial_rho_accessibility`
  (Spearman of rank residuals after linear regression on the accessibility rank), and
  `stratified_rho_accessibility` (the site-weighted mean Spearman within 5 accessibility
  quantiles).
- Reference rows: `ref_accessibility`, `ref_motif_score` and `ref_accessibility_motif` (rank
  sum).
- The motif score is not a trivial baseline, since the models never see a PWM (user
  correction), but it is a useful reference.
- Limits:
  - Five strata control a strong driver only partly. In a test with a shared driver,
    stratification lowers but does not remove the correlation; the partial version does.
  - `analysis_2026-11-04` shows NaN for `ref_accessibility` and `ref_accessibility_motif` on
    the raw metric, because sites with no nearby call have undefined accessibility. This was
    fixed by making `spearman` skip missing pairs, and the corrected values are in the table
    above. `analysis_2026-11-04_context3` was run after the fix.

**Prior-only references** (`--prior-only` in `smf_hmm.cli.infer_posterior`, commit `8626d7f`):
the prior state marginals of each window, written for every read in the posterior format so
the analysis scores them like any other model.

**Context-offset model** (commit `d160bfc`):
- `ContextLogitOffsets` in `smf_hmm/emission_modules.py`: a learned `nnx.Param` over all
  4^(2k+1) ±k contexts, zero-initialised, centred over the 1,792 callable contexts (for k = 3).
- `HMM(context_offsets=...)` adds it to every state's logit of emitting 1 in `region_terms`, so
  training, inference and plots all use it. Checkpoints record `{k, assay_type}`.
- Enabled with `emission.context_offset_k=3`. The run otherwise uses `neural_global`'s exact
  overrides; the resolved configs differ only in unused fields.
- Design reasoning:
  - **Why logit-additive.** `context_bias_forms.py` compared logit-additive and efficiency
    forms on held-out chromosomes, using local accessibility as a proxy for the state mix.
    Logit-additive fit slightly better at every width. The gain over local accessibility is
    only 0.0037–0.0039 bits/call at ±3–4, about 12 times smaller than the GM12878 m6A bias.
  - **Why learned rather than fixed.** It needs no external fit, and its lookups stay
    consistent under RC augmentation.
- **What it learned:**
  - offsets with SD 0.50 over callable contexts;
  - a call-weighted correlation of 0.36 with the external ±3 fit (0.49 excluding CGC contexts;
    that fit used merged-CGC data);
  - most protected: CpGs next to G/C, e.g. `AAACGGC` +1.63;
  - least protected: GpCs in `AGCT` or `GGCC`, e.g. `AAGCTTA` −1.36.

**Interpretation (inferred):**
- The ±3 table cannot explain the OCT4 floor. The model learned no special offset at OCT4
  motif calls, and the floor is unchanged.
- **Revised hypothesis:** with KL weight 3, the prior-occupancy KL pulls prior TF occupancy to
  4 × 0.015 = 0.06. Read evidence supports about 0.02, judging by the MAP.
  - The prior must place the extra mass somewhere. The cheapest place is a weak-TF state
    (P(protected) = 0.2, nearly open), at sequences the CNN recognises, such as motifs.
  - A context table cannot relieve that pressure.
  - The same mechanism fits the call-free 8 bp TF segments in the
    [merged-CGC run](../2026-11-03/merged-cgc-147bp-histone-run.md).
- **Caveat:** posterior-marginal occupancy matches the prior (about 0.072 total TF). Weak TF
  is close to open per read, so the reads hardly contradict it on average. The MAP rate is a
  hard-path summary that undercounts diffuse states, so it is a lower bound on what the data
  supports.

**MAP-target experiment** (`2026-11-04_fixed_blocks_global_kl3_map_targets`).

Setup:
- KL targets are the control's MAP occupancy: [open, TF_short_strong, TF_short_weak,
  TF_long_strong, TF_long_weak, other_long, histone] = [0.3833, 0.007, 0.0001, 0.0127, 0.0005,
  0.102, 0.4944]. Targets must sum to 1.
- KL weight 3; everything else identical to the control.
- Trained 100 epochs. Validation NLL 12.351 (control 12.356); prior-occupancy KL 0.007.
- Learned P(protected): TF_short_weak 0.51, TF_long_weak 0.74. Strong groups, other_long and
  histone stay at 0.99; open at 0.01.

Predictions made before the result:
1. the OCT4 floor drops;
2. footprint reads and CTCF keep their P(TF);
3. OCT4's raw Spearman falls, while its partial Spearman given accessibility holds up.

Predictions 1 and 2 held. Prediction 3 held only in part: the partial Spearman fell too. So
the floor carried motif-prior information beyond accessibility as well.

- **Contribution to site-mean motif P(TF) by read category**, low → high ChIP quartile,
  chr8 / chr16:

| | control | MAP targets |
|---|---|---|
| OCT4, all-accessible reads | 0.021→0.063 / 0.021→0.092 | **0.001→0.002 / 0.002→0.002** |
| OCT4, no motif call | 0.023→0.055 / 0.019→0.051 | 0.005→0.008 / 0.005→0.007 |
| OCT4, protected motif call | 0.016→0.030 / 0.016→0.034 | 0.015→0.028 / 0.015→0.031 |
| CTCF, protected motif call | 0.023→0.521 / 0.024→0.525 | 0.021→0.524 / 0.021→0.527 |

- **ChIP agreement**, chr8 / chr16:

| | CTCF | ESRRB | KLF4 | OCT4 |
|---|---|---|---|---|
| **raw Spearman** | | | | |
| control (`neural_global_rerun`) | 0.66 / 0.68 | 0.28 / 0.30 | 0.24 / 0.20 | 0.48 / 0.57 |
| MAP targets | 0.65 / 0.66 | −0.00 / 0.01 | 0.15 / 0.11 | 0.22 / 0.23 |
| MAP targets, prior only | 0.69 / 0.68 | −0.07 / −0.11 | 0.18 / 0.15 | 0.22 / 0.21 |
| single_tf | 0.64 / 0.67 | 0.27 / 0.30 | 0.13 / 0.11 | 0.34 / 0.41 |
| **partial Spearman, given accessibility** | | | | |
| control | 0.63 / 0.65 | 0.16 / 0.16 | 0.02 / 0.04 | 0.31 / 0.41 |
| MAP targets | 0.64 / 0.65 | −0.05 / −0.05 | 0.01 / 0.02 | 0.15 / 0.20 |
| single_tf | 0.62 / 0.65 | 0.16 / 0.18 | −0.01 / 0.01 | 0.14 / 0.22 |
| motif score (reference) | 0.53 / 0.52 | 0.26 / 0.28 | 0.17 / 0.17 | 0.28 / 0.41 |

**Interpretation:**
- With a TF target near what the reads support, P(TF) means footprint evidence: high on reads
  protected at the motif, near zero otherwise.
- The ChIP benchmark then credits only CTCF strongly, plus a residual OCT4 signal (partial
  0.15–0.20). That residual fits the small binding-specific OCT4 footprint (about +0.02–0.03
  over context) measured in the [ChIP evaluation
  note](../2026-11-03/chip-eval-single-tf-vs-multi-tf.md).
- ESRRB and KLF4 have no detectable footprint here. Their earlier correlations came from the
  prior's motif-plus-accessibility placement.
- **The trade-off is a modelling choice.**
  - If P(TF) should be footprint evidence, the MAP-target model is the honest one, and it
    fits the reads at least as well. The benchmark should then be read through partial or
    posterior-minus-prior metrics.
  - If P(TF) should be binding probability, it has to come from sequence plus accessibility.
    That is better modelled explicitly than through an occupancy target.
- **Caveat:** MAP undercounts diffuse states, so the MAP-rate targets are a lower bound. An
  intermediate target, e.g. halfway between the MAP rates and 0.015 per TF group, would show
  how gradual the trade-off is. Not yet run (user asked to hold).

# Commands

All runs used Slurm job 50835752 (akundaje partition, 1 L40S, 4 CPUs, 24 GB). The GPU is in
exclusive-process mode, so prior-only inference for the older checkpoints ran on CPU while
training held the GPU, loading GPU-saved checkpoints with `HostCheckpointer`.

```bash
cd /home/users/diamant/repos/smf_hmm
M=/scratch/users/diamant/mesc_validation/covariation_windows
bash /scratch/users/diamant/2026-11-03_neural_global_pair.sh                          # context3 + rerun training, ~2.3 h each
bash $M/infer_prior_only_2026-11-04.sh                                                # CPU prior-only for neural_global, single_tf, single_tf_147_merged_cgc
bash $M/infer_neural_global_pair_2026-11-04.sh                                        # GPU posteriors + prior-only for context3, rerun
PY=/home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python
for CHR in chr8 chr16; do D=$M/$CHR-512-valid90-reads10; P=$D/posteriors
  $PY -m evaluation.tf_chip.analyze_covariation $D/dataset.h5 \
    --prediction phase_type=$P/phase_type.h5 \
    --prediction neural_global=$P/neural_global.h5 --prediction neural_global_prior=$P/neural_global_prior.h5 \
    --prediction neural_global_rerun=$P/neural_global_rerun.h5 --prediction neural_global_rerun_prior=$P/neural_global_rerun_prior.h5 \
    --prediction neural_global_context3=$P/neural_global_context3.h5 --prediction neural_global_context3_prior=$P/neural_global_context3_prior.h5 \
    --prediction single_tf=$P/single_tf.h5 --prediction single_tf_prior=$P/single_tf_prior.h5 \
    --output-dir $D/analysis_2026-11-04_context3
  $PY -m evaluation.tf_chip.motif_call_diagnostics $D/dataset.h5 $D/analysis_2026-11-04_context3 \
    --model neural_global_rerun --model neural_global_context3
done
bash /scratch/users/diamant/2026-11-04_fixed_blocks_global_kl3_map_targets/train.sh   # MAP-target training, job 50835752
# The allocation ended after training finished; the rest ran in Slurm job 51132065 (sh04-06n09, 1 L40S).
bash $M/infer_map_targets_2026-11-04.sh                                               # posteriors + prior-only, chr8 and chr16
for CHR in chr8 chr16; do D=$M/$CHR-512-valid90-reads10; P=$D/posteriors
  $PY -m evaluation.tf_chip.analyze_covariation $D/dataset.h5 \
    --prediction neural_global=$P/neural_global.h5 --prediction neural_global_prior=$P/neural_global_prior.h5 \
    --prediction neural_global_rerun=$P/neural_global_rerun.h5 --prediction neural_global_rerun_prior=$P/neural_global_rerun_prior.h5 \
    --prediction neural_global_context3=$P/neural_global_context3.h5 \
    --prediction neural_global_map_targets=$P/neural_global_map_targets.h5 --prediction neural_global_map_targets_prior=$P/neural_global_map_targets_prior.h5 \
    --prediction single_tf=$P/single_tf.h5 \
    --output-dir $D/analysis_2026-11-04_map_targets
  $PY -m evaluation.tf_chip.motif_call_diagnostics $D/dataset.h5 $D/analysis_2026-11-04_map_targets \
    --model neural_global_rerun --model neural_global_map_targets
done
```

- One-off snippets, not saved scripts:
  - the learned-table summary and its comparison with `logit_table_k3.npz`;
  - the per-read floor and the offsets at motif calls;
  - the control model's MAP, posterior and prior occupancy.
- `analysis_2026-11-04` has the raw-metric NaN bug for two reference rows; see Details.

# Related Notes

- [ChIP evaluation: single-TF vs multi-TF](../2026-11-03/chip-eval-single-tf-vs-multi-tf.md): The OCT4 weak-TF floor, the context-bias hypothesis tested here, and the reference-score idea implemented here.
- [Merged-CGC assay and 147 bp histone run](../2026-11-03/merged-cgc-147bp-histone-run.md): Call-free TF segments, an earlier sign that the occupancy KL pushes TF mass where no reads contradict it.
- [Fiber-seq HMM absorbs m6A context bias](../2026-10-30/fiber-hmm-enzyme-bias-diagnostic.md): The case where context bias matters a lot (Fiber-seq), unlike this mESC assay.

# Open Questions

- How gradual is the trade-off between the original targets and the MAP-rate targets? An
  intermediate run would answer this.
- Should binding probability (sequence plus accessibility) be modelled explicitly, separately
  from footprint evidence?
- In the MAP-target model, how are MAP and posterior-marginal TF occupancy now related?
- For CTCF, reads add only 0.01–0.02 over the prior in site ranking. Do they matter more for
  per-molecule calls, which ChIP cannot score?
- Should the ChIP benchmark's headline metric be the partial Spearman given accessibility, or
  the posterior-minus-prior gain?

# Sources

- Analyses: `/scratch/users/diamant/mesc_validation/covariation_windows/{chr8,chr16}-512-valid90-reads10/analysis_2026-11-04{,_context3,_map_targets}/`
- Posteriors: `.../posteriors/{neural_global,single_tf,single_tf_147_merged_cgc,neural_global_rerun,neural_global_context3,neural_global_map_targets}_prior.h5`, `neural_global_{rerun,context3,map_targets}.h5`
- Runs: `/scratch/users/diamant/2026-11-03_fixed_blocks_global_kl3_{context3,rerun}`, `/scratch/users/diamant/2026-11-04_fixed_blocks_global_kl3_map_targets`
- smf_hmm `nd/experiments`: `d160bfc` (context offsets), `8626d7f` (prior-only), `4436579` (example plots), `53b9552` (controlled metrics), `61f2e83` (diagnostics)
