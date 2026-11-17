---
title: Ablation round 1 from the KL 300 GELU model — TF 15-40 is the best single change and tying TF potentials helps; a wider histone flank is mildly positive; narrowing other_long to 45-80 and widening the open band do not help
date: 2026-11-16
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - Runs /scratch/users/diamant/2026-11-16_fiber_abl2_{tie,tf15,flank,ol45_80,flank_ol45_80,open065}_kl300_alphagenome_bf16 (train.sbatch, code_state.patch); base /scratch/users/diamant/2026-11-16_fiber_gelu_B_kl300_alphagenome_bf16
  - Slurm training 52574423/39/43/45/47/49, inference 52574431/40/44/46/48/50, evaluation 52574451 (fiber_eval_work/analyze_r1_val.sbatch), Viterbi 52574464-69
  - Outputs fiber_eval_work/{idr_summary_r1_val.tsv, kernel_ref/specificity_r1/, prior_log_odds/val_r1_priors.tsv}; tf_mechanism/{scorecard_r1.tsv, viterbi/viterbi_segments_r1_*_x1.tsv.gz, open_overdispersion_*_x1.tsv, open_protection_r1.txt}
tags: [smf-hmm, fiberseq, ablation, tf-duration, histone, other-long, open-state, prior-kl]
---

# Summary

**Design.** The base is GELU design B with prior occupancy KL 300 (`fix_kl300`): spacer + TF/3, no TF
chaining, other_long 40-100 bp, 147 bp histone with a 25 bp taper, open band [0.05, 0.45]. Each run
makes one change, except A3 + A4 together. Seed 0, val only, prior ×1.
- **A1 tie:** TF and other_long chain states share one prior potential each.
- **A2 TF length:** TF 15-40 bp.
- **A3 flank:** histone taper 40 bp with edge floor 0.2.
- **A4 other_long:** 45-80 bp, occupancy target 0.5%.
- **A3 + A4.**
- **A6 open:** open band [0.05, 0.65]. This requires closed band [0.65, 0.995], since open's upper bound
  must not exceed closed's lower; the histone edge floor is kept at 0.45.

**Results.** Effects are small; with GELU seed spread about 0.001 AUROC and 0.02 log-odds, differences of
≥ 0.005 AUROC are probably real.
- **A2 (TF 15-40) is the best single change:**
  - best AUROC (0.756 vs 0.744), accessible-read AUROC (0.674) and ρ within bound sites (0.227 vs
    0.208);
  - AUPRC 0.527;
  - 44% fewer TF runs, 92% kernel-supported (vs 82%);
  - CREM AUROC 0.740 → 0.782, NRF1 0.840 → 0.862.
  - The minimum-length pile-up moves with the minimum (57% of Viterbi TFs exactly 15 bp): the model fills
    short protected patches with the shortest TF it is allowed.
- **A1 (tie) is clearly positive:**
  - the share of Viterbi TFs at the minimum falls 40% → 19% (median 9 → 11 bp);
  - best prior TF log-odds, mean 0.83 → 0.94 (CTCF 1.97 → 2.54);
  - AUROC +0.007, kernel support 86%.
- **A3 (flank) is mildly positive:** AUROC +0.007, off-motif TF 6.18 → 5.84, fewer other_long at the
  100 bp cap (58% → 49%), ρ −0.01.
- **A4 (other_long 45-80) hurts:** AUPRC −0.008, off-motif TF +0.5. other_long re-piles at its new 80 bp
  cap (56-60%), so 80-146 bp stretches still lack a home, even with the wider flank (A3 + A4: AUROC
  0.755 but AUPRC 0.518, ρ 0.189).
- **A6 (open to 0.65) does not help.** Open learns 0.488, not the new ceiling (all other runs pin at
  0.45). Off-motif TF is best (5.68), but AUROC, AUPRC and prior log-odds drop (NRF1 prior 1.71 → 1.29).
- **No change models the clustered protection in open DNA.** Runs of ≥ 8 protected calls are 1.26-1.40×
  the open emission's expectation in every run (base 1.31). Open is also over-predicted as protected at
  ×1: 0.40 observed vs 0.43.

**Caveat.** With A1 or A2, NRF1's median Viterbi TF length at bound motifs drops from 25 to 16-17 bp
while its AUROC rises.

**Proposed round 2:** full = A1 + A2 + A3 (other_long 40-100 kept, open unchanged), plus leave-one-out of
each.

# Key Points

| run | AUROC | accessible-read AUROC | ρ within bound | matched AUPRC | on-motif share | off-motif TF / read-kb | prior TF log-odds (mean) |
|---|---|---|---|---|---|---|---|
| base KL 300 | 0.744 | 0.666 | 0.208 | 0.525 | 0.331 | 6.18 | 0.83 |
| A1 tie | 0.751 | 0.670 | 0.202 | 0.525 | 0.337 | 6.19 | 0.94 |
| A2 TF 15-40 | 0.756 | 0.674 | 0.227 | 0.527 | 0.343 | 6.06 | 0.90 |
| A3 flank | 0.751 | 0.668 | 0.198 | 0.526 | 0.336 | 5.84 | 0.85 |
| A4 OL 45-80 | 0.749 | 0.667 | 0.200 | 0.517 | 0.340 | 6.66 | 0.84 |
| A3 + A4 | 0.755 | 0.671 | 0.189 | 0.518 | 0.348 | 6.31 | 0.88 |
| A6 open 0.65 | 0.739 | 0.664 | 0.204 | 0.521 | 0.344 | 5.68 | 0.79 |

| run | TF runs / 1,000 reads | Viterbi TFs at the minimum | kernel-supported | Viterbi TF median | CTCF / NRF1 median at bound motifs | OL at its cap | open overdispersion ≥ 8 / ≥ 12 |
|---|---|---|---|---|---|---|---|
| base | 269 | 40% (8 bp) | 82% | 9 | 30 / 25 | 58% (100) | 1.31 / 1.52 |
| A1 | 234 | 19% | 86% | 11 | 29 / 17 | 47% | 1.28 / 1.43 |
| A2 | 151 | 57% (15 bp) | 92% | 15 | 31 / 16 | 59% | 1.40 / 1.40 |
| A3 | 252 | 39% | 80% | 9 | 29 / 24 | 49% | 1.27 / 1.47 |
| A4 | 264 | 36% | 84% | 10 | 31 / 22 | 60% (80) | 1.30 / 1.46 |
| A3 + A4 | 245 | 36% | 83% | 10 | 31 / 15 | 56% (80) | 1.29 / 1.45 |
| A6 | 245 | 37% | 84% | 9 | 30 / 21 | 59% | 1.26 / 1.43 |

- **Per-TF AUROC** for the TFs with the largest changes:

  | TF | base | A1 tie | A2 TF 15-40 | A3 flank | A4 OL 45-80 | A3 + A4 | A6 open |
  |---|---|---|---|---|---|---|---|
  | CREM | 0.740 | 0.766 | 0.782 | 0.748 | 0.764 | 0.769 | 0.731 |
  | NRF1 | 0.840 | 0.848 | 0.862 | 0.841 | 0.835 | 0.848 | 0.842 |
  | RFX5 | 0.738 | 0.743 | 0.757 | 0.753 | 0.756 | 0.749 | 0.738 |
  | YY1 | 0.590 | 0.606 | 0.607 | 0.595 | 0.603 | 0.612 | 0.581 |

# Details

- **Smoke tests.** 1-epoch smoke runs passed for A6 and A2 before submission.
- **Overdispersion per run.** It uses that run's Viterbi open segments, learned open protection (wandb
  summary) and context offsets.
- **Readouts.** The scorecard and Viterbi metrics come from `nub_scorecard.py` (every 4th window) and
  `viterbi_lengths.py` (every 8th window).

```bash
S=/scratch/users/diamant; W=$S/fiber_eval_work
for v in tie tf15 flank ol45_80 flank_ol45_80 open065; do sbatch $S/2026-11-16_fiber_abl2_${v}_kl300_alphagenome_bf16/train.sbatch; done
# + infer_variant_val_extra_any.sbatch r1_<v> ...; analyze_r1_val.sbatch; MODELS=r1_<v> viterbi_lengths.sbatch
python $W/tf_mechanism/nub_scorecard.py $W/tf_mechanism/scorecard_r1.tsv r1_tie r1_tf15 r1_flank r1_ol45_80 r1_flank_ol45_80 r1_open065
CONFIG=r1_<v>_x1 RUN=<run dir> P_OPEN=<learned> python $W/tf_mechanism/open_overdispersion.py
```

# Related Notes

- [TF nub mechanism overview](fiber-tf-nub-mechanism-overview.md): Why these changes were tested.
- [Short TF calls: sparsity vs dissociation](fiber-short-tf-calls-sparsity-vs-dissociation.md): Why a 15 bp
  minimum was expected to cost little.
- [Open-state overdispersion](fiber-open-state-overdispersion.md): The clustered protection none of these
  changes models.

# Open Questions

- Do A1 + A2 + A3 combine additively (round 2)?
- What should hold the 80-146 bp protected stretches? Other_long re-piles at whatever cap it has.
- Why does NRF1's bound-motif TF length shrink under tying or a 15 bp minimum while its AUROC improves?

# Sources

- Runs, jobs and outputs in the frontmatter.
