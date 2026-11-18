---
title: TF 15-60 on the variable-length nucleosome — small gains (AUROC +0.003, ρ within bound +0.017) from removing TF's 40 bp cap, but TF almost never uses 41-60 bp; other_long keeps those patches because TF's zero-bias duration prior makes them 4-5× less likely
date: 2026-11-17
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - Run /scratch/users/diamant/2026-11-17_fiber_hist190_tf60_kl300_alphagenome_bf16 (train.sbatch, code/); comparison r3_hist100_190 (variable nucleosome) and r2_full
  - Slurm training 52722574, val inference 52722576, Viterbi 52722581, evaluation 52722585 (fiber_eval_work/analyze_r4_val.sbatch), readouts 52755268
  - Outputs fiber_eval_work/{idr_summary_r4_val.tsv, kernel_ref/specificity_r4/, prior_log_odds/val_r4_priors.tsv}; tf_mechanism/{viterbi_stats_r4.txt, scorecard_r4.tsv, short_ol/r4_hist190_tf60_x1/}
tags: [smf-hmm, fiberseq, tf-duration, other-long, ctcf, prior-log-odds]
---

# Summary

**Change.** TF's maximum goes from 40 to 60 bp; nothing else changes. The base is histone 100-190
([variable nucleosome](fiber-variable-length-nucleosome.md)). The goal was for one TF state to absorb
the isolated 41-60 bp other_long patches, many of them at bound CTCF
([short other_long](fiber-short-other-long-calls.md)). Seed 0, val chromosomes only, prior ×1.

**Results.**
- **Small improvements on most metrics:**
  - AUROC 0.765 → 0.768;
  - accessible-read AUROC 0.687 → 0.688;
  - ρ within bound 0.215 → 0.232 (best so far);
  - matched AUPRC 0.535 → 0.538;
  - off-motif TF 4.78 → 4.69 per read-kb.

  Mean prior TF log-odds is flat (1.09 → 1.08): NRF1 1.81 → 2.05 and YY1 0.51 → 0.59, while CREM
  0.89 → 0.63 and USF2 1.26 → 1.04. CTCF AUROC 0.861 → 0.865; CTCF matched AUPRC 0.823 → 0.827.
- **TF barely uses the new range.**
  - Only 0.7% of Viterbi TF calls exceed 40 bp (1.0 per 1,000 reads).
  - TF calls at the old 40 bp cap fall from 2.5% to 0.2%, so the gains probably come from removing the
    cap's censoring (inferred).
  - other_long is nearly unchanged: 126 vs 129 segments per 1,000 reads; 1,329 vs 1,369 short calls,
    93% of them still 41-60 bp patches.
- **Why other_long keeps them** (inferred from the model's structure):
  - **Emissions are not the reason:** both runs learned TF 0.94 and other_long 0.956 protection.
  - **TF's zero-bias duration is a truncated geometric decaying from 15 bp.** At 41-59 bp it is
    4.1-5.1× less likely than other_long's (e.g. 50 bp: 0.0063 vs 0.0296). After per-position entry
    rates (TF 2.3e-4 vs other_long 1.7e-4), other_long is still about 3-4× more likely a priori
    (about 1.3 nats).
  - **TF must also be flanked by 10 bp open spacers.** 36% of the remaining ≤ 60 bp other_long calls sit
    within 12 bp of a nucleosome.

**Interpretation.** Raising the cap alone cannot move these patches into TF. TF would need a flatter
duration prior over 15-60 bp, e.g. the `duration_mode` / `duration_scale` form, or a uniform one.
Until then, TF 15-60 is a modest improvement from removing the 40 bp censoring. It is a reasonable
new base if its prior changes (CREM, USF2) are acceptable.

# Key Points

| model | AUROC | accessible-read AUROC | protected-read AUROC | ρ within bound | matched AUPRC | on-motif share | off-motif TF / read-kb | prior TF log-odds |
|---|---|---|---|---|---|---|---|---|
| r2_full (histone 147, TF 15-40) | 0.764 | 0.681 | 0.780 | 0.213 | 0.528 | 0.353 | 5.65 | 1.01 |
| histone 100-190, TF 15-40 | 0.765 | 0.687 | 0.781 | 0.215 | 0.535 | 0.355 | 4.78 | 1.09 |
| **histone 100-190, TF 15-60** | **0.768** | **0.688** | **0.783** | **0.232** | **0.538** | **0.362** | **4.69** | 1.08 |

| model | TF / 1,000 reads | TF at 15 bp | TF median | TF > 40 bp | TF at its max | other_long / 1,000 reads | CTCF / NRF1 median at bound motifs | off-motif TF runs / 1,000 reads | kernel-supported |
|---|---|---|---|---|---|---|---|---|---|
| histone 100-190, TF 15-40 | 152 | 21% | 18 | – | 2.5% (40) | 130 | 27 / 21 | 116 | 86% |
| histone 100-190, TF 15-60 | 152 | 17% | 19 | 0.7% | 0.2% (60) | 126 | 29 / 21 | 107 | 86% |

- **Per-TF AUROC:**

  | TF | r2_full | 100-190 / TF 40 | 100-190 / TF 60 |
  |---|---|---|---|
  | CTCF | 0.868 | 0.861 | 0.865 |
  | NRF1 | 0.868 | 0.859 | 0.866 |
  | CREM | 0.825 | 0.808 | 0.803 |
  | USF2 | 0.788 | 0.783 | 0.782 |
  | RFX5 | 0.753 | 0.771 | 0.771 |
  | SPI1 | 0.746 | 0.753 | 0.755 |
  | EBF1 | 0.731 | 0.729 | 0.727 |
  | PAX5 | 0.697 | 0.702 | 0.706 |
  | YY1 | 0.599 | 0.617 | 0.633 |

- **Prior TF log-odds:**

  | TF | r2_full | 100-190 / TF 40 | 100-190 / TF 60 |
  |---|---|---|---|
  | CTCF | 2.56 | 2.70 | 2.77 |
  | NRF1 | 1.68 | 1.81 | 2.05 |
  | USF2 | 1.08 | 1.26 | 1.04 |
  | CREM | 0.87 | 0.89 | 0.63 |
  | RFX5 | 0.75 | 0.85 | 0.86 |
  | SPI1 | 0.66 | 0.73 | 0.67 |
  | EBF1 | 0.54 | 0.56 | 0.58 |
  | PAX5 | 0.57 | 0.51 | 0.51 |
  | YY1 | 0.37 | 0.51 | 0.59 |
  | mean | 1.01 | 1.09 | 1.08 |

- **Zero-bias duration probability** (`_duration_exit_probabilities`, 5% at the cap):

  | length | TF 15-60 | other_long 40-100 | ratio |
  |---|---|---|---|
  | 41 bp | 0.0114 | 0.0463 | 4.1 |
  | 50 bp | 0.0063 | 0.0296 | 4.7 |
  | 55 bp | 0.0045 | 0.0230 | 5.1 |

  Mean TF duration is 28.8 bp; other_long 58.6 bp.

# Details

- **Learned emissions** (wandb summaries):
  - TF 0.941 (TF 40) and 0.937 (TF 60);
  - other_long 0.956 in both;
  - histone 0.987 in both;
  - open 0.44 in both.
- **Remaining other_long.** Of calls ≤ 60 bp, 70% have an m6A inside, and the median protected fraction
  is 0.95.
- **No smoke run.** The code was identical to the 100-190 run, which had a smoke test.

```bash
S=/scratch/users/diamant; W=$S/fiber_eval_work; T=$W/tf_mechanism
sbatch $S/2026-11-17_fiber_hist190_tf60_kl300_alphagenome_bf16/train.sbatch   # TF {min 15, max 60}
sbatch --dependency=afterok:<train> $W/infer_variant_val_extra_any.sbatch r4_hist190_tf60 $S/2026-11-17_fiber_hist190_tf60_kl300_alphagenome_bf16
MODELS=r4_hist190_tf60 sbatch --dependency=afterok:<train> $T/viterbi_lengths.sbatch
sbatch --dependency=afterok:<infer> $W/analyze_r4_val.sbatch
python $T/viterbi_stats.py r2_full_x1 r3_hist100_190_x1 r4_hist190_tf60_x1
CONFIG=r4_hist190_tf60_x1 POSTERIOR=r4_hist190_tf60 N_READS=0 python $T/short_other_long.py
python $T/nub_scorecard.py $T/scorecard_r4.tsv r4_hist190_tf60
```

# Related Notes

- [Variable-length nucleosome](fiber-variable-length-nucleosome.md): The base for this run.
- [Short other_long calls](fiber-short-other-long-calls.md): The 41-60 bp patches this run targeted.
- [Spacer and duration shaping](../2026-11-16/fiber-tf-nubs-spacer-and-duration-shaping.md): Earlier
  evidence that TF's duration shaping drives where calls land.

# Open Questions

- Would a flatter TF duration prior over 15-60 bp move the patches into TF, and would it bring back
  long TF chains or nubs?
- Why does CREM's prior log-odds fall (0.89 → 0.63) when TF can be longer?

# Sources

- Run, jobs and outputs in the frontmatter.
