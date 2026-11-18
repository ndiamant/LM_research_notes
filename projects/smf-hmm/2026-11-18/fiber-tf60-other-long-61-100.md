---
title: TF 15-60 with other_long 61-100 — TF now takes the 41-60 bp patches (59% on bound motifs, mostly CTCF); best AUROC (0.774), ρ within bound (0.246) and prior log-odds (1.19), but matched AUPRC falls (0.538 → 0.520, CREM and YY1) and off-motif TF rises 31%; other_long now piles at 61 bp
date: 2026-11-18
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - Run /scratch/users/diamant/2026-11-17_fiber_hist190_tf60_ol61_kl300_alphagenome_bf16 (train.sbatch, code/); comparison r4_hist190_tf60 (TF 15-60, other_long 40-100)
  - Slurm training 52756811, val inference 52756812, Viterbi 52756813, evaluation 52756814 (fiber_eval_work/analyze_r5_val.sbatch), readouts 52855408
  - Outputs fiber_eval_work/{idr_summary_r5_val.tsv, kernel_ref/specificity_r5/, prior_log_odds/val_r5_priors.tsv}; tf_mechanism/{viterbi_stats_r5.txt, scorecard_r5.tsv, short_ol/r5_tf60_ol61_x1/}
tags: [smf-hmm, fiberseq, tf-duration, other-long, ctcf, prior-log-odds]
---

# Summary

**Change.** other_long's minimum goes from 40 to 61 bp, on top of TF 15-60 and histone 100-190
([TF 15-60](../2026-11-17/fiber-tf-15-60.md)). That run showed that other_long's zero-bias duration
out-weighed TF's at 41-60 bp, so TF barely used the range. Now TF is the only state for 41-60 bp
protected patches. Seed 0, val chromosomes only, prior ×1.

**Results.**
- **TF takes the 41-60 bp patches, and they are good calls.**
  - 9.1% of Viterbi TF calls exceed 40 bp (15.4 per 1,000 reads), up from 0.7%.
  - Of these, 59% sit on IDR-bound motifs, 49% on bound CTCF alone; TF calls of 15-40 bp manage 18%.
  - The median TF length at bound CTCF motifs goes from 29 to 35 bp.
- **Best so far on ranking and the prior:**
  - AUROC 0.774 (TF 15-60: 0.768);
  - protected-read AUROC 0.789;
  - ρ within bound 0.246 (0.232);
  - mean prior TF log-odds 1.19 (1.08), with CTCF 3.07 (2.77) and YY1 0.65.
  - CTCF AUROC 0.880 (0.865), CTCF matched AUPRC 0.849 (0.827); PAX5 AUROC +0.025.
- **Worse on specificity:**
  - **Matched AUPRC falls from 0.538 to 0.520,** mostly from CREM (0.668 → 0.566) and YY1
    (0.325 → 0.287). CREM AUROC is also lower (0.803 → 0.793), although its prior log-odds recovers
    (0.63 → 0.87).
  - **Off-motif TF rises from 4.69 to 6.12 per read-kb** (scorecard version 5.87 → 7.02), and TF
    occupies 1.07% of positions vs 0.82%.
    - Part of the rise may be a measurement effect (inferred): the metric counts TF bases more than
      5 bp from a motif, so a 50 bp TF over a 19 bp CTCF motif adds about 20 off-motif bases.
    - The scorecard's off-motif TF runs barely change (107 → 112 per 1,000 reads), and kernel support
      improves (86% → 90%).
- **other_long piles at its new minimum.**
  - 12% of other_long segments are exactly 61 bp (7.7% at 40 bp before).
  - Calls of 61-63 bp have protected spans of 64-72 bp (IQR); 84% contain an m6A, and 43% sit next to
    a nucleosome.
  - other_long falls to 105 segments per 1,000 reads (126).

**Interpretation.**
- Separating TF and other_long by length works as intended for CTCF-like footprints. It gives the largest AUROC and prior gain of the
  nucleosome and TF-length runs.
- The CREM and YY1 AUPRC losses and the off-motif rise are real costs. Whether they come from longer TF
  calls at weak-footprint TFs, or from the metric counting motif flanks, is untested.
- other_long's new 61 bp pile-up repeats the pattern seen in every range tried: it collects whatever
  its neighbours cannot hold, here 62-72 bp patches.

# Key Points

| model | AUROC | accessible-read AUROC | protected-read AUROC | ρ within bound | matched AUPRC | on-motif share | off-motif TF / read-kb | prior TF log-odds |
|---|---|---|---|---|---|---|---|---|
| histone 100-190, TF 15-40 | 0.765 | 0.687 | 0.781 | 0.215 | 0.535 | 0.355 | 4.78 | 1.09 |
| + TF 15-60 | 0.768 | 0.688 | 0.783 | 0.232 | **0.538** | 0.362 | **4.69** | 1.08 |
| **+ other_long 61-100** | **0.774** | **0.690** | **0.789** | **0.246** | 0.520 | **0.377** | 6.12 | **1.19** |

| model | TF / 1,000 reads | TF > 40 bp | TF median | TF at 15 bp | other_long / 1,000 reads | other_long at its min / max | other_long median | CTCF / NRF1 median at bound motifs | kernel-supported |
|---|---|---|---|---|---|---|---|---|---|
| TF 15-60, other_long 40-100 | 152 | 0.7% | 19 | 17% | 126 | 7.7% (40) / 9.1% | 56 | 29 / 21 | 86% |
| TF 15-60, other_long 61-100 | 171 | 9.1% | 22 | 13% | 105 | 12.2% (61) / 7.9% | 70 | 35 / 16 | 90% |

- **Per TF:**

  | TF | AUROC, TF 15-60 | AUROC, + other_long 61-100 | matched AUPRC, TF 15-60 | matched AUPRC, + other_long 61-100 | prior log-odds, TF 15-60 | prior log-odds, + other_long 61-100 |
  |---|---|---|---|---|---|---|
  | CTCF | 0.865 | 0.880 | 0.827 | 0.849 | 2.77 | 3.07 |
  | NRF1 | 0.866 | 0.877 | 0.791 | 0.785 | 2.05 | 1.95 |
  | USF2 | 0.782 | 0.785 | 0.632 | 0.622 | 1.04 | 1.28 |
  | CREM | 0.803 | 0.793 | 0.668 | 0.566 | 0.63 | 0.87 |
  | RFX5 | 0.771 | 0.766 | 0.473 | 0.463 | 0.86 | 0.96 |
  | SPI1 | 0.755 | 0.763 | 0.413 | 0.403 | 0.67 | 0.77 |
  | EBF1 | 0.727 | 0.737 | 0.388 | 0.392 | 0.58 | 0.54 |
  | PAX5 | 0.706 | 0.731 | 0.326 | 0.315 | 0.51 | 0.62 |
  | YY1 | 0.633 | 0.637 | 0.325 | 0.287 | 0.59 | 0.65 |

- **Where the new long TF calls land.** Viterbi TF calls by overlap with IDR-labelled motifs:

  | TF calls | n | on a bound motif | on an unbound motif | off motifs |
  |---|---|---|---|---|
  | 15-40 bp | 10,004 | 18% | 6% | 77% |
  | 41-60 bp | 1,000 | 59% (49% CTCF) | 5% | 36% |

# Details

- **Learned emissions** were not checked for this run; the TF 15-60 run had TF 0.937 and other_long
  0.956.
- **No smoke run.** The code was identical to the earlier runs.
- **Visual check.** Per-read montages, comparing this run with TF 15-60 / other_long 40-100: 4 reads at
  the strongest IDR-bound site per TF and val chromosome, for 8 TFs. They are in
  `/scratch/users/diamant/gm12878_validation/ol61_vs_ol40_panels/montages/`, made with
  `evaluation.tf_chip.plot_model_comparison`; reads were chosen from the new model.
  - **Short patches move from other_long to TF everywhere.** Protected patches of about 35-60 bp that
    the old model called other_long (posterior about 0.6-0.9, TF 0.2-0.4) are now called TF with
    posterior about 0.9-1.0.
  - **This happens at motifs and away from them.** Examples:
    - on-motif: the CREM motif in chr20:57,398,498 read 8; the YY1 motif in chr8:22,245,102 read 96;
    - off-motif: CREM read 23 (bp ~430-465) and read 11 (bp ~150-175); YY1 read 65 (bp ~380-420);
      CTCF chr8:22,864,961 read 104 (bp ~100-130); NRF1 chr8:10,054,697 reads 97-98.
  - **Length, not sequence, decides.** The prior's TF potential is near zero at these off-motif
    patches.
  - **So the off-motif TF rise is not only the motif-flank counting effect.** Off-motif 41-60 bp
    patches become TF too, which likely explains the CREM and YY1 AUPRC loss (inferred from examples,
    not quantified).

```bash
S=/scratch/users/diamant; W=$S/fiber_eval_work; T=$W/tf_mechanism
sbatch $S/2026-11-17_fiber_hist190_tf60_ol61_kl300_alphagenome_bf16/train.sbatch   # other_long {min 61, max 100}
sbatch --dependency=afterok:<train> $W/infer_variant_val_extra_any.sbatch r5_tf60_ol61 $S/2026-11-17_fiber_hist190_tf60_ol61_kl300_alphagenome_bf16
MODELS=r5_tf60_ol61 sbatch --dependency=afterok:<train> $T/viterbi_lengths.sbatch
sbatch --dependency=afterok:<infer> $W/analyze_r5_val.sbatch
python $T/viterbi_stats.py r3_hist100_190_x1 r4_hist190_tf60_x1 r5_tf60_ol61_x1
CONFIG=r5_tf60_ol61_x1 POSTERIOR=r5_tf60_ol61 N_READS=0 SHORT_MAX=63 python $T/short_other_long.py
python $T/nub_scorecard.py $T/scorecard_r5.tsv r5_tf60_ol61
```

# Related Notes

- [TF 15-60](../2026-11-17/fiber-tf-15-60.md): Why raising TF's cap alone did not move the patches.
- [Short other_long calls](../2026-11-17/fiber-short-other-long-calls.md): The 41-60 bp patches, many at
  CTCF.
- [Fiber-seq model report](../2026-11-17/fiberseq-model-report.md): The running summary of all
  changes.

# Open Questions

- What are the new off-motif TF bases: flanks of long on-motif calls, or new off-motif calls?
  Measuring off-motif TF runs per call length would tell.
- Why do CREM and YY1 lose matched AUPRC while their prior log-odds rise?
- What should hold the 62-72 bp patches now piled at other_long's minimum?

# Sources

- Run, jobs and outputs in the frontmatter.
