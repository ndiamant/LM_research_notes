---
title: Which part of the prior makes the TF nubs? Flattening the prior's TF placement along each window removes 86% of off-motif TF calls with no AUROC or AUPRC change; flattening its open/histone landscape instead raises AUROC by 0.016 — a diagnostic, not a fix
date: 2026-11-16
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /scratch/users/diamant/fiber_eval_work/tf_mechanism/{infer_flattened_prior.py, infer_flattened_prior.sbatch (Slurm 52509610), nub_scorecard.py} -> posteriors gelu_flatTF.h5, gelu_flatLand.h5; scorecard_base.tsv, scorecard_flat.tsv
  - /scratch/users/diamant/fiber_eval_work/analyze_flat_val.sbatch (Slurm 52509613) -> idr_summary_flat_val.tsv, kernel_ref/specificity_flat/
  - Model: gelu_s0 (GELU-head ablation B seed 0), all at prior x1
tags: [smf-hmm, fiberseq, prior, tf-calls, diagnostic, chip-benchmark]
---

# Summary

**Question.** The sequence prior both places TF (its TF potentials along the sequence) and shapes the
open/histone/other_long landscape. Which part produces the off-motif TF nubs, and which part costs
posterior AUROC relative to a constant prior?

**Diagnostic.** gelu_s0 is inferred with one part of its prior potentials replaced by its mean over each
window's positions. Emissions, transitions and the window's overall level are unchanged.
- **flatTF:** the TF microstates are flattened.
- **flatLand:** all other microstates are flattened.

This is a post hoc manipulation used only to attribute effects, not proposed as a fix.

**Results (prior ×1).**
- **Flattening TF placement removes almost all of the off-motif TF excess:**
  - 14.4 → 2.0 off-motif TF bases per read-kb, and TF runs 670 → 102 per 1,000 reads;
  - kernel support of TF runs 59% → 90%;
  - the poorly-labelled vs well-labelled rate ratio 3.3 → 1.4;
  - mean AUROC (0.741 → 0.740) and matched AUPRC (0.488 → 0.483) are essentially unchanged.
  - The prior's TF placement helps CTCF, NRF1, USF2 and CREM (their AUROC drops when it is flattened)
    and hurts EBF1, PAX5, RFX5 and SPI1 (theirs rises).
- **Flattening the landscape raises AUROC by +0.016** (0.741 → 0.757) and matched AUPRC by +0.005, with
  little change in off-motif TF (12.7).
- **Neither reaches the constant-prior run** (AUROC 0.796 at ×1). That run's emissions were trained
  without a sequence prior; the gelu_s0 emissions co-adapted with its prior, so part of the gap sits in
  the emissions (inferred).
- **8 bp runs remain 22%** after flattening TF (from 32%). Flattening keeps each TF microstate's
  window-mean potential, so the duration cliff at TF position 9 survives
  ([duration-shaping note](fiber-tf-nubs-spacer-and-duration-shaping.md)).

# Key Points

| model (prior ×1) | AUROC | accessible-read AUROC | protected-read AUROC | ρ within bound | matched AUPRC | off-motif TF / read-kb |
|---|---|---|---|---|---|---|
| gelu_s0 | 0.741 | 0.649 | 0.765 | 0.126 | 0.488 | 14.4 |
| gelu_s0, TF flattened | 0.740 | 0.643 | 0.761 | 0.174 | 0.483 | 2.0 |
| gelu_s0, landscape flattened | 0.757 | 0.664 | 0.779 | 0.142 | 0.493 | 12.7 |
| collapsed (constant prior) | 0.796 | 0.663 | 0.808 | 0.180 | 0.496 | 6.5 |

- **Scorecard** (every 4th val window; off-motif TF runs = P(TF) > 0.5 runs away from benchmark and strong
  footprinting motifs):

  | model | off-motif TF / read-kb | worst / best efficiency decile | TF runs / 1,000 reads | exactly 8 bp | beside a nucleosome | kernel-supported |
  |---|---|---|---|---|---|---|
  | gelu_s0 ×1 | 18.2 | 3.3 | 670 | 32% | 71% | 59% |
  | gelu_s0 ×0.5 | 7.1 | 2.1 | 189 | 10% | 69% | 91% |
  | gelu_s0, TF flattened | 2.4 | 1.4 | 102 | 22% | 67% | 90% |
  | gelu_s0, landscape flattened | 15.6 | 2.7 | 600 | 33% | 67% | 59% |
  | collapsed ×1 | 7.4 | 0.34 | 260 | 17% | 66% | 83% |
  | ReLU B ×1 | 15.6 | 3.4 | 589 | 33% | 71% | 64% |
  | ReLU B + efficiency ×1 | 11.0 | 2.5 | 440 | 33% | 71% | 72% |

- **Per-TF AUROC** (gelu_s0 / TF flattened / landscape flattened / collapsed):

  | TF | gelu_s0 | TF flattened | landscape flattened | collapsed |
  |---|---|---|---|---|
  | CTCF | 0.862 | 0.834 | 0.877 | 0.879 |
  | NRF1 | 0.832 | 0.813 | 0.833 | 0.897 |
  | USF2 | 0.774 | 0.757 | 0.790 | 0.816 |
  | CREM | 0.770 | 0.742 | 0.796 | 0.805 |
  | SPI1 | 0.721 | 0.747 | 0.744 | 0.776 |
  | PAX5 | 0.655 | 0.669 | 0.688 | 0.749 |
  | YY1 | 0.599 | 0.592 | 0.611 | 0.717 |

  (EBF1 and RFX5 are in idr_summary_flat_val.tsv.)

# Details

- **Implementation.** `infer_flattened_prior.py` wraps the model's `region_terms` so the chosen columns
  of the per-position prior bias are replaced by their window mean, then calls
  `smf_hmm.posterior_inference.infer_posterior` unchanged. Spacer states share the open potential, so
  they flatten with the landscape.

```bash
T=/scratch/users/diamant/fiber_eval_work/tf_mechanism
sbatch $T/infer_flattened_prior.sbatch                                         # 52509610
sbatch --dependency=afterok:52509610 /scratch/users/diamant/fiber_eval_work/analyze_flat_val.sbatch   # 52509613
python $T/nub_scorecard.py $T/scorecard_flat.tsv gelu_flatTF gelu_flatLand
```

# Related Notes

- [TF nub mechanism overview](fiber-tf-nub-mechanism-overview.md): All hypotheses and fix trials.
- [Spacer and duration shaping](fiber-tf-nubs-spacer-and-duration-shaping.md): The surviving 8 bp
  preference.
- [Prior ablation and ReLU-head collapse](../2026-11-15/fiber-prior-ablation-and-relu-head-collapse.md):
  The collapsed run's higher AUROC.

# Open Questions

- Can training keep the useful TF placement (CTCF, NRF1, USF2, CREM) without nubs elsewhere? The tied
  potentials trial removes duration shaping but not placement.

# Sources

- Scripts and outputs in the frontmatter.
