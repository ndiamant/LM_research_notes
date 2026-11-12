---
title: Train-chromosome copy of the IDR benchmark — no memorization in neural posteriors (train − val gap ≈ baselines'), a small train advantage for the prior alone, and 2× tighter CIs confirming phase-type's lead
date: 2026-11-12
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - smf_hmm repo, branch nd/experiments (uncommitted) — evaluation/tf_chip/configs/gm12878_train.yaml
  - Train benchmark /scratch/users/diamant/gm12878_train/ (benchmark_cache/, covariation_windows/{chr2,chr11,chr12,chr17,chr19}-512-valid90-reads10-flank128/ with posteriors/ and analysis_idr/)
  - Pooled summaries /scratch/users/diamant/fiber_eval_work/idr_summary_val_train_ptref.tsv (reference phase-type) and idr_summary_val_train.tsv (reference geometric)
  - Job scripts /scratch/users/diamant/fiber_eval_work/{train_prepare.sbatch, train_baselines.sbatch, train_neural.sbatch, analyze_idr_train.sbatch, summary_val_train_pt.sbatch}
  - Slurm 52032085 (prep), 52032109 / 52032113 (geometric / phase-type), 52032114 (neural), 52038652 (analysis), 52057480 (summary)
tags: [smf-hmm, fiberseq, gm12878, chip, idr, benchmark, memorization, train-split, alphagenome, phase-type]
---

# Summary

**Setup.** We built the [IDR-labelled benchmark](fiber-idr-chip-benchmark.md) again on five of
the models' training chromosomes (chr2, chr11, chr12, chr17, chr19): 18.5k windows against 4.5k on
val, and about 4× the bound sites (CTCF 3,013, NRF1 1,305, USF2 430, CREM 460, RFX5 124). The test
chromosomes (chr1, chr3, chr6) were not used.

**Memorization check: none visible in the posteriors.**
- Every model is about 0.01 higher in AUROC on train than on val, and that includes the
  sequence-free baselines (phase-type +0.009, geometric +0.011 raw AUROC). The neural posteriors
  are not higher (AG KL10 +0.012, AG ol001 +0.010, CNN KL10 +0.010). The shared shift is a
  chromosome difference, not fitting.
- On the four footprinting TFs, all posteriors move by −0.009 to +0.002.
- The **prior alone** (DNA only) shows a small train advantage on footprinting TFs: +0.019 raw /
  +0.025 protected AUROC for AG KL10, +0.010 / +0.014 for AG ol001, against −0.009 to −0.002 for
  the baselines. That is consistent with mild locus memorization in the DNA prior that the reads
  wash out in the posterior (inferred; the size is near per-TF noise).

**So the train chromosomes are usable for power on posterior metrics.** Report them beside val,
and treat prior-only numbers on train with care.

**With about 2× tighter CIs**, the val conclusions hold and sharpen (scale-1 posteriors; the
prior-scale variants were run on val only):
- Phase-type is significantly above every neural posterior on protected-read AUROC for every TF.
- With accessibility controlled, only CNN KL10 ties phase-type on CTCF (0.80).
- Within bound sites, the neural posteriors tie phase-type on CTCF, NRF1 and USF2. They are
  significantly below on CREM, EBF1, PAX5, SPI1, YY1 and RFX5; on val, CREM and EBF1 had been
  ties for some models.
- The AG priors rank bound CTCF sites significantly better than phase-type (0.29 vs 0.27).

# Key Points

- **Train − val gap, mean over the 9 TFs (footprinting TFs CTCF, NRF1, USF2, CREM in brackets):**

  | model | raw AUROC | AUROC within accessibility quintiles | protected-read AUROC |
  |---|---|---|---|
  | phase-type | +0.009 (−0.009) | +0.006 (−0.005) | +0.011 (−0.007) |
  | geometric | +0.011 (−0.002) | +0.008 (+0.006) | +0.012 (−0.003) |
  | AG KL10 | +0.012 (+0.000) | +0.010 (−0.000) | +0.015 (−0.002) |
  | AG ol001 | +0.010 (−0.006) | +0.008 (−0.008) | +0.014 (−0.005) |
  | CNN KL10 | +0.010 (+0.002) | +0.009 (+0.004) | +0.013 (+0.000) |
  | AG KL10 prior only | +0.009 (+0.019) | +0.006 (+0.014) | +0.012 (+0.025) |
  | AG ol001 prior only | +0.010 (+0.010) | +0.009 (+0.005) | +0.012 (+0.014) |

  CNN KL10's prior-only marginals were run on both splits, but the val `analysis_idr` predates
  them, so CNN has no prior-only gap here.
- **The within-bound ChIP rank is not usable as a memorization check.** The train − val
  differences are +0.05 to +0.09 for every model, mostly from RFX5 (val −0.30 on 23 bound sites,
  train +0.20 on 124). Its CIs on val are 0.15-0.76 wide.
- **CI widths (train / val, protected-read AUROC):** CTCF 0.017 / 0.035, NRF1 0.025 / 0.058, USF2
  0.046 / 0.109, CREM 0.052 / 0.100, RFX5 0.084 / 0.184.

# Details

## Train results ("−" / "+" significantly below / above phase-type, 95% paired bootstrap)

| TF | metric | phase-type | geometric | AG KL10 | AG ol001 | CNN KL10 | AG KL10 prior | AG ol001 prior | CNN prior |
|---|---|---|---|---|---|---|---|---|---|
| CTCF | AUROC within accessibility quintiles | 0.80 | 0.73− | 0.76− | 0.77− | 0.80 | 0.78− | 0.78− | 0.76− |
| CTCF | protected-read AUROC | 0.90 | 0.85− | 0.86− | 0.86− | 0.88− | 0.83− | 0.84− | 0.81− |
| CTCF | within-bound ρ | 0.27 | 0.21− | 0.25− | 0.25 | 0.28 | 0.29+ | 0.29+ | 0.26 |
| NRF1 | protected-read AUROC | 0.92 | 0.91− | 0.89− | 0.89− | 0.90− | 0.82− | 0.83− | 0.82− |
| NRF1 | within-bound ρ | 0.40 | 0.33− | 0.40 | 0.41 | 0.42 | 0.40 | 0.39 | 0.41 |
| USF2 | AUROC within accessibility quintiles | 0.70 | 0.64− | 0.64− | 0.66− | 0.66− | 0.69 | 0.70 | 0.65− |
| USF2 | protected-read AUROC | 0.89 | 0.87− | 0.82− | 0.83− | 0.84− | 0.76− | 0.78− | 0.75− |
| USF2 | within-bound ρ | 0.29 | 0.24 | 0.31 | 0.32 | 0.30 | 0.30 | 0.31 | 0.24 |
| CREM | protected-read AUROC | 0.83 | 0.82 | 0.75− | 0.76− | 0.77− | 0.72− | 0.72− | 0.74− |
| CREM | within-bound ρ | 0.43 | 0.33− | 0.33− | 0.34− | 0.34− | 0.31− | 0.36− | 0.29− |
| EBF1 | within-bound ρ | 0.25 | 0.14− | 0.19− | 0.21− | 0.16− | 0.09− | 0.13− | −0.03− |
| SPI1 | AUROC within accessibility quintiles | 0.67 | 0.67 | 0.59− | 0.60− | 0.62− | 0.64− | 0.66− | 0.66 |
| SPI1 | within-bound ρ | 0.24 | 0.21− | 0.11− | 0.13− | 0.18− | 0.19− | 0.18− | 0.21 |
| YY1 | protected-read AUROC | 0.83 | 0.82− | 0.67− | 0.68− | 0.71− | 0.57− | 0.61− | 0.61− |
| PAX5 | within-bound ρ | 0.06 | 0.06 | 0.00− | −0.00− | −0.00− | −0.04− | −0.08− | −0.06− |

All TFs, raw AUROC and CIs are in `idr_summary_val_train_ptref.tsv`.

**Differences from val:**
- CREM's accessibility-controlled AUROC is lower on train for every model (phase-type 0.59 vs
  0.65).
- YY1 is higher (0.67 vs 0.58) and RFX5's within-bound ρ is positive (0.20 vs −0.30).
- These shifts are shared by the baselines, so they reflect chromosome composition or the small
  val counts, not the models.

## Command provenance

```bash
cd /home/users/diamant/repos/smf_hmm; W=/scratch/users/diamant/fiber_eval_work
sbatch $W/train_prepare.sbatch                                              # 52032085: calls TSVs, cache, window selection
sbatch --dependency=aftercorr:52032085 --job-name=geometric  $W/train_baselines.sbatch geometric   # 52032109
sbatch --dependency=aftercorr:52032085 --job-name=phase_type $W/train_baselines.sbatch phase_type  # 52032113 (chr2 4.6 h)
sbatch --dependency=aftercorr:52032085 $W/train_neural.sbatch               # 52032114 (AG KL10, AG ol001, CNN KL10, each + prior-only)
sbatch --dependency=afterok:52032113:52032114 $W/analyze_idr_train.sbatch  # 52038652: analyze_idr.sh per chromosome + summary (reference geometric)
sbatch $W/summary_val_train_pt.sbatch                                       # 52057480: summary with phase-type as reference
```

The first summary used geometric as its reference, because geometric came first in the job's
model list. The phase-type-reference rerun is the one quoted. The train − val gaps do not depend
on the reference.

# Related Notes

- [IDR-labelled benchmark](fiber-idr-chip-benchmark.md): The val results and the memorization
  rationale this tests.
- [Prior-scale ablation](fiber-prior-scale-ablation.md): Val-only. The half-strength prior that
  ties phase-type there has not been run on train.
- [AG phase-type-structure run](fiber-ag-phase-type-structure-run.md): Val-only.

# Open Questions

- **Prior scale on train:** does the half-strength prior (AG ol001, s = 0.5) hold its val
  ties and wins on train? That would take one GPU job (about 40 min) plus analysis.
- **Pooling:** should the paper pool val + train for posterior metrics, given the gaps match the
  baselines'? A conservative choice is val as primary and train as a replication.
- **Prior-only memorization:** is the +0.02 prior-only train advantage on footprinting TFs real?
  Its per-TF bootstrap CIs (about ±0.02-0.05) overlap zero.

# Sources

- Outputs and scripts in the frontmatter.
