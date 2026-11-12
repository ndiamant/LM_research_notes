---
title: Scaling the AG prior at inference — weaker priors steadily recover phase-type's bound-vs-unbound AUROC, a half-strength prior ranks bound sites best, and the residual gap at prior 0 is structure or emissions
date: 2026-11-12
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - smf_hmm repo, branch nd/experiments (uncommitted) — smf_hmm/posterior_inference.py (InferenceConfig.prior_scale), smf_hmm/cli/infer_posterior.py (--prior-scale), tests/test_posterior_inference.py
  - Val posteriors /scratch/users/diamant/gm12878_validation/covariation_windows_v2/{chr20,chr8}-512-valid90-reads10-flank128/posteriors/{ag_kl10,ag_ol001}_s{0,025,05}.h5 (Slurm 52036336)
  - IDR analysis analysis_idr_prior_scale/ (Slurm 52038648); pooled summary /scratch/users/diamant/fiber_eval_work/idr_summary_prior_scale_val.tsv
  - Scripts /scratch/users/diamant/fiber_eval_work/{infer_prior_scale_val.sbatch, analyze_prior_scale_val.sbatch, analyze_idr.sh, idr_split_summary.py}
tags: [smf-hmm, fiberseq, gm12878, alphagenome, prior, ablation, idr, benchmark]
---

# Summary

`--prior-scale s` multiplies the sequence prior's per-position terms at inference:
- s = 1 is the trained posterior.
- s = 0 drops the DNA prior and keeps the transitions, the learned emissions and the ±3 bp context
  offsets.
- The initial distribution is unchanged; it only affects the first few bases.

On val, with IDR labels ([benchmark note](fiber-idr-chip-benchmark.md)):

- **Bound vs unbound improves monotonically as the prior weakens**, but at s = 0 it is still
  below phase-type on most TFs. On protected-read AUROC, AG ol001 at s = 0 is close on NRF1 (0.91
  vs 0.92), CREM (0.87 vs 0.88) and RFX5 (0.79 vs 0.79), but behind on USF2 (0.80 vs 0.87), SPI1
  (0.76 vs 0.83) and YY1 (0.67 vs 0.76). So part of the gap is the prior, and the rest comes from
  the states and emissions (other_long, fixed 147 bp histone, learned emissions).
- **Ranking bound sites by ChIP is best at an intermediate scale.** AG ol001 at s = 0.5 is
  significantly above phase-type on CTCF (0.28 vs 0.24) and PAX5 (0.08 vs −0.07). AG KL10 at
  s = 1 is below phase-type on CREM, EBF1, SPI1 and YY1; s = 0.25 removes all four deficits
  and s = 0.5 all but YY1.
- **With accessibility controlled, AG ol001 at s = 0.5 ties phase-type on bound vs unbound** for
  CTCF 0.81, NRF1 0.75, USF2 0.64, CREM 0.65, PAX5 0.64 and EBF1 0.65. It remains below on SPI1
  (0.67 vs 0.71) and YY1 (0.44 vs 0.58).

So the trained prior is too strong relative to the reads. A half-strength prior keeps its ranking
benefit and loses most of its cost. This supports the "limit what the prior controls / prior
dropout" ideas, and gives a cheap inference-time setting to report.

# Key Points

AUROC within accessibility quintiles and protected-read AUROC (bound vs unbound), and partial
Spearman within bound sites. "+" / "−" mark a significant difference from phase-type (95% paired
bootstrap).

| TF | metric | phase-type | AG KL10 s0 | s0.25 | s0.5 | s1 | AG ol001 s0 | s0.5 | s1 |
|---|---|---|---|---|---|---|---|---|---|
| CTCF | AUROC within accessibility quintiles | 0.80 | 0.78− | 0.79 | 0.79 | 0.77− | 0.79 | 0.81 | 0.78− |
| CTCF | protected-read AUROC | 0.89 | 0.87− | 0.87− | 0.87− | 0.86− | 0.88− | 0.88− | 0.87− |
| CTCF | within-bound ρ | 0.24 | 0.26 | 0.27 | 0.27 | 0.25 | 0.27+ | 0.28+ | 0.25 |
| NRF1 | protected-read AUROC | 0.92 | 0.89− | 0.89− | 0.88− | 0.87− | 0.91− | 0.90− | 0.89− |
| NRF1 | within-bound ρ | 0.45 | 0.44 | 0.46 | 0.47 | 0.45 | 0.48 | 0.48 | 0.44 |
| USF2 | protected-read AUROC | 0.87 | 0.80− | 0.81− | 0.81− | 0.78− | 0.80− | 0.81− | 0.80− |
| CREM | within-bound ρ | 0.38 | 0.35 | 0.38 | 0.33 | 0.25− | 0.37 | 0.37 | 0.28 |
| PAX5 | within-bound ρ | −0.07 | 0.03+ | 0.07+ | 0.05+ | −0.07 | 0.02+ | 0.08+ | −0.07 |
| SPI1 | AUROC within accessibility quintiles | 0.71 | 0.65− | 0.66− | 0.65− | 0.61− | 0.65− | 0.67− | 0.62− |
| YY1 | protected-read AUROC | 0.76 | 0.66− | 0.63− | 0.61− | 0.59− | 0.67− | 0.61− | 0.60− |

The full tables (all TFs, raw AUROC, CIs) are in `idr_summary_prior_scale_val.tsv`. The CNN
KL10 prior alone was also scored: it is below phase-type on bound vs unbound everywhere, but its
within-bound ρ is 0.51 on NRF1 (phase-type 0.45, not significant).

The [read m6A-rate note](fiber-read-m6a-rate-evaluation.md) shows the s = 0.5 model still
shifts `other_long` with molecule labeling rate, just less.

## Command provenance

```bash
W=/scratch/users/diamant/fiber_eval_work
sbatch $W/infer_prior_scale_val.sbatch     # 52036336, akundaje GPU: --prior-scale {0, 0.25, 0.5} (AG KL10), {0, 0.5} (AG ol001)
sbatch --dependency=afterok:52036336 $W/analyze_prior_scale_val.sbatch   # 52038648: analyze_idr.sh + idr_split_summary.py
```

The `--prior-scale` option and its two tests (`tests/test_posterior_inference.py`, 10 passing)
are uncommitted on `nd/experiments`.

# Related Notes

- [IDR-labelled benchmark](fiber-idr-chip-benchmark.md): The metrics and the s = 1 results.
- [Read m6A-rate evaluation](fiber-read-m6a-rate-evaluation.md): Uses the s = 0 and s = 0.5
  posteriors.
- [Fiber neural HMMs vs phase-type](../2026-11-11/fiber-neural-vs-phase-type-nucleosome-free.md):
  The earlier finding that AG's posterior follows its prior.

# Open Questions

- **Phase-type-structure run:** does the AG model trained with phase-type states and a TF cap
  (training run 2026-11-12_fiber_single_tf_phasetype_tf7to40_hist130plus_kl10_context3_alphagenome_bf16,
  Slurm 52036412, running at writing) close the s = 0 gap? And what is its best prior scale?
- **Training-time tempering:** would prior dropout or tempering during training give the same
  benefit as tempering at inference?
- **Train chromosomes:** do the same scales hold there (train IDR analysis pending)?

# Sources

- Outputs and code in the frontmatter.
