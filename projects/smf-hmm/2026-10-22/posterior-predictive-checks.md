---
title: Posterior predictive checks of the fixed-blocks SMF HMM
date: 2026-10-22
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - smf_hmm repo, smf_hmm/posterior_checks/ (uncommitted at time of writing)
  - smf_hmm repo, evaluation/posterior_checks/ (uncommitted; README.md documents every command)
  - Run /scratch/users/diamant/2026-09-29_fixed_blocks_global_kl3 (checkpoint, .hydra/config.yaml)
  - Outputs under /scratch/users/diamant/2026-09-29_fixed_blocks_global_kl3/posterior_checks/
  - TF/ChIP cache /scratch/users/diamant/mesc_validation/benchmark_cache/mesc_tf_chip_validation_v1
  - SMFNet_public, src/smf_net/evaluation/flank_baselines.py and data_preparation/fiberseq/README.md
  - User discussion, 2026-10-21 to 2026-10-22
tags: [smf-hmm, posterior-predictive-checks, misspecification, nucleosomes, tf, held-out, sequence-context, enzyme-bias, fiberseq]
---

# Summary

We built a posterior predictive check (PPC) toolkit for the SMF HMM and ran it
on `fixed_blocks_global_kl3` (fixed-block transitions, CNN prior with global
normalization, occupancy KL weight 3) using the mESC val chromosomes (chr8,
chr16). The main methodological finding: **conditional replicate checks at
observed sites are nearly circular for this model.** Posterior states are fit
to the same calls, and the emissions sit at 0.01 / 0.99, so replicates
x_rep ~ p(x | z) copy the observed call at 97.8% of sites. Checks that compare
latent statistics (prior vs pooled posterior) and **held-out site predictions**
(each site predicted from all others) are the informative ones.

Main model findings:

- **Every learned emission sits on its constraint bound.** open 0.010,
  strong/other_long/histone 0.990, weak TF 0.200 (protected probability).
- **Weak TF calls are almost data-free.** 66-70% of weak elements cover no
  observed site, and of those with one site, 81% of sites are accessible.
- **Nucleosome spacing is misfit even though lengths match.** The prior
  (nucleosome + geometric open gap) has no repeat length; the posterior has an
  excess at 160-220 bp (median 188 vs 194 bp, full val, 2114 bp windows).
- **Learned durations pile up against their bounds** (histone at the 130 bp
  minimum with a learned tail mean of ~8.5 bp vs the configured 21;
  TF_short_strong 50% at 8 bp; TF_long_strong 25% at 40 bp; other_long 29% at
  120 bp; 121-129 bp is reachable by no group).
- **Held-out calibration is overconfident at both extremes** and **nucleosome
  leakage differs by site type**: confidently-histone GpC sites are methylated
  2.55% vs CpG 1.25% (one histone emission predicts 2.11% / 1.65%).
- **Local sequence context adds little beyond the model on this SMF data**
  (best held-out deviance gain 0.002 per site; +-5 bp no better than +-3 bp),
  but dense C clusters (e.g. `CG[C]NN`) are methylated more than predicted in
  open DNA.

# Key Points

- **Three comparisons, not one.** A: prior predictive (z ~ p(z | DNA),
  x ~ p(x | z)); B: posterior with observed x (z ~ p(z | x_obs, DNA));
  C: posterior replicate (B's z, x_rep ~ p(x | z)). B vs C is the classic
  conditional PPC and is near-tautological here. A vs B (pooled posterior vs
  prior) is valid but attenuated, because z moves to absorb misfit; the
  absorbed part shows up in latent statistics (occupancy, lengths, spacing).
- **Held-out predictions remove the circularity.** Sites are split into K = 5
  folds of contiguous blocks; each fold is predicted from exact posterior
  marginals with that fold hidden. Tested for leakage, agreement with masked
  marginals, and calibration on model-simulated data. Block size is a knob:
  1 for SMF, larger for Fiber-seq (kinetic artefacts shared by neighbors).
- **Thresholding on states inferred with the site is itself circular.**
  Confident (>0.9) full-data calls give 1.000 / 0.001 accessible rates, more
  extreme than the emissions; with the site hidden, open is 97.3-98.5% and
  histone 1.2-2.6% accessible, and context differences appear.
- **Center-to-center spacing is the right nucleosome array statistic**:
  symmetric edge trimming changes called lengths but not centers.
- **Motif-anchored profiles** (ctcf, esrrb, klf4, oct4, ChIP quantile strata)
  show real footprints (CTCF top-10% ChIP: 0.61 -> 0.37 accessible over the
  motif), but observed vs replicate agreement there is uninformative (the copy
  effect). A held-out version is the planned fix.

# Details

## Code layout (smf_hmm repo, uncommitted)

- `smf_hmm/posterior_checks/` (library, tested in `tests/test_posterior_checks.py`, 26 tests):
  `replicates.py` (FFBS posterior paths + replicates; prior paths),
  `segments.py`, `edge_profile.py`, `element_statistics.py` (length, spacing),
  `site_counts.py`, `occupancy.py`, `anchor_profile.py`, `heldout.py`,
  `context_calibration.py`.
- `evaluation/posterior_checks/` (CLIs with `compute` / `plot` or `fit`;
  shared run loading in `run_data.py`; usage in `README.md`):
  `edge_profile.py`, `element_statistics.py`, `site_counts.py`,
  `occupancy.py`, `motif_profile.py`, `context_calibration.py`.
- Defaults: seeded random subset of 500 val regions (`--num-regions`, -1 for
  all), data config read from the run's `.hydra/config.yaml`.
- `statsmodels` added to `sherlock_env.yaml` / `local_env.yaml` (already in the
  installed `hmm` env). Clustered covariances are computed directly because
  statsmodels does not fully support `cov_type="cluster"` with `var_weights`.

## Findings by check

| Check | Result |
| --- | --- |
| Edge profile, histone (B vs C, full val) | Called edges sit next to the last accessible site (0.99 accessible at offset -1). Small residual: +0.3 pp accessibility 5-40 bp inside each edge, -0.4 pp in the core (edge ~1.3% vs core ~0.6%). |
| Lengths, histone (A vs B) | Posterior = prior (median 135, 12% at 130 bp). The learned prior shortened the tail (per-bp stay 0.883 vs configured 0.952). |
| Spacing, histone (A vs B, 2114 bp windows, full val) | Prior median 194 / mode 147; posterior median 188 / mode 148, excess at 160-220 bp, deficit at 250-400 bp. Suggests a linker macrostate. |
| Site counts (B vs C) | Weak TF: 66-70% of elements have n = 0 sites. Mild overdispersion at n = 2-3 for weak groups; strong groups match trivially. |
| Occupancy (A vs B) | Posterior/prior within 0.96-1.08; data want slightly more open and TF_long, less histone; prior already off the KL targets in the same direction. |
| TF lengths (A vs B) | Posterior = prior for every TF group (low power), but learned priors pile up at duration bounds. |
| TF edge profiles (B vs C) | Flat (strong ~1% accessible, weak ~81% vs 80%); boundaries absorb any center-vs-edge effect. |
| Copy rate (B vs C) | x_rep = x at 97.8% of observed sites (100 val regions, 476k sites). |
| Single-site leverage | Share of >0.9 calls still >0.9 with the site hidden: open 61%, histone 67%, TF_long_strong 57%, TF_short_strong 19%. |
| Held-out reliability | Well calibrated in the middle; predicted 99% accessible -> ~97% observed, predicted 1% -> 1.5-2%. |
| Context, confident held-out calls (p >= 0.9) | See table below. |

Minority-outcome rates for confident held-out calls (2000 val regions, 2.47M
site-reads):

| Site | State | Outcome | Observed | Expected | O/E |
| --- | --- | --- | --- | --- | --- |
| CpG | open | unmethylated | 2.11% | 2.30% | 0.92 |
| GpC | open | unmethylated | 2.74% | 3.24% | 0.84 |
| GCG | open | unmethylated | 1.50% | 2.07% | 0.73 |
| CpG | histone | methylated | 1.25% | 1.65% | 0.76 |
| GpC | histone | methylated | 2.55% | 2.11% | 1.20 |
| GCG | histone | methylated | 1.86% | 1.72% | 1.08 |

Context regression (offset logistic, flank encoding of SMFNet's
`flank_baselines.py`, region-clustered SEs): largest effects are a C at -2 for
GpC in open DNA (+0.35 log-odds; `CG[C]` misses about half as often) and C at
-2 / +2 for GCG in open DNA (+0.27 / +0.21); C at +2 lowers CpG/GCG leakage in
nucleosomes (-0.16). Units are log-odds of methylation relative to the
held-out prediction (exp gives an odds ratio). Cross-chromosome deviance
(fit chr8 -> score chr16 and reverse) improves by at most 0.002 per site.
Caveat: the SMF calls are folded onto the plus strand, so orientation-specific
effects are averaged; enzyme preference and ONT calling bias are not
separable without controls.

## Implications for Fiber-seq (PacBio, Hia5 m6A)

- With constant per-state emissions, enzyme bias is absorbed into states:
  poorly-methylated contexts in open DNA become short footprints, and the
  DNA-conditioned CNN prior can learn those k-mers as "binding". In this SMF
  run the prior already puts 1.8x more TF_short_strong occupancy on CpG than
  GpC sites (cause unknown).
- Plan: held-out checks as default; per-k-mer emission offsets in the model
  (shared across states, or split open / protected); within-molecule
  opposite-strand contrasts as a control (each HiFi molecule reports both
  strands; T sites are bottom-strand adenines); same-strand vs
  opposite-strand neighbor-call correlation to detect kinetic crosstalk;
  larger held-out blocks. No naked-DNA control is available yet.

## Commands

Environment: `PY=/home/groups/btrippe/diamant/miniforge/envs/hmm/bin/python`,
run from the smf_hmm repo root with `MPLCONFIGDIR=$PWD/.mplconfig`, on a GPU
allocation (sh04-03n07, H100; Slurm jobs 48546604 and 49031037).
`RUN=/scratch/users/diamant/2026-09-29_fixed_blocks_global_kl3`.

```bash
"$PY" -m evaluation.posterior_checks.edge_profile compute --run-dir "$RUN" --output-dir "$RUN/posterior_checks" --num-regions -1
"$PY" -m evaluation.posterior_checks.element_statistics compute --run-dir "$RUN" --output-dir "$RUN/posterior_checks/windows_2114" --statistics length spacing --output-len 2114 --num-regions -1
"$PY" -m evaluation.posterior_checks.site_counts compute --run-dir "$RUN" --output-dir "$RUN/posterior_checks/tf_500"
"$PY" -m evaluation.posterior_checks.occupancy compute --run-dir "$RUN" --output-dir "$RUN/posterior_checks/tf_500"
"$PY" -m evaluation.posterior_checks.element_statistics compute --run-dir "$RUN" --output-dir "$RUN/posterior_checks/tf_500"
"$PY" -m evaluation.posterior_checks.edge_profile compute --run-dir "$RUN" --output-dir "$RUN/posterior_checks/tf_edge_all" --groups TF_short_strong TF_short_weak TF_long_strong TF_long_weak --flank 20 --depth 40 --num-regions -1
"$PY" -m evaluation.posterior_checks.motif_profile compute --run-dir "$RUN" --output-dir "$RUN/posterior_checks/motif_profiles" --num-regions -1
"$PY" -m evaluation.posterior_checks.context_calibration compute --run-dir "$RUN" --output-dir "$RUN/posterior_checks/context" --num-regions 2000
JAX_PLATFORMS=cpu "$PY" -m evaluation.posterior_checks.context_calibration fit "$RUN/posterior_checks/context/context_sites.npz" --title "kl3 val, 2000 regions"
```

Deviations: the first edge-profile and smoke runs used the then-default
first-N regions (since replaced by a seeded random subset); the first
context `fit` scored +-3 bp fits against +-5 bp tokens (fixed, regression
test added). Ad hoc diagnostics (copy rate, leverage, A vs B by context) were
temporary scripts under `/tmp` and are not preserved; their numbers are
recorded above.

Outputs live on `$SCRATCH` (purged after 90 days without modification); copy
anything to keep to `$GROUP_HOME` or Oak.

# Related Notes

- [Posterior-mean TF/ChIP evaluation](../2026-09-23/posterior-mean-tf-chip-evaluation.md): the ChIP-based evaluation these checks complement; same mESC val setup.
- [Multi-TF motif caches](../2026-09-16/multi-tf-motif-cache-and-state-call-visualization.md): source of the motif sites and ChIP strata used by the motif-anchored profiles.
- [SMFNet flank baseline](../../smfnet/2026-10-05/flank-baseline-logistic-regression.md): the context parametrization reused in the context-calibration regression.

# Open Questions

- Does a linker macrostate in `FixedBlockTransition` reproduce the ~185 bp
  spacing, and does a 147 bp histone minimum (with or without edge-specific
  emissions) fit better than 130 bp? Needs training runs.
- Held-out motif-anchored profiles: does the model under-predict protection at
  motif centers of bound (high-ChIP) sites?
- Is the CNN prior's CpG enrichment of TF_short_strong biology or learned
  assay bias?
- Should emissions be per site type (GpC vs CpG nucleosome leakage differs
  ~2x), and should the emission bands be widened given every parameter is at
  its bound?
- Molecule-level methyltransferase efficiency (per-read overdispersion of
  confident-open misses) is not yet checked.

# Sources

- smf_hmm repo: `smf_hmm/posterior_checks/`, `evaluation/posterior_checks/README.md`, `tests/test_posterior_checks.py`
- Run: `/scratch/users/diamant/2026-09-29_fixed_blocks_global_kl3`
- Outputs: `/scratch/users/diamant/2026-09-29_fixed_blocks_global_kl3/posterior_checks/` (`histone_edge_profile.pdf`, `windows_2114/histone_{length,spacing}.pdf`, `tf_500/{site_counts,occupancy}.pdf`, `tf_edge_all/`, `motif_profiles/motif_profiles.pdf`, `context/context_{residuals,reliability,effects_flank3,effects_flank5}.pdf`, `context/context_deviance.json`)
- SMFNet_public: `src/smf_net/evaluation/flank_baselines.py`, `data_preparation/fiberseq/README.md`
