---
title: What causes the spurious TF-length footprints? Mostly the AG prior — it adds minimum-length (8 bp) TF nubs beside nucleosomes and on poorly labelled sequence, where reads are least informative; emissions are calibrated, other TFs explain little, and the short histone is a plot artifact
date: 2026-11-16
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /scratch/users/diamant/fiber_eval_work/tf_mechanism/ (gpu_diagnostics.py and .sbatch; causal_chain.py; segments.py; heldout_calibration.py; plot_paths.py; bm/tf_by_bm.py; edges/prior_vs_edges.py; make_figures.py; out/, paths/, figures/)
  - /scratch/users/diamant/fiber_eval_work/prior_tf_bumps/ (bias_vs_bumps.py, efficiency_tracking.py, homopolymer_test.py, bump_sources.py)
  - Slurm: 52505020/52505021 (structure, synthetic, examples), 52505759 (gelu_s0 reads + heldout), 52505761 (collapsed reads); fix trials 52509099/52509102/52509105 (+ inference 52509100/03/06, evaluation 52509108)
  - User requests and observations, 2026-11-16
tags: [smf-hmm, fiberseq, tf-calls, prior, nucleosome, labelling-efficiency, overview]
---

# Summary

**Question** (user): what causes the TF-length footprints away from benchmark motifs? They are
suspiciously often in low-efficiency stretches, and many look like histones split into other_long + TF.

**Models.**
- **gelu_s0:** GELU-head ablation B, seed 0 (spacer + TF/3, no TF chaining, other_long 40-100 bp, 147 bp
  tapered histone, additive context offsets, KL 10).
- **gelu_s0 at prior ×0:** the same emissions and transitions with no sequence prior.
- **collapsed:** the ReLU seed-2 run whose prior collapsed to a constant.

**Data.** 259k val reads with posterior path samples, plus 8.2M held-out calls.

| hypothesis | verdict | note |
|---|---|---|
| "Very short histones" violate the 147 bp minimum | **No: plot artifact.** Every sampled interior histone is 147 bp; argmax-of-marginals backgrounds make short histone patches where histone and other_long + TF compete | [short histone](fiber-short-histone-plot-artifact.md) |
| Insufficient m6A bias correction in low-efficiency stretches | **Largely refuted.** Held-out calls are calibrated within ±0.5 points across efficiency; without the prior, TF calls *fall* in poorly labelled DNA; the prior multiplies them 11× there. Long A/T runs are slightly under-corrected (minor) | [bias correction](fiber-tf-nubs-bias-correction.md) |
| The prior's TF bumps are other TFs | **Mostly refuted.** No enrichment at footprinting-motif matches at strict or relaxed thresholds; prior-added calls fall on motifs below chance; a few real exceptions (e.g. an E-box) | [other TFs](fiber-tf-nubs-other-tf-motifs.md) |
| Histones split into other_long + TF; prior puts TF at nucleosome edges | **Supported.** 72% of off-motif ×1 TF calls lie within 30 bp of a nucleosome (376 vs 54 per 1,000 reads without the prior); mostly 8 bp nubs 11-30 bp beyond a 147 bp histone; prior TF rises with nucleosome-edge density; noisy nucleosome-length stretches also fragment under the strict histone emission | [split nucleosomes](fiber-tf-nubs-split-nucleosomes-and-edges.md) |
| Deterministic spacer-TF-spacer transitions are favourable | **Refuted as stated** (an excursion costs about 8.5 nats). **Refined:** the prior's per-microstate potentials penalize TF position 9 (−1.3 nats), making exactly-8 bp TFs the preferred length (39% vs 9% under transitions alone) | [spacer and durations](fiber-tf-nubs-spacer-and-duration-shaping.md) |
| Which part of the prior: TF placement or the open/histone landscape? (diagnostic) | **TF placement.** Flattening it along each window removes 86% of off-motif TF (14.4 → 2.0 per read-kb) with unchanged AUROC/AUPRC; flattening the landscape instead raises AUROC by 0.016 | [prior flattening](fiber-tf-nubs-prior-flattening.md) |
| Under-methylated reads produce the TF calls | **Refuted for TF** (TF calls rise with labelling). Under-methylated reads over-use other_long (58 vs 25 per read-kb), and without a b_m term the model is miscalibrated by ±3-4 points | [under-methylated reads](fiber-tf-nubs-undermethylated-reads.md) |

**Overall picture.**
- **Where the excess comes from.** It is not in the emissions: GELU ×0 and the collapsed model call
  off-motif TF at roughly the model-free kernel's rate.
- **What the AG prior adds:**
  - about 2.3× more off-motif TF at every labelling level;
  - concentrated next to nucleosome ends and on A/T-rich, poorly labelled sequence;
  - shaped by its own per-microstate potentials into 8 bp nubs.
- **The nubs are weakly supported.**
  - Only 52-61% are supported by the read's own calls (kernel), against 93% at ×0.
  - The median TF-called base buys 0.19 nats; forbidding TF turns 69% of them into open.
  - Held-out calls inside confident ones are mostly protected (91.5%), so they are short protected
    remnants (edge or partial-nucleosome protection, inferred) rather than accessible DNA.

**Fix trials, round 1: none removes the excess.** GELU, design B, seed 0, about 3.5-4 h each; Slurm 52509099/102/105,
evaluation 52509108.
1. **tied:** TF and other_long chain states share one prior potential each (`transition.tied_blocks`).
2. **eff_filter:** labelling-efficiency emission (ε 0.003) and the 2% molecule-rate filter.
3. **tied_eff_filter:** both.

Results:
- **Tying barely changes the nubs.** The 8 bp share drops 32% → 26%, off-motif TF 18.2 → 17.2 per read-kb
  (scorecard), with the same AUROC (0.742) and slightly better prior log-odds (mean 0.79 vs 0.74). The
  duration cliff is a minor contributor; the nubs come mainly from TF placement (prior-flattening note).
- **The efficiency emission trims off-motif TF** (14.1 scorecard; 11.3 per read-kb in the specificity
  metric vs 14.4) and raises the on-motif share (0.30 vs 0.28). As with the ReLU head, it costs about
  0.03 AUROC (0.709) with the same matched AUPRC (0.492).
- **Both together:** off-motif 13.4 (scorecard), AUROC 0.713.

Round 2 (submitted): prior occupancy KL 300 and 1000 with the GELU head. It targets the excess TF mass
the prior keeps where reads cannot contradict it; the earlier KL 300-3000 runs used the collapsing
ReLU head. Slurm 52521012 / 52521015, evaluation 52521017.

| model (prior ×1) | AUROC | matched AUPRC | off-motif TF / read-kb | prior log-odds (mean) | 8 bp share of off-motif TF runs | kernel-supported runs |
|---|---|---|---|---|---|---|
| gelu_s0 | 0.741 | 0.488 | 14.4 | 0.74 | 32% | 59% |
| fix: tied | 0.742 | 0.491 | 13.4 | 0.79 | 26% | 64% |
| fix: efficiency + filter | 0.709 | 0.492 | 11.3 | 0.72 | 32% | 66% |
| fix: both | 0.713 | 0.487 | 10.6 | 0.78 | 28% | 72% |
| diagnostic: TF placement flattened | 0.740 | 0.483 | 2.0 | — | 22% | 90% |
| collapsed (constant prior) | 0.796 | 0.496 | 6.5 | 0.00 | 17% | 83% |

(Off-motif TF per read-kb is from the specificity metric; the 8 bp share and kernel support are from
nub_scorecard.py, every 4th window.)

# Key Points

- **What the sequence prior changes** (gelu_s0 ×1 vs ×0, confident TF segments per 1,000 reads):
  - adjacent to a nucleosome: 376 vs 54;
  - isolated without kernel support: 69 vs 2.2;
  - on benchmark motifs: 183 vs 63.
- **Off-motif TF calls by local efficiency:** ×1 rises from 1.05% to 3.53% of positions from the best-
  to the worst-labelled decile; ×0 falls from 0.54% to 0.32%.
- **Held-out calibration:**
  - by efficiency: within ±0.5 points;
  - by b_m: +3.0 to −3.9 points;
  - by state: off-motif TF 91.5% protected, histone 96.5%, open 38%.
- **Synthetic stretches:**
  - clean protected stretches get the "right" state at every length;
  - 10%-m6A stretches fragment, e.g. 180-210 bp → histone + gap + 25-33 bp TF.

# Details

- **The GPU diagnostics script** (`gpu_diagnostics.py`, one model per process) produces:
  - segment-duration distributions from the transition matrix;
  - Viterbi decoding of synthetic reads;
  - 50 posterior path samples for 44 example reads under ×1, ×0, collapsed and TF-forbidden;
  - one sample per read for 259k reads;
  - read log-likelihoods with and without TF;
  - split-half held-out predictions.
- **The collapsed ReLU checkpoint** loads only with the pre-GELU code, so it runs from a git worktree of
  f0401ff (`tf_mechanism/smf_hmm_head`) via `PYTHONPATH`.
- **Path-view figures** (`paths/`, one per example) show calls as ticks, the kernel and prior, and 8
  sampled paths per setting with segment lengths labelled.

```bash
T=/scratch/users/diamant/fiber_eval_work/tf_mechanism
MODELS=gelu_s0 sbatch $T/gpu_diagnostics.sbatch structure synthetic examples
MODELS=collapsed CODE_DIR=$T/smf_hmm_head sbatch $T/gpu_diagnostics.sbatch synthetic examples
MODELS=gelu_s0 sbatch $T/gpu_diagnostics.sbatch reads heldout
MODELS=collapsed CODE_DIR=$T/smf_hmm_head sbatch $T/gpu_diagnostics.sbatch reads
python $T/causal_chain.py; python $T/segments.py; python $T/heldout_calibration.py; python $T/bm/tf_by_bm.py
python $T/edges/prior_vs_edges.py; python $T/plot_paths.py; python $T/make_figures.py
# fix trials: /scratch/users/diamant/2026-11-16_fiber_gelu_B_{tied,eff_filter,tied_eff_filter}_alphagenome_bf16/train.sbatch
```

# Related Notes

- [Prior ablation and ReLU-head collapse](../2026-11-15/fiber-prior-ablation-and-relu-head-collapse.md):
  The GELU head and the collapsed model used here.
- [Prior TF log-odds and design takeaways](../2026-11-14/fiber-prior-tf-log-odds-and-design-takeaways.md):
  The prior metric used to judge the fix trials.
- [Nucleosome lengths and edges](../2026-11-13/fiber-nucleosome-lengths-and-edges.md): Wide nucleosome
  edges and shortening with labelling, background for the remnants.

# Open Questions

- KL 300 / 1000 fix-trial results (to be added).
- What are the protected remnants beside nucleosomes (chromatosome, partial unwrapping, real TFs at
  nucleosome borders)? Should they get their own state, or a longer/variable histone?

# Sources

- Scripts, jobs and outputs in the frontmatter.
