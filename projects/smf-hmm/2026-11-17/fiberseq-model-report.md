---
title: Building the GM12878 Fiber-seq HMM — components, the evidence behind each choice, and metric progression (living report draft, source for a LaTeX/PDF report)
date: 2026-11-17
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - The smf-hmm research notes from 2026-10-22 to 2026-11-17 (index README.md)
  - /scratch/users/diamant/fiber_eval_work/report_metrics/*.tsv (metric tables collected from idr_summary_*_val.tsv, kernel_ref/specificity*/, prior_log_odds/val_*priors.tsv)
  - Run configs /scratch/users/diamant/2026-10-22_fiber_fixed_blocks_global_kl3 to 2026-11-17_fiber_hist190_tf60_kl300_alphagenome_bf16
tags: [smf-hmm, fiberseq, report, model-design, metrics]
---

# Building the GM12878 Fiber-seq HMM: components, evidence, and metrics

*Living report, the source for a LaTeX/PDF version. Started 2026-11-17; figures go in the companion folder `fiberseq-model-report/`. It deviates from the standard note body structure because it is a report.*

## Status and conventions

**Status.** This is a draft. Numbers are copied from the research notes and evaluation tables listed in
[Sources](#sources). The latest model, TF 15-60 (`hist190_tf60`), is included. Comments in
*[brackets]* mark editorial to-dos for the LaTeX version.

**Conventions.**
- **Chromosomes.** Every metric is on the validation chromosomes (chr8, chr20) unless marked "train"
  (chr2, chr11, chr12, chr17, chr19).
- **Test set.** The test chromosomes (chr1, chr3, chr6) were held out for the whole project and never
  evaluated.
- **Prior scale.** All posterior metrics use the trained prior as is (prior ×1); we make no post hoc
  prior adjustments.

---

## 1. Summary

*[One-paragraph abstract for the LaTeX version.]*

**The current model** is `hist100_190`; run directory `2026-11-17_fiber_hist100_190_kl300_alphagenome_bf16`.
It is a fixed-topology HMM with four state groups: open, TF, other_long and histone.
- **TF block:** 15-40 bp, with 10 bp open spacers and no chaining.
- **other_long block:** 40-100 bp, occupancy 1%.
- **Nucleosome:** variable length, 100-190 bp. Lengths are peaked at 147 bp, and the 40 bp tapered
  edges sit at both ends whatever the length.
- **Emissions:** banded per-group protection probabilities, plus learned ±3 bp m6A context offsets.
- **Prior:** sequence-conditioned. A frozen AlphaGenome encoder feeds a residual GELU head, which
  outputs one tied potential per state group, under global normalization.
- **Training objective:** read likelihood plus a weight-300 KL from the prior's occupancy to the
  target occupancy.

**Evolution.** It grew out of the first Fiber-seq model (2026-10-22), which reused the mESC SMF design
unchanged:
- six blocks, four of them TF groups;
- a CNN prior with KL weight 3;
- no context correction.

**Headline progress.** The table compares three stages:
- the first neural model on the IDR benchmark, AG KL 10, before the 2026-11-11 other_long fix;
- the best model of the old design, `ag_sptf`, 2026-11-12;
- the current model.

All are at prior ×1, with the phase-type baseline for reference.

| metric | phase-type (baseline) | AG KL 10 | `ag_sptf` | current `hist100_190` |
|---|---|---|---|---|
| IDR AUROC | **0.809** | 0.700 | 0.772 | 0.765 |
| AUROC, accessible reads | 0.685 | 0.629 | 0.663 | **0.687** |
| AUROC, protected reads | **0.826** | 0.744 | 0.790 | 0.781 |
| ρ within bound sites | 0.165 | 0.111 | 0.169 | **0.215** |
| accessibility-matched AUPRC | 0.485 | – | – | **0.535** |
| off-motif TF per read-kb | 38.6 | – | – | **4.8** |
| prior TF log-odds (9-TF mean) | – | 0.38 | 0.80 | **1.09** |

**Caveat to state prominently.** On plain bound-vs-unbound AUROC, nothing has beaten phase-type
(0.809), and the collapsed constant-prior model also scores higher (0.796). Those AUROCs reward
over-calling TF: phase-type calls 38.6 off-motif TF per read-kb. With bias-corrected emissions it
falls to 0.734. Our progress shows up as:
- accessible-read AUROC and within-bound ranking;
- precision at matched recall (AUPRC);
- far fewer spurious calls (off-motif TF);
- a sequence prior that separates bound from unbound motifs (prior log-odds).

The case for these metrics is in §2.4.

---

## 2. Data, splits, and evaluation

### 2.1 Data
- **Reads.** GM12878 Fiber-seq (PacBio, m6A), from
  `/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/active.h5`; reference hg38.
- **Assay `fiberseq_m6a`.** Every A and T is a callable site. Training applies reverse-complement
  augmentation with probability 0.5.
- **Training windows.** 512 bp crops from 2,114 bp regions; 32 reads per region and 128 regions per
  batch. A read must have at least 20 calls in the window.
- **Context bias.** m6A calling has a strong sequence-context bias.
  - Per-site m6A rates run from 0.002 (`AA[A]AA`) to 0.24 (`CC[A]GC`), so the 95th-percentile context
    has 64× the odds of the 5th.
  - ±3 bp captures 95% of the best held-out gain.
  - One strand-symmetric table suffices (A-site vs T-site correlation 0.9998).
  - Source: 2026-10-30 context-bias note.
- **Blind positions.** In `AAA[A]AAA` m6A is 0.3% even in the most open reads (2026-11-10). Long A/T
  runs therefore look protected whatever their occupancy.

### 2.2 Splits
- **Train:** chr2, chr4, chr5, chr7, chr9-19, chr21, chr22, chrX.
- **Validation:** chr8, chr20; used for checkpoint selection and every reported metric.
- **Test:** chr1, chr3, chr6; held out and never evaluated.
- **Replication:** the train-chromosome IDR benchmark (chr2, chr11, chr12, chr17, chr19) found no
  posterior memorization. The train − val raw AUROC gap is +0.009 to +0.012 for every model, the
  phase-type baseline included. Prior-only scores show a small train advantage of about +0.02
  (2026-11-12).

### 2.3 TF benchmark
- **TF selection.** A genome-wide masked footprint screen of 91 motifs (2026-11-10)
  split the 9 benchmark TFs into two groups:
  - **Footprinting TFs:** CTCF and NRF1 strongly; USF2 (E-box) and CREM (CRE) weakly; RFX5 weakly and
    exploratory, with only 23 bound val sites.
  - **ChIP'd but no detectable footprint:** EBF1, PAX5, SPI1, YY1.
- **Site labels (2026-11-12).** From ENCODE conservative IDR peaks:
  - **bound:** a summit within 100 bp of the motif;
  - **unbound:** no summit within 1 kb;
  - everything in between is excluded.
- **Windows.** 512 bp windows with at least 90% valid positions and at least 10 reads.

### 2.4 Metrics, and why each was added

| metric | definition | added | reason |
|---|---|---|---|
| Spearman ρ with pooled ChIP | read-averaged motif P(TF) vs ChIP | carried over from SMF | first benchmark |
| accessibility-controlled ρ, references | partial ρ given accessibility (±25 bp); accessibility / motif-score reference rows | 2026-11-04 (SMF) | the benchmark mostly measured prior plus accessibility |
| protected-read metric | mean P(TF) over reads with a protected call in the motif | 2026-11-05 (SMF) | separate footprint evidence from accessibility |
| P(TF+open) | nucleosome-free score | 2026-11-11 | neural models called bound but unprotected DNA open |
| IDR AUROC (raw, accessible-read, protected-read), ρ within bound, paired bootstrap CIs | bound vs unbound AUROC; within-bound ranking | 2026-11-12 | ChIP-quantile labels sat in the noise mode; IDR labels are standard and untuned |
| train-chromosome replication | same benchmark on 5 train chromosomes | 2026-11-12 | more power; memorization check |
| TF run / chain statistics | runs ≥ 35 bp; chains > 40, 60 or 100 bp | 2026-11-12 | detect TF chaining |
| accessibility-matched AUPRC | strong bound sites (≥ median ChIP) vs accessibility-matched unbound sites, 3:1, 100 resamples | 2026-11-13 | the AUROCs reward over-calling |
| off-motif TF per read-kb, on-motif share | TF-called bases more than 5 bp from benchmark motifs | 2026-11-13 | count spurious calls directly |
| debiased kernel P(bound) | model-free, ±7 bp, ε = 0.02, fit on train chromosomes | 2026-11-13 | model-free footprint reference (CTCF bound / unbound confident-footprint ratio 37×) |
| prior TF log-odds, prior-only AUROC | logit(median prior P(TF) at bound motifs) − logit(median at unbound) | 2026-11-14 | the read-dominated posterior hid large prior changes, e.g. a collapsed prior still scored 0.709 AUROC; standing metric for every new model |
| `prior_potential_sd/<group>`, seed SD | spread of prior potentials during training; spread across seed replicates | 2026-11-15 | detect a dead (constant) prior head |
| Viterbi length histograms, nub scorecard, open overdispersion | segment lengths per state; off-motif TF runs, minimum-length share, kernel support; clustered protection in open DNA | 2026-11-16 | diagnose spurious short TF calls and duration pile-ups |

*[LaTeX: state that AUROC differences of ≥ 0.005 are likely real; the GELU seed SD is about 0.001
AUROC and 0.02 log-odds (2026-11-15/16).]*

---

## 3. Starting point

### 3.1 Inherited from the mESC SMF work
- **Fixed-block HMM** (`FixedBlockTransition`, commit 7a3e2b3, 2026-09-28):
  - an open hub;
  - one chain of microstates per protected element, with frozen capped durations;
  - destination biases from a sequence prior.
- **Ordered-group emissions with absolute protection bands** (`OrderedGroupSMFEmissions`, bec22ba).
- **Global prior normalization** (b2b0d78). With zero potentials, the model equals the unconditioned
  HMM.
- **Prior occupancy KL objective** (bd301aa).
- **Posterior predictive checks** (2026-10-22). These found every emission on its band edge and the
  histone pinned at its minimum length, and predicted that Fiber-seq context bias would need explicit
  modelling.
- **147 bp fixed histone with an invaded-edge taper** (2026-11-03, d0d5c32). The design reasoning was
  to keep the histone fixed-length with deep edge invasion, rather than make it variable-length.
- **Learned ±k context logit offsets** (2026-11-04, d160bfc).
- **The occupancy-KL lesson** (2026-11-04). A TF target about 3.5× the rate the reads support created
  a weak-TF floor; this was the actual cause of a ChIP gain that had been put down to context bias.
- **Frozen AlphaGenome encoder prior** (2026-11-06, 4c239f7 / babfb52):
  - per-scale LayerNorm adapters;
  - gradient clipping 1.0, weight EMA 0.999, bf16 encoder;
  - validation NLL 12.208 vs 12.360 for a matched CNN.

### 3.2 First Fiber-seq model (2026-10-22, `fiber_fixed_blocks_global_kl3`)
- **Transition** `fixed_blocks`, six blocks:
  - TF_short strong and weak, 8-25 bp;
  - TF_long strong and weak, 20-40 bp;
  - other_long, 50-120 bp, occupancy 0.08;
  - histone ≥ 130 bp, mean 150, occupancy 0.55.
- **Emissions:**
  - open 0.3 in band [0.05, 0.55];
  - closed band [0.55, 0.995];
  - weak TF groups 0.75, strong 0.95, other_long 0.9, histone 0.97;
  - no context offsets, no taper.
- **Prior:** one-hot CNN `[128,128,128]` with ReLU, kernel 9, global normalization.
- **Objective:** KL weight 3. Targets: open 0.31, each TF group 0.015, other_long 0.08, histone 0.55.
- **Optimizer:** Adam lr 1e-3, no clipping, no EMA.

### 3.3 What was wrong with it (2026-10-30)
- **The DNA-only prior learned the m6A context bias.** Prior open rises from 0.000 to 0.262 across
  labelling-efficiency deciles.
- **States were repurposed to model labelling.** TF_short_weak sat exactly on its band floor (0.550),
  and open learned 0.091.
- **Benchmark (2026-11-10, CTCF ρ chr8 / chr20):** 0.53 / 0.52, vs 0.58 / 0.53 for phase-type.
  Accessibility alone beats every model on EBF1, PAX5, SPI1 and YY1.

---

## 4. Changes from the first model to the current one

### 4.1 Adopted changes (chronological)

| # | date | change | main evidence | source note |
|---|---|---|---|---|
| 1 | 11-09 | learned ±3 bp context offsets (`context_offset_k=3`) | the first model absorbed the context bias into states and prior (§3.3); offsets cost 16,384 parameters; SMF validation NLL 12.337 vs 12.356 | 10-30 bias diagnostic; 11-04 |
| 2 | 11-09 | one TF group instead of four (`fixed_blocks_single_tf`); TF protection floor 0.85 | the weak TF groups were labelling sinks; the SMF multi-TF ChIP gain came from the KL target, not the weak groups | 10-30; 11-04 |
| 3 | 11-10 | fixed 147 bp histone, 25 bp invaded-edge taper; other_long 40-100 bp | SMF: with the taper, the histone edge learned freely inside its band (0.27; was at its floor); the variable histone had pinned near its 130 bp minimum (PPC) | 11-03; 10-22 |
| 4 | 11-10 | bands: open [0.05, 0.45], closed [0.45, 0.995]; TF floor as `min_protected_probability` 0.85 | emissions sat on band edges (PPC); keep TF distinct from open | 10-22; a51fbef |
| 5 | 11-10 | AlphaGenome frozen-encoder prior (3 blocks, adapters 64, bf16), clip 1.0, EMA 0.999 | SMF: better likelihood and sequence prior; clipping prevents potentials of about 3×10⁵ | 11-06 |
| 6 | 11-11 | other_long occupancy 0.10 → 0.01, histone 0.60 → 0.69 (`ol001`) | other_long absorbed 16% of positions and many CTCF footprints (footprint-read P(TF) 0.49 vs 0.72 for phase-type); after: P(TF) 0.54, IDR metrics +0.00 to 0.03 | 11-11 other_long |
| 7 | 11-12 | IDR benchmark; train-chromosome replication | see §2.3-2.4 | 11-12 IDR; train vs val |
| 8 | 11-12 | 10 bp open spacers around TF (`spacer_length=10`) and TF target 0.02 → 0.0067 | largest single gain in the lineage: AUROC 0.716 → 0.772, ρ within bound 0.106 → 0.169, mean prior log-odds 0.44 → 0.80; confirmed on train; prior log-odds CTCF 0.87 → 1.94 from spacers; TF bases in chains > 40 bp 8% vs 26% for phase-type | 11-12 spacer |
| 9 | 11-13/14 | no TF chaining (`spacer_repeat_probability=0`) | without other_long, TF chains replaced it; forbidding chaining recovered most of the loss (`ag_nool` AUROC 0.633 → `ag_nochain` 0.715; off-motif 28.8 → 13.3) | 11-13 no-other_long; 11-14 |
| 10 | 11-14 | prior TF log-odds as a standing metric; no post hoc prior scaling | the posterior benchmark hid a collapsed prior | 11-14 |
| 11 | 11-15 | residual GELU prior head, h + conv(gelu(RMSNorm(h))), zero-init | the ReLU head died: 77-100% dead first-layer units; one seed collapsed to a constant prior (log-odds 0.00); GELU: 3 seeds, log-odds 0.74 / 0.71 / 0.71, SD 0.02 | 11-15 |
| 12 | 11-16 | prior occupancy KL weight 10 → 300 | KL 10 left the prior's TF occupancy free; at KL 300 AUPRC 0.488 → 0.525, off-motif 14.4 → 6.2, log-odds 0.74 → 0.83, AUROC 0.741 → 0.744 | 11-16 nub mechanism and spacer/duration notes |
| 13 | 11-17 | tied TF and other_long potentials (`tied_blocks`) | the prior's per-microstate potentials made 8 bp TF the preferred length (−1.3 nat penalty at state 9); tying cut TF calls at the minimum from 40% to 19%; log-odds 0.83 → 0.94 | 11-16 spacer/duration; ablation round 1 |
| 14 | 11-17 | TF minimum 8 → 15 bp | 8 bp TF "nubs" patched clustered protection in open DNA; short calls at weak-footprint TFs look like background; AUROC 0.744 → 0.756, ρ 0.227, kernel support 82% → 92% | 11-16 overdispersion; sparsity vs dissociation; round 1 |
| 15 | 11-17 | histone taper 25 → 40 bp with edge floor 0.2 | nucleosome edges are about 35 bp wide (11-13); round 1: AUROC +0.007, off-motif 6.18 → 5.84, other_long at the 100 bp cap 58% → 49% | 11-13 nucleosome edges; round 1 |
| 16 | 11-17 | combined 13 + 14 + 15 (`r2_full`) | AUROC 0.764, AUPRC 0.528, off-motif 5.65, log-odds 1.01; each part contributes (leave-one-out) | ablation round 2 |
| 17 | 11-17 | variable-length nucleosome, 100-190 bp: P(d) ∝ exp(−\|d − 147\|/15), 40-state exit tail, tied histone potentials | other_long held short nucleosomes (100 bp cap; spans 104-113 bp) and over-long ones (other_long at its minimum beside a histone; spans about 160-190 bp); result: AUROC 0.765, AUPRC 0.535, off-motif 4.78, log-odds 1.09, other_long 323 → 130 per 1,000 reads | short other_long; variable nucleosome |
| 18 | 11-17 | TF maximum 40 → 60 bp (`hist190_tf60`) | aimed at the isolated 41-60 bp other_long patches (27% on bound CTCF); result: AUROC 0.768, ρ within bound 0.232, AUPRC 0.538, off-motif 4.69, log-odds 1.08; but only 0.7% of TF calls exceed 40 bp (TF's zero-bias duration prior makes 41-59 bp 4-5× less likely than other_long's), so the gains likely come from removing the 40 bp cap's censoring | TF 15-60 note |

*[LaTeX: one short subsection per row, with a figure where one exists (§8).]*

### 4.2 Tried and not adopted

| change | result | source |
|---|---|---|
| cross-read kernel denoising before the HMM | own-read smoothing beats cross-read borrowing (AUROC 0.890 vs 0.878); the kernel biases the m6A rate | 10-30 kernel |
| sequence-dependent TF emission head | learned a CGC artifact (74% of its variance), not enzyme bias or TF identity | 11-02 (SMF) |
| phase-type state structure with the AG prior (`ag_pt`: TF 7-40, no other_long) | the prior chained capped TFs through short open gaps: 22% of TF bases in chains > 100 bp; NRF1 AUROC 0.74 vs 0.89 | 11-12 ag_pt |
| per-read m6A-rate filter 2-8% | drops about 75% of reads; hurts within-bound ranking (USF2 ρ 0.29 → −0.03) | 11-12 rate |
| molecule-rate filter (background rate < 2%) | drops 3.4% of reads; small effect in the ablation (0.769 vs 0.773); not in the current model | 11-13 cuts; 11-15 |
| labelling-efficiency emissions (per-read b_m) | better within-bound ρ (0.195 vs 0.169) but worse AUROC (0.730 vs 0.772; −0.03 in the ablation); degraded the prior (RFX5 log-odds 0.43 → 0.11); per-read b_m miscalibrated by ±3-4 points | 11-12 efficiency; 11-14; 11-15; 11-16 undermethylated reads |
| removing other_long | TF chains replace it (`ag_nool` AUROC 0.633, off-motif 28.8) | 11-13 |
| KL 3000 | prior collapsed: log-odds about 0 for every TF while the posterior still scored 0.709 | 11-13; 11-14; 11-15 |
| constant (sequence-free) prior | diagnostic: AUROC 0.790-0.796 but no prior information; explains why an informative prior can lower AUROC | 11-13; 11-15 |
| other_long narrowed to 45-80 or 50-90 bp | re-piles at its new ends; AUPRC −0.008 / −0.009 | rounds 1-2 |
| open band raised to 0.65 | open learns 0.488; AUROC, AUPRC and log-odds drop | round 1 |
| histone 100-147 only | AUROC 0.743 (−0.021), every TF lower; over-long nucleosomes remain | variable nucleosome |

---

## 5. Current model specification

*[LaTeX: one table per group, with a diagram of the transition topology.]*

**Data.** `fiberseq_m6a`; splits as §2.2; at least 20 calls per read; no molecule filter; no labelling
efficiency.

**Transition** (`FixedBlockTransition`, frozen edges, 351 states):

| block | duration | occupancy | zero-bias duration model |
|---|---|---|---|
| open | – | 0.2933 | hub state, plus spacer states |
| TF | 15-40 bp (15-60 in `hist190_tf60`) | 0.0067 | truncated geometric, 5% at the cap; mean 22.5 bp |
| other_long | 40-100 bp | 0.01 | truncated geometric; mean 58.6 bp |
| histone | 100-190 bp | 0.69 | P(d) ∝ exp(−\|d − 147\|/15); mean 146.7 bp; last 40 states form an exit tail |

- **Spacers.** TF is entered and left only through 10-state open spacers, with no re-entry
  (`spacer_repeat_probability=0`). Two TF calls are therefore at least 21 bp apart.
- **Tied potentials.** Each block's states share one potential (`tied_blocks=[TF,other_long,histone]`),
  and spacers use the open potential. The prior therefore outputs four distinct potentials per
  position.

**Emissions** (`OrderedGroupSMFEmissions`):
- **Bands.** Open in [0.05, 0.45]. TF in [0.85, 0.995]. other_long and histone in [0.45, 0.995].
- **Histone edges.** A 40 bp log-odds taper at both ends; the learned edge value lies between 0.2 and
  the core value.
- **Context offsets.** A ±3 bp table (4⁷ entries), centred and added to every state's logit.

**Prior** (`AlphaGenomePriorArchitecture`):
- **Encoder.** Frozen AlphaGenome embedder plus 3 blocks, run at 1 bp in bf16. It gives 4 scales
  (768-1152 channels).
- **Adapters.** LayerNorm, then linear to 64 per scale, concatenated to 256.
- **Head.** A linear layer to 128, then two residual GELU blocks (kernel 9, zero-init), then a final
  RMSNorm + GELU, then a zero-init linear output.
- **Normalization.** Global.

**Objective.** Mean read NLL + 300 × KL(target ‖ prior occupancy). Targets for open / TF / other_long /
histone are 0.2933 / 0.0067 / 0.01 / 0.69, equal to the zero-bias stationary occupancies.

**Optimizer.** Adam lr 1e-3, global-norm clipping 1.0, weight EMA 0.999 (used for validation and
checkpointing), 100 epochs, seed 0.

---

## 6. Metrics progression

### 6.1 Main lineage, prior ×1

*AUPRC and off-motif TF were not computed for models before 2026-11-13.*

| model (date) | key change | AUROC | AUROC acc. | AUROC prot. | ρ bound | AUPRC | off-motif TF / kb | prior log-odds |
|---|---|---|---|---|---|---|---|---|
| phase-type (baseline) | fixed phase-type HMM, no prior | 0.809 | 0.685 | 0.826 | 0.165 | 0.485 | 38.58 | – |
| phase-type, efficiency-corrected | bias-corrected emissions | 0.734 | 0.636 | 0.758 | 0.162 | 0.469 | 17.29 | – |
| CNN KL 10 (11-10) | single TF, CNN prior | 0.729 | 0.647 | 0.770 | 0.133 | | | 0.35 |
| AG KL 10 (11-10) | AlphaGenome prior, 147 bp tapered histone | 0.700 | 0.629 | 0.744 | 0.111 | | | 0.38 |
| AG ol001 (11-11) | other_long 1% | 0.716 | 0.639 | 0.756 | 0.106 | | | 0.44 |
| AG spacer (11-12) | 10 bp spacers | 0.769 | 0.653 | 0.789 | 0.138 | | | 0.76 |
| AG TF/3 (11-12) | TF target 0.0067 | 0.707 | 0.639 | 0.747 | 0.095 | | | 0.41 |
| `ag_sptf` (11-12) | spacers + TF/3 | 0.772 | 0.663 | 0.790 | 0.169 | | | 0.80 |
| `ag_nochain` (11-13) | no chaining, no other_long, efficiency, KL 10 | 0.715 | 0.639 | 0.740 | 0.163 | 0.475 | 13.25 | 0.65 |
| base B, ReLU (11-14) | `ag_sptf` + no chaining + other_long 40-100 | 0.751 | 0.655 | 0.771 | 0.137 | 0.494 | 12.37 | 0.74 |
| base B, collapsed seed (11-14) | dead head = constant prior | 0.796 | 0.663 | 0.808 | 0.180 | 0.496 | 6.50 | 0.00 |
| GELU B (11-15) | residual GELU head | 0.741 | 0.649 | 0.765 | 0.126 | 0.488 | 14.36 | 0.74 |
| GELU B, KL 300 (11-16) | KL 10 → 300 | 0.744 | 0.666 | 0.766 | 0.208 | 0.525 | 6.18 | 0.83 |
| + tie (11-16) | tied potentials | 0.751 | 0.670 | 0.770 | 0.202 | 0.525 | 6.19 | 0.94 |
| + TF 15-40 (11-16) | TF minimum 15 | 0.756 | 0.674 | 0.773 | 0.227 | 0.527 | 6.06 | 0.90 |
| + flank (11-16) | taper 40, floor 0.2 | 0.751 | 0.668 | 0.770 | 0.198 | 0.526 | 5.84 | 0.85 |
| `r2_full` (11-17) | tie + TF 15-40 + flank | 0.764 | 0.681 | 0.780 | 0.213 | 0.528 | 5.65 | 1.01 |
| histone 100-147 (11-17) | short nucleosomes only | 0.743 | 0.674 | 0.761 | 0.222 | 0.532 | 5.18 | 1.10 |
| `hist100_190` (11-17) | variable nucleosome 100-190 | 0.765 | 0.687 | 0.781 | 0.215 | 0.535 | 4.78 | 1.09 |
| **`hist190_tf60` (11-17)** | **TF 15-60** | **0.768** | **0.688** | **0.783** | **0.232** | **0.538** | **4.69** | **1.08** |

### 6.2 Per-TF AUROC, prior ×1

| TF | phase-type | AG KL 10 | AG ol001 | `ag_sptf` | GELU B | GELU B KL 300 | `r2_full` | `hist100_190` |
|---|---|---|---|---|---|---|---|---|
| CTCF | 0.884 | 0.840 | 0.850 | 0.876 | 0.862 | 0.866 | 0.868 | 0.861 |
| NRF1 | 0.894 | 0.786 | 0.818 | 0.886 | 0.832 | 0.840 | 0.868 | 0.859 |
| USF2 | 0.847 | 0.718 | 0.748 | 0.796 | 0.774 | 0.793 | 0.788 | 0.783 |
| CREM | 0.844 | 0.736 | 0.744 | 0.803 | 0.770 | 0.740 | 0.825 | 0.808 |
| RFX5 | 0.792 | 0.710 | 0.706 | 0.763 | 0.748 | 0.738 | 0.753 | 0.771 |
| EBF1 | 0.719 | 0.678 | 0.682 | 0.725 | 0.706 | 0.734 | 0.731 | 0.729 |
| PAX5 | 0.756 | 0.626 | 0.659 | 0.708 | 0.655 | 0.674 | 0.697 | 0.702 |
| SPI1 | 0.820 | 0.673 | 0.689 | 0.745 | 0.721 | 0.722 | 0.746 | 0.753 |
| YY1 | 0.726 | 0.537 | 0.552 | 0.645 | 0.599 | 0.590 | 0.599 | 0.617 |

### 6.3 Per-TF prior TF log-odds

| TF | AG KL 10 | AG ol001 | `ag_sptf` | GELU B | GELU B KL 300 | `r2_full` | `hist100_190` |
|---|---|---|---|---|---|---|---|
| CTCF | 0.87 | 0.90 | 2.02 | 1.84 | 1.97 | 2.56 | 2.70 |
| NRF1 | 0.85 | 1.11 | 1.55 | 1.40 | 1.71 | 1.68 | 1.81 |
| USF2 | 0.30 | 0.42 | 0.79 | 0.65 | 0.99 | 1.08 | 1.26 |
| CREM | 0.29 | 0.35 | 0.60 | 0.77 | 0.78 | 0.87 | 0.89 |
| RFX5 | 0.32 | 0.22 | 0.43 | 0.57 | 0.43 | 0.75 | 0.85 |
| EBF1 | 0.20 | 0.20 | 0.40 | 0.33 | 0.48 | 0.54 | 0.56 |
| PAX5 | 0.20 | 0.21 | 0.42 | 0.38 | 0.35 | 0.57 | 0.51 |
| SPI1 | 0.25 | 0.32 | 0.48 | 0.48 | 0.49 | 0.66 | 0.73 |
| YY1 | 0.14 | 0.20 | 0.49 | 0.27 | 0.28 | 0.37 | 0.51 |

### 6.4 State-level diagnostics (Viterbi, every 8th val window)

| model | off-motif TF runs / 1,000 reads | kernel-supported | TF at its minimum | other_long / 1,000 reads | other_long at its cap | histone length |
|---|---|---|---|---|---|---|
| GELU B, KL 300 | 269 | 82% | 40% (8 bp) | 395 | 58% | 147 (fixed) |
| `r2_full` | 123 | 92% | 20% (15 bp) | 323 | 40% | 147 (fixed) |
| `hist100_190` | 116 | 86% | 21% (15 bp) | 130 | 9% | median 142 (100-190) |

### 6.5 Reading the progression
*[Prose for LaTeX.]*
- **Raw AUROC.** It rose from 0.700 (AG KL 10) to 0.772 with spacers + TF/3. After that it fell when
  the GELU head and the informative prior arrived (0.741), then recovered to 0.765 through the
  duration and nucleosome changes. It has not passed `ag_sptf` at ×1, and it remains below phase-type.
- **Specificity improved steadily.**
  - Off-motif TF per read-kb: 13.3 (`ag_nochain`) → 14.4 (GELU B) → 6.2 (KL 300) → 4.8 (`hist100_190`); phase-type 38.6.
  - Matched AUPRC: 0.475 (`ag_nochain`) → 0.525 (KL 300) → 0.535; phase-type 0.485.
- **Prior quality** (mean log-odds) went 0.38 → 0.80 → 1.09, with CTCF at 2.70.
- **Within-bound ranking** (ρ 0.111 → 0.215) and **accessible-read AUROC** (0.629 → 0.687) are now at
  or above phase-type.

---

## 7. Open issues and next steps
- **Nucleosome lengths follow the duration prior.** Viterbi median 142 vs a prior median of 147, with
  about the same spread. Options: a tighter scale (about 8 bp), or a "spike at 147 + broad tail"
  prior.
- **TF 15-60** (`hist190_tf60`) gives small gains, but TF barely uses 41-60 bp. A flatter TF duration prior would be needed to move the CTCF-like patches out of other_long.
- **Removing other_long.** The remaining 61-99 bp other_long calls look like sub-nucleosomal particles,
  so a histone minimum of about 60-70 bp may let other_long go. That conflicts with a tight length
  prior.
- **Open-state overdispersion.** Clustered protection in accessible DNA is unmodelled: runs of ≥ 8
  protected calls are 1.3-1.5× the open emission's expectation.
- **The informative-prior AUROC gap.** A constant prior scores AUROC 0.790-0.796, while informative
  priors score 0.74-0.77. This is unexplained.
- **Per-read labelling efficiency.** Off-motif TF calls concentrate in poorly labelled DNA, and per-read
  b_m is miscalibrated.

---

## 8. Figure plan for the LaTeX version
*[Copy each figure into `report/figures/`; all source files are on scratch.]*

1. Transition topology diagram: open hub, spacers, TF, other_long, histone with exit tail. *[to draw]*
2. m6A context bias: rate by context, held-out gain by context width (10-30 note figures).
3. Footprint screen and blind positions (11-10 note figures).
4. ReLU head collapse vs GELU: dead units, prior log-odds by seed (11-15 note figures;
   `fiber_eval_work/prior_log_odds/prior_log_odds_*.png`).
5. Viterbi length histograms through rounds 1-3 (`tf_mechanism/figures/viterbi_lengths_{kl,r2,r3}_proportion.png`).
6. Short other_long examples, isolated vs nucleosome-adjacent (2026-11-17 note figures).
7. Per-read posteriors at bound motifs for the lineage (`gm12878_validation/r2_ablation_panels/`;
   `evaluation.tf_chip.plot_model_comparison`).
8. Metric progression chart: AUROC, AUPRC, off-motif TF and prior log-odds by model
   (`fiber_eval_work/report_metrics/lineage_metrics.tsv`). *[to make]*

---

## Sources

- **Research notes:** `LM_research_notes/projects/smf-hmm/`, 2026-10-22 to 2026-11-17. The index is
  `README.md`.
- **Metric tables** in `/scratch/users/diamant/fiber_eval_work/`:
  - `idr_summary_*_val.tsv`;
  - `kernel_ref/specificity*/{auprc_matched,off_motif}.tsv`;
  - `prior_log_odds/val_*priors.tsv`.

  They are collected in `report_metrics/{lineage_metrics,idr_by_model,auprc_by_model,off_motif_by_model,prior_log_odds_by_model,per_tf_auroc,per_tf_prior_log_odds}.tsv`.
- **Viterbi diagnostics:** `fiber_eval_work/tf_mechanism/{viterbi_stats_r3.txt, scorecard_*.tsv, short_ol/}`.
- **Runs:** `/scratch/users/diamant/2026-10-22_fiber_fixed_blocks_global_kl3` (first) to
  `2026-11-17_fiber_hist100_190_kl300_alphagenome_bf16` (current). The header comments in the current
  `train.sbatch` files are stale; the flags and `.hydra/config.yaml` are authoritative.
- **Code:** smf_hmm branch `nd/experiments`.
  - **Committed:** see §3.1, plus the TF floor (a51fbef), spacers (d128fc1) and the edge floor
    (f0401ff).
  - **Uncommitted:** GELU head, `tied_blocks`, `ConstantGroupPrior`, the variable-length histone and
    `prior_tf_log_odds.py`.
