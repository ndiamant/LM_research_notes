---
title: Swiss-Prot INDIGO vs AR — memorization confounds, ELBO advantage, and a shared data ceiling
date: 2026-09-18
project: arid
agent: claude
status: draft
sources:
  - /home/users/diamant/repos/ARID/examples/swissprot_indigo_pretrain_finetune.py
  - /home/users/diamant/repos/ARID/examples/swissprot_indigo_sample.py
  - /home/users/diamant/repos/ARID/examples/swissprot_indigo_trajectory_stats.py
  - /home/users/diamant/repos/ARID/examples/swissprot_chunk_coherence.py
  - /home/users/diamant/repos/ARID/examples/swissprot_evaluate.py
  - /home/users/diamant/repos/ARID/arid/schedules.py
  - /home/users/diamant/repos/ARID/arid/proteins.py
  - /scratch/users/diamant/arid_runs/indigo_partition3_ar_pretrained
  - /scratch/users/diamant/arid_runs/indigo_partition3_random_init
  - /scratch/users/diamant/arid_runs/indigo_partition3_random_init-100epcs
  - /home/users/diamant/repos/ARID/outputs/chunk_coherence/results.json
  - /home/users/diamant/repos/ARID/outputs/chunk_coherence/results_novel_only.json
tags:
  - indigo
  - swissprot
  - autoregressive
  - memorization
  - evaluation
  - elbo
---

# Summary

On Swiss-Prot at matched size (20.4M parameters) and data (174k train
sequences), tokenwise INDIGO with a random three-chunk partition is a **better
density model than the architecture-matched GPT-2 AR baseline** and a **worse
sampler**. Its ELBO on held-out proteins is 2.287 nats/residue against AR's
exact 2.383, while its novel-sample ESM-C pseudo-perplexity is 9.71 against
AR's 8.55.

Four comparisons in this investigation were confounded by memorization and
reversed or shrank substantially once it was controlled. AR emits 10.0% exact
training copies and 25.1% near-copies; INDIGO emits 1.6% and 10.3%. Because a
memorized sample *is* a real protein, ESM-C scores it like one, so raw PPPL,
k-mer JS divergence, and the block-coherence statistic all reward copying. The
apparent AR advantage roughly halves after filtering (0.406 → 0.183
nats/residue on whole sequences).

Both models saturate at val ≈ 2.35 nats/residue with train–val gaps near 0.9,
and 4x the training epochs bought zero novel-sample quality. The binding
constraint appears to be **data, not formulation or capacity** — though this is
inferred from generalization gaps, not yet measured; data-scaling runs are in
flight.

A separate structural finding: with three chunks, **97.8% of pointer steps are
trivial continuations**, so `position_acc` is uninformative (99.12% is
achievable while getting every free choice wrong) and the partition resampling
provides much weaker augmentation than a true any-order objective would.

# Key Points

- **INDIGO wins on held-out likelihood.** ELBO 2.287 nats/residue (25 epochs)
  vs AR's exact 2.383. The ELBO is an upper bound, so the true marginal is
  better still. Roughly 15 nats/sequence.
- **AR wins on sampling, by half as much as it first appeared.** Novel-only
  PPPL 8.55 vs 9.71; the unfiltered numbers were 6.46 vs 8.30.
- **Memorization is the dominant evaluation confound.** AR: 10.0% exact, 25.1%
  near-duplicate, mean 20-mer containment 0.241. Held-out real proteins sit at
  0.028 containment, so containment is a usable memorization scale.
- **The 100-epoch run's apparent improvement was entirely memorization.**
  All-sample PPPL improved 9.50 → 8.71, but novel-only PPPL was flat
  (10.50 → 10.36) while exact copies rose 0.6% → 1.8% and near-duplicates
  5.8% → 10.4%.
- **Both formulations hit the same data ceiling.** AR reaches a 0.94 train–val
  gap in 25 epochs; INDIGO takes 100 epochs to reach 0.90. The any-order
  resampling slows memorization roughly 4x but does not prevent it.
- **INDIGO's block coherence is genuinely worse than AR's**, and this survives
  novelty filtering: context benefit Δ = +0.011 vs AR's +0.063, against bounds
  of +0.210 (real) and −0.130 (chimeras). On that scale INDIGO is 41% and AR
  57% of the way from chimeric to coherent.
- **The pointer head is at its information-theoretic floor** and is only 2.8%
  of the objective, so Rao-Blackwellizing it is not worth doing under this
  schedule.
- **The models' generation order matches the training schedule closely**
  (1.09–1.12 jumps/sequence vs the schedule's 1.16), so the generation gap is
  not caused by off-distribution trajectories.

# Details

## Runs

All runs use `examples/swissprot_indigo_pretrain_finetune.py`, seed 7, batch 64,
lr 3e-4, `max_content_length` 256, the random three-chunk partition hardcoded in
`ProteinIndigoCollator`, and the homology-aware split at
`/scratch/users/diamant/datasets/swiss_prot/splits/splits.json`
(174,015 train / 12,251 val).

```bash
# AR-pretrained INDIGO (25 + 25 epochs)
python examples/swissprot_indigo_pretrain_finetune.py \
  --pretrain-epochs 25 --finetune-epochs 25 \
  --output-dir /scratch/users/diamant/arid_runs/indigo_partition3_ar_pretrained \
  --wandb --wandb-name indigo-partition3-ar-pretrained \
  --wandb-tags partition3_random ar_pretrained pretrain25 finetune25

# random init (0 + 25), and the 100-epoch extension
python examples/swissprot_indigo_pretrain_finetune.py \
  --pretrain-epochs 0 --finetune-epochs 100 \
  --output-dir /scratch/users/diamant/arid_runs/indigo_partition3_random_init-100epcs \
  --wandb --wandb-name indigo-partition3-random-init-100epcs \
  --wandb-tags partition3_random random_init finetune100
```

**The INDIGO finetune phase requires an 80 GB GPU.** On a 48 GB card
FlexAttention fails in `indigo_flex_attention_forward` with
`OutOfResources: shared memory, Required: 184352, Hardware limit: 101376`.
Use `-C "GPU_MEM:80GB"`, not the `48GB|80GB` constraint copied from the
reverse-edit submit scripts. Finetune costs ~4 min/epoch, pretrain ~2.1 min.

## Teacher-forced losses

| run | train span | val span | gap | val position |
|---|---|---|---|---|
| AR GPT-2, 25ep | 1.430 | 2.368 | +0.94 | — |
| INDIGO ar_pretrained, 25ep | — | 2.306 | — | 0.067 |
| INDIGO random_init, 25ep | 1.700 | 2.292 | +0.59 | 0.066 |
| INDIGO random_init, 100ep | 1.450 | 2.353 | +0.90 | 0.067 |

AR's `val/loss` is per token over N+1 targets (residues + EOS); INDIGO's
`val/span_loss` is per residue. INDIGO's pooled `val/loss` (~1.18) is **not**
comparable to AR's 2.368 — it averages over ~2N+1 decisions (one value and one
pointer per residue, plus DONE), so it is roughly half the per-residue NLL by
construction.

## The ELBO comparison

Converting both to nats per sequence over the val set (12,251 sequences, mean
length 152.7, mean order entropy H(q) = ln C(L−1,2) + ln 6 = 10.92 nats):

| | nats/sequence | per residue |
|---|---|---|
| AR, exact −log p(x) | 363.9 | 2.383 |
| INDIGO 25ep, joint −log p(x, order) | 360.1 | — |
| INDIGO 25ep, ELBO on −log p(x) | **349.2** | **2.287** |
| INDIGO 100ep, joint | 369.5 | — |
| INDIGO 100ep, ELBO | 358.6 | 2.348 |

The ELBO bounds INDIGO against *itself*, not against AR, so INDIGO's joint
falling below AR's exact NLL is not a violation — it means INDIGO assigns
held-out proteins higher probability.

Notable: the pointer term (10.1 nats/sequence) nearly cancels the entropy
correction (10.92), so the ELBO per residue ≈ `span_loss`. *Caveat:* H(q) is an
**upper** bound on the true order entropy, because the trajectory→latent map is
not injective (a left-to-right permutation produces the same pointer sequence
for every cut choice). The correction therefore over-credits INDIGO by up to
~1 nat/sequence out of a 14-nat margin.

## Memorization and the novelty-filtered evaluation

T = 1.0, 1024 samples per source, ESM-C 600M. Novelty = not an exact train
match **and** ≤50% of the sample's 20-mers present in a 20-mer index built from
50k training sequences. PPPL columns are 384-sample subsets, so absolute levels
are noisier than the full 1024-sample evaluations and should not be mixed with
them; within-table comparisons are consistent.

| source | exact% | containment | near-dup% | PPPL all | PPPL novel |
|---|---|---|---|---|---|
| real (val) | 0.0 | 0.028 | 1.2 | 2.827 | 2.833 |
| AR model | 10.0 | 0.241 | 25.1 | 6.463 | **8.553** |
| INDIGO ar_pretrained 25ep | 1.6 | 0.101 | 10.3 | 8.298 | **9.713** |
| INDIGO random_init 25ep | 0.6 | 0.061 | 5.8 | 9.497 | **10.501** |
| INDIGO random_init 100ep | 1.8 | 0.105 | 10.4 | 8.705 | **10.355** |

Exact match badly undercounts: AR is 10% exact but 25% near-duplicate. The
composition metrics are confounded in the same direction — AR's `kmer3_js`
*worsens* from 0.0417 to 0.0462 once copies are removed, because memorized
samples are real proteins with perfect k-mer statistics.

`examples/swissprot_evaluate.py` now reports every metric for both subsets.
**Breaking format change:** `metrics.json` is now
`{"config":…, "metrics": {"all": {...}, "novel": {...}}}`, so existing tooling
reading `metrics.<key>` needs `metrics.all.<key>`. Evaluation directories
written before 2026-09-18 are in the old flat format.

```bash
python examples/swissprot_evaluate.py --no-esmc \
  --sequences-file <run>/generations_temp_1.0/generated.txt \
  --splits /scratch/users/diamant/datasets/swiss_prot/splits/splits.json \
  --output-dir outputs/<run>/evaluation
```

AR pretraining's benefit survives the correction (novel PPPL 9.71 vs 10.50) and
is worth roughly 75 extra finetune epochs: the 100-epoch random-init run reaches
the 25-epoch pretrained run's all-sample quality.

## Chunk coherence

`examples/swissprot_chunk_coherence.py` measures the context benefit

```text
delta = PLL(whole sequence) − PLL(its blocks scored separately)
```

per residue over the same residues, against two bounds: real proteins cut the
same way (coherent) and blocks taken from three *different* length-matched real
proteins (chimeric). Blocks for INDIGO come from the model's own recorded
trajectory; for other sources they are sampled from the same schedule.

| source | Δ all | Δ novel only | 95% CI (novel) |
|---|---|---|---|
| real proteins | +0.210 | +0.210 | [0.190, 0.230] |
| chimeras | −0.130 | −0.130 | [−0.145, −0.115] |
| AR model | +0.093 | **+0.063** | [0.047, 0.079] |
| INDIGO ar_pretrained | +0.011 | **+0.011** | [0.002, 0.021] |
| INDIGO random_init | +0.001 | **+0.000** | [−0.008, 0.009] |

Memorization inflated AR's coherence by a third; INDIGO's was unaffected. On
novel samples the whole-sequence gap decomposes as 0.132 nats (72%) from blocks
being worse *in isolation* and 0.051 nats (28%) from worse coherence — so local
fragment quality dominates, and the original "fragments fine, assembly bad"
hypothesis is only partly right.

*Caveats:* the real and chimera rows are not novelty-filtered (both derive from
held-out val under a homology-aware split; val measures 1.2% near-duplicate).
The chimera bound is not perfectly level-matched — its blocks score −1.338
against the real set's −1.197 despite both being real fragments — but Δ is a
within-sequence difference, so the offset largely cancels.

AR reaching only 57% of the way from chimeric to coherent suggests a substantial
part of this is what 20M-parameter samples look like at this data scale, with
INDIGO worse rather than qualitatively broken. *This interpretation is
speculative.*

## Why the pointer head is uninformative under this schedule

Measured on 200 val proteins by replaying the collator's own trajectories:

| pointer step kind | share |
|---|---|
| continuation (target = one slot right of the last insertion) | 97.79% |
| chunk start with a free choice of gap | 0.88% |
| forced (step 1, only one gap exists) | 0.67% |
| DONE | 0.67% |

A model that always continues and never learns chunk placement scores 99.12%.
Measured `position_acc` is 98.6%. With three chunks there are at most two
informative decisions per protein, and adjacency collapses that to ~1.3.

The pointer is also near its information-theoretic floor: `position_loss` is
0.067 nats/step against a chain-rule bound of 0.071 (itself an over-estimate of
the true floor), and the head is only 2.8% of the pooled objective. **Therefore
Rao-Blackwellizing the pointer loss — replacing the one-hot target with the
exact posterior over next insertion positions — is not worth doing under the
3-chunk schedule.** The recoverable KL is 0.002–0.027 nats/step on 2.8% of the
loss. It becomes worthwhile only with schedules that have real order entropy
(more chunks, or a Markov jump process).

## Ruled out

- **Early stopping / bad checkpoint.** No `EarlyStopping` callback exists.
  `write_generations` samples the final in-memory model, never the best
  checkpoint, so generations always reflect the last epoch. In the 100-epoch run
  best val/loss was epoch 56 (1.1960) vs 1.2059 at epoch 100 — a real but small
  0.010 nats difference.
- **Dataloader failing to resample trajectories.** Verified empirically: 31/32
  sequences receive a different trajectory between consecutive epochs, and 4.66
  distinct trajectories per sequence over 5 epochs (max 5). The persistent
  worker's cached generator keeps advancing, as intended.
- **Off-distribution generation order.** Model trajectories match the schedule:
  1.09 / 1.12 jumps per sequence against 1.16, blocks 2.27 / 2.30 against 2.32.
  Only 3 sequences in 1024 per model exceed 2 jumps.
- **Unfinished LR annealing in the 100-epoch run.** The cosine reaches
  1.1e-5 by epoch 90 and 3.0e-6 at epoch 100, the same floor as the 25-epoch run.

## Why the any-order augmentation is weak here

*Inferred mechanism, not directly measured.* Resampling the partition changes
*which context* a residue is predicted from but never *what* it predicts: at a
continuation step the target is "the residue after this one", identical to AR's
target. All epochs present the same ~33M residue-prediction targets, and only
~1% of steps (chunk starts) receive genuinely new targets. A second-order effect
may push the wrong way: the extra context is a *retrieval cue* that helps
identify which training protein is being reproduced.

This contrasts sharply with the synthetic island task, where the position head
is ~68% of the objective (see related note below). The difference is the chunk
count: 3 contiguous chunks over a 150-residue protein is nearly left-to-right,
whereas the synthetic schedule inserts throughout.

# Related Notes

- [INDIGO autoregressive pretraining status and Swiss-Prot plan](../2026-09-16/indigo-ar-pretraining-and-protein-plan.md):
  This note executes that plan's proposed next experiment — randomly
  initialized vs AR-pretrained INDIGO on Swiss-Prot at matched capacity. Its
  synthetic prediction that pretraining helps the value head is confirmed at
  protein scale, though the margin is small.
- [Cost of any-order INDIGO generation, and learned deletions via tombstones](indigo-any-order-cost-and-tombstone-deletions.md):
  Same-day synthetic results. Its headline finding that order entropy dominates
  the objective (position head ~68%) is the opposite of what the 3-chunk protein
  schedule shows (2.8%), and the chunk count explains the difference. Its
  entropy-probe method is what makes the pointer-floor calculation here cheap.
- [Swiss-Prot temperature-1 results and append-3 supervision scaling](../2026-09-15/swissprot-temperature-1-and-supervision-scaling.md):
  Establishes the T=1.0 evaluation convention used throughout.
- [Swiss-Prot partition-schedule results](../2026-09-10/swissprot-partition-schedule-results.md):
  Source of the `partition_3_random` schedule these runs match.

# Open Questions

- **Is the data ceiling real?** Two matched-compute data-scaling runs are in
  flight (43.5k sequences × 100 epochs and 87k × 50 epochs, both ~86k optimizer
  steps, matching the existing full-data 25-epoch run). If val loss is still
  descending steeply from 87k → 174k, the data-limited reading is confirmed.
- **Would a smaller model generalize better?** If genuinely data-limited at
  20.4M parameters, halving `d_model` should *improve* val loss. Untested, and
  independent of the scaling runs.
- **Does more data close the coherence gap?** Δ is a property of the generation
  process, not obviously of scale. INDIGO's 41% vs AR's 57% may be structural.
- **Why are INDIGO's blocks worse in isolation (72% of the novel gap)?** Within
  a block INDIGO performs the same left-to-right continuation AR does, and ties
  AR on teacher-forced loss, yet its sampled blocks score 0.132 nats worse. This
  is exposure bias *within* a fragment and is unexplained.
- **Is distance-blindness responsible?** The relative-position encoding is
  `sign(distance) + 1` — three codes, no distance — and GPT-2 absolute positions
  index insertion time, not sequence position. So a gap's right neighbour sits
  at unknown range. `relative_position_embeddings="distance_buckets"` already
  exists and would test this in one run. *Hypothesis, untested.*
- **Should model selection use val loss at all?** Across these runs val loss and
  sample quality moved in opposite directions (mediated by memorization), which
  makes `ModelCheckpoint(monitor="val/loss")` a questionable selection rule.

# Sources

- Training entry point: [`examples/swissprot_indigo_pretrain_finetune.py`](https://github.com/ndiamant/ARID/blob/main/examples/swissprot_indigo_pretrain_finetune.py)
  (wandb logging and a per-epoch `GenerationCallback` added 2026-09-17; finetune
  metrics log under bare `train/`/`val/` so they share panels with the
  reverse-edit runs, pretraining under `pretrain_*`)
- Re-sampling at arbitrary temperature, with trajectory recording:
  `examples/swissprot_indigo_sample.py`
- Generation-order statistics: `examples/swissprot_indigo_trajectory_stats.py`
- Coherence measurement: `examples/swissprot_chunk_coherence.py`
- Shared trajectory helpers (`insertion_positions_from_partition`, `chunk_ids`,
  `blocks_from_positions`) in `arid/schedules.py`; k-mer and containment helpers
  in `arid/proteins.py`
- Run artifacts under `/scratch/users/diamant/arid_runs/indigo_partition3_*`
- Coherence results: `outputs/chunk_coherence/results.json` and
  `results_novel_only.json`
