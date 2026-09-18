---
title: Cost of any-order INDIGO generation, and learned deletions via tombstones
date: 2026-09-18
project: arid
agent: claude
status: draft
sources:
  - /home/users/diamant/repos/ARID/arid/indigo.py
  - /home/users/diamant/repos/ARID/arid/island_task.py
  - /home/users/diamant/repos/ARID/arid/island_plots.py
  - /home/users/diamant/repos/ARID/arid/schedules.py
  - /home/users/diamant/repos/ARID/examples/indigo_insertion_synthetic.py
  - /home/users/diamant/repos/ARID/examples/indigo_tombstone_synthetic.py
  - /home/users/diamant/repos/ARID/examples/tombstone_entropy_probe.py
  - /home/users/diamant/repos/ARID/tests/test_indigo.py
  - /home/users/diamant/repos/ARID/outputs/indigo_d192_L8
  - /home/users/diamant/repos/ARID/outputs/indigo_mdl_1616
  - /home/users/diamant/repos/ARID/outputs/indigo_mdl_88
  - /home/users/diamant/repos/ARID/outputs/indigo_tombstone_1x1
  - /home/users/diamant/repos/ARID/outputs/indigo_tombstone_mix
tags:
  - indigo
  - any-order
  - tombstones
  - deletions
  - synthetic
  - instrumentation
---

# Summary

Two results on the binary-island synthetic task. First, the cost of any-order
INDIGO generation is dominated by *order entropy*, not by model capacity: the
position head accounts for about 68% of the objective, and a left-to-right
control reaches 0.945 validity in 10 epochs where the any-order model needs
thousands. Second, learned deletions via tombstones work: the trajectory stays
append-only, one-shot log probabilities still match cached sampling, and the
model uses the delete action at close to the correct rate, lifetime, and count.

The tombstone cost was predicted before training. An entropy probe forecast
+25% objective entropy for one junk token; the measured pointer loss came in at
+25.1%. That probe is now the cheap way to price schedule changes.

The largest single win came from a training-distribution fix rather than the
model. With junk in every trajectory, 15.3% of samples took an off-distribution
shortcut and were only 0.229 valid. Putting zero-junk trajectories into the
mixture at 50/50 raised those rows to 0.840 and overall validity from 0.711 to
0.789 (best 0.828, still rising at 500 epochs).

Three interventions were tried and are reported honestly: the gap-window
pointer helped (1.85x lower hazard), scaling helped modestly (1.53x at 7.5x
parameters, with no generalization gap, so capacity was partly binding), and a
left-to-right-to-any-order curriculum **actively hurt** through a diagnosable
mechanism. The framing that survives all of it: maximum likelihood reproduces
the training order distribution and does not discover executable orders.

# Key Points

- Order entropy, not capacity, is the dominant cost. At the default
  `(min,max)=(1,4)` schedule, position loss is about 68% of the objective, and
  its target is drawn independently of the data, so most of it is irreducible.
- Any-order generation must fit a much larger conditional than the data
  suggests: the task has 8 targets and 108 left-to-right-reachable states, but
  **1289** insertion-reachable states (any subsequence of a target), a 12x
  larger domain with far higher-entropy targets.
- `min_delete_length` was added to `DeletionOnlySchedule`. Setting
  `min = max = sequence_length` pins every trajectory to one whole-sequence
  span, which reverses into *exact* left-to-right generation (verified: pointer
  targets are `1..16` then DONE on every draw, 1 unique order in 2000 samples,
  0.00 bits of order entropy). This is the true AR control; `max=16` alone is
  not, because span lengths are drawn uniformly.
- Measured conditional order entropy H(pos | prev) across the ladder:
  `(1,1)` 3.40, `(2,2)` 2.16, `(1,4)` 1.95, `(4,4)` 1.05, `(8,8)` 0.30,
  `(16,16)` 0.00 bits. Unique orders per 2000 samples: 2000, 1998, 1999, 567,
  9, 1.
- Models reproduce their training order entropy almost exactly and show no
  concentration onto orders they execute well: `(1,4)`-trained models
  self-generate at H=1.97 with all 2048 sampled orders distinct, while failing
  8-13% of the time. **MLE cannot discover executable orders**; an objective
  that can prefer them (rejection-sampling finetuning, RL on validity) would be
  needed for a learned-order claim.
- The pointer head was information-starved: it keyed each gap by its left
  neighbour alone. Keying by both bounding slots (`--pointer-context
  gap_window`) gave 1.85x lower off-manifold hazard. Because the key is
  *additive*, the score decomposes into gathers over scalar tables, which also
  cut the pointer's one-shot activation memory from O(L^2 * d) to O(L^2)
  (measured 3.8x less at L=256, and the ratio grows with L).
- Scaling from 0.51M to 3.79M parameters (d96/L4 to d192/L8) gave +4.6 validity
  points and 1.53x lower hazard, with *no* train/val gap at either size, so
  capacity was partly binding. Width is nearly free in wall clock here; only
  depth costs time.
- Curriculum learning (anneal span floor from sequence length down to the
  target) **hurt monotonically**: hazard rose 16%/27%/41% for 50/100/200
  curriculum epochs. Mechanism identified below.
- Tombstones preserve the invariant they exist for. A deletion appends a marker
  instead of removing anything, retired tokens keep their order slot so settled
  sign relations never change, and a parity test confirms one-shot scoring
  matches cached sampling on deletion trajectories to ~1e-7.
- The delete head is **ambiguity-limited, not supervision-limited**: halving its
  supervision (50/50 mixture) left `delete_recall` at 0.212 versus 0.216 and
  `delete_acc` at 0.106 versus 0.108.

# Details

## Instrumentation built first

Prior comparisons used 32 samples for validity, giving +-7% binomial error --
too noisy to resolve the effects being chased. The following now exist and are
prerequisites for trusting any of the numbers above:

- **Batched sampler** (`IndigoInsertionModel.sample_batch`). State is kept in
  chronological coordinates with a rank tensor recovering the current layout;
  every row is stepped every iteration with finished rows fed padding, keeping
  buffers rectangular. 512 samples now cost 1.29s versus 6.16s for 32
  one-at-a-time, roughly 76x per sample. `VALIDITY_SAMPLES` raised to 512.
- **Off-manifold hazard** (`off_manifold_mass`). Insertions only add, so a
  partial state can still reach a target exactly when it is a subsequence of
  one -- an exact, cheap oracle. The metric sums sampler probability over
  actions that leave the manifold, counting only decisions from states that are
  still reachable and not yet complete. Relative standard deviation 1.8% versus
  validity's 4.5% at 512 samples (29.4% at 32).
- **Diversity** (`target_tv_distance`) as total variation between the data
  distribution and the valid samples, reported beside its perfect-sampler
  floor. The floor is large and non-obvious: 8 targets at 512 samples reads
  0.046 even for a flawless model. Analytic form validated against Monte Carlo
  to three decimals.
- **Temperature default moved 0.8 to 1.0.** At 0.8 the hazard was understated
  1.7-1.8x, because the metric consumes the temperature-adjusted distribution.
  All earlier validity numbers were flattered; comparisons to entropy floors at
  T<1 are a category error.
- `run_config.json` written per run, since nothing previously recorded which
  temperature or span setting produced an output.
- Batched loader support (`--num-workers`). The collator was ~60% of epoch wall
  clock. Four workers plus an algorithmic fix took epochs from 6.0s to 1.8s.

Two bugs found by this instrumentation are worth recording. The collator held a
seeded `torch.Generator` built in the parent process, so all loader workers
replayed the *same* deletion patterns (4 distinct trajectories in 16 instead of
16) -- silent augmentation collapse. And `append_relative_position` built the
full k x k relation matrix and discarded all but the last row, making trajectory
construction cubic in length; computing the row directly halved collator cost
and matters much more at protein lengths.

## Loss is a poor proxy for sample quality, in both directions

Across the 5000-epoch runs, validation loss moved only 0.035-0.040 nats from
epoch 200 to 5000 while validity went 0.74 to 0.98. Position accuracy sat at
0.733-0.743 throughout. Conversely, a single-epoch training spike in the
`(8,8)` run (train loss 0.1302 to 0.1437, +10%) collapsed validity from 0.994
to 0.639 and raised hazard 56x, recovering fully within two epochs.

Consequences: (1) run comparisons must use tail averages, not final-epoch
values -- that run would have recorded 0.64 for a model actually at 0.998; (2)
the hazard and TV metrics caught an event that loss barely registered.

## Why the curriculum failed

`val/position_loss` is measured against the target distribution throughout, so
it shows what the left-to-right phase does to fitness for the real objective:

| epoch | span floor | curr200 pos_loss | control pos_loss | curr200 validity |
|---:|---:|---:|---:|---:|
| 1 | 16 | 2.51 | 1.55 | -- |
| 25 | 11 | **6.23** | 1.05 | -- |
| 50 | 8 | 5.82 | 0.95 | **0.826** |
| 100 | 4 | 3.45 | 0.87 | **0.082** |
| 200 | 1 | 0.90 | 0.86 | 0.387 |
| 500 | 1 | 0.85 | 0.84 | 0.834 |

A uniform random pointer over ~17 candidates scores about 2.8 nats. The
curriculum drove it to 6.23 -- the head is not uninformed but *confidently
wrong*, having learned "always insert at the right end", which assigns near-zero
probability to the positions the target distribution requires. It takes ~150
epochs to unlearn. Note epoch 50: validity 0.826 with target position loss
5.82, i.e. an excellent left-to-right generator that cannot do any-order.

This is the failure mode InDIGO's own paper warns about ("the position
prediction module learns much faster ... and quickly latches onto spurious
correlations"), except InDIGO used left-to-right pretraining as the *cure* for
beam-searched orders, a narrow target distribution. Against a broad random-order
target it becomes the cause. It is also consistent with sigma-GPT's curriculum
being safe only because its order is an *input*, so there is no position head to
corrupt. **Inferred, not proven:** the mechanism is consistent with the data but
was not isolated by an intervention.

Budget-adjusted, the curriculum was still behind at matched epochs-at-target
(50 target epochs: 0.111 versus the control's 0.244), so it is not merely spent
budget.

## Tombstone design and the entropy probe

`examples/tombstone_entropy_probe.py` prices a schedule before training. For
one junk span of length one at random forward placement it forecast +25%
objective entropy; the trained model's pointer loss came in at +25.1% (1.0235
versus the insertion-only 0.8182). Other probe findings:

- Placement of the forward insertion controls the *reverse lifetime* of a junk
  token: 9.2 steps at placement 0.0 versus 2.8 at 0.9. Randomising placement
  costs only ~2 bits over the cheapest fixed placement but produces a wide
  lifetime spread (5.8 +- 4.1), so a random-placement run cannot isolate whether
  long lifetimes are the failure mode.
- Span source barely matters **on this task**: donor spans copied from another
  sample cost 2.75 bits at length 3 versus 3.00 uniform, under 1% of the total.
  The binary alphabet leaves no headroom (uniform 1.00 bit, data marginal
  0.954). **This task cannot answer the span-source question**; proteins, with
  a 20-token alphabet, plausibly would.
- Multiple junk spans are *sub-additive*: 11.09 bits per junk token at k=1
  falling to 9.02 at k=4. Longer spans are cheaper per deleted token (8.07 at
  span 2) because "where is this span" is paid once per span. Six junk tokens
  against 16 real roughly doubles the objective.
- The delete target is substantially ambiguous: 4-11 of the available single
  token deletions keep a target reachable. Much of the delete entropy therefore
  sits on distinctions that do not affect validity -- the same pattern as order
  entropy.

Forward coordinates were chosen over specifying the reverse trajectory
directly. With k spans of variable length, forward coordinates make trajectory
validity free by construction (run edits until the sequence is empty), whereas
reverse coordinates would require hand-guaranteeing that every junk token
inserted is later deleted and that the trajectory ends at the clean data.
Lifetime is still controllable in forward coordinates, since it equals
`t_delete - t_insert`.

## Tombstone results

Both runs are d192/L8, 500 epochs, T=1.0, `max_span=1`, donor spans, random
placement.

| | always-1-junk | 50/50 mixture | insertion-only d192/L8 |
|---|---:|---:|---:|
| validity (tail-10) | 0.7107 | **0.7891** | 0.9053 |
| best validity | -- | 0.8281 (ep 500) | -- |
| `rows_deleting` | 0.847 | 0.447 | -- |
| `delete_recall` | 0.216 | 0.212 | -- |
| `delete_acc` | 0.108 | 0.106 | -- |
| position loss | 1.0235 | 0.9658 | 0.8182 |
| span loss | 0.4502 | 0.4399 | 0.4096 |
| TV excess | -0.003 | -0.001 | +0.002 |

The mechanics are learned nearly exactly: 1.00-1.01 deletes per deleting row
(training 1.00), retired-token lifetime 5.26-5.47 (training 5.8 +- 4.1), 99.9%
of deleting rows using the exact training trajectory shape, and full diversity
throughout.

**The two-population finding.** In the always-1-junk run, splitting samples by
whether they deleted:

| | share | edits | validity |
|---|---:|---:|---:|
| retired its junk | 84.7% | 18 (100%) | 0.798 |
| skipped the junk | 15.3% | 16 (99.7%) | **0.229** |

The 15.3% were not failed deletions but a clean shortcut: 16 insertions, no
junk, correct final length. That trajectory shape appears in **zero** training
trajectories, so the model was extrapolating off-distribution. This is exposure
bias in the *trajectory shape* dimension rather than the usual state dimension.

**The fix.** `JunkConfig` now takes `span_counts` and
`span_count_probabilities`, mirroring `ChunkPartitionSchedule`'s chunk counts,
defaulting to `(0, 1)` uniform. After the fix the populations converged and the
ordering flipped, non-deleting rows now being the better population because a
no-junk trajectory has 16 decisions to get right instead of 18:

| | before | after |
|---|---:|---:|
| deleted at least once | 0.798 | 0.729 |
| never deleted | 0.229 | **0.840** |

**A predicted risk that did not materialise.** Halving the delete supervision
was expected to degrade the delete head. It did not move (`delete_recall` 0.212
versus 0.216), which is consistent with the probe's ambiguity finding: the head
is limited by genuine ambiguity in the target, not by example count. The
planned 0.25/0.75 tilt arm was therefore dropped as uninformative.

**Action-space overhead.** The probe predicts +12.5% for a 50/50 mixture (half
of +25%), but observed pointer loss is +18% over insertion-only. The excess is
plausibly the action-space expansion itself: on a zero-junk trajectory the
pointer must still rule out every delete candidate. **Speculative**, but if
real it is a fixed cost for *having* the capability that scales with live
sequence length, and so matters more at protein scale.

## What this says about the pretraining-compute question

Framing "any-order learns slower than AR" as one number conflates three things,
and only the third is a defect worth fixing:

1. **Order-entropy tax** (position loss 0.82 versus 0.00 nats). Not a
   deficiency -- it is the price of a latent order, and it is what is being
   bought.
2. **Reduced conditioning** (span loss 0.41 versus 0.13, about 32% worse
   per-token perplexity). Predicting from two gap neighbours is intrinsically
   harder than from a full prefix, and is exactly what makes the model useful
   for infilling.
3. **Sample efficiency** -- the extra epochs. The only real defect.

Because (1) and (2) are structural, any-order will never match AR on loss or
validity however well optimised, so "close enough to AR" cannot mean "same
loss".

This synthetic task also **overstates** the penalty by construction. Data
entropy is ~0.19 bits/token across 8 targets, so the order tax is 68% of the
objective; protein per-residue entropy of 2-3 nats would push the same tax to
roughly 25%. And all-or-nothing validity over 16 tokens maximally punishes
compounding error, where real objectives are per-token and locally redundant.
The >50x epoch multiplier seen here is close to a worst case. sigma-GPT's real
text numbers (30.43 versus 18.14 perplexity, largely closed by curriculum) are
a better guide to the actual penalty.

# Related Notes

- [INDIGO autoregressive pretraining status and Swiss-Prot plan](../2026-09-16/indigo-ar-pretraining-and-protein-plan.md): Uses the same `min/max-delete-length` schedule knobs recorded here, and independently observes the early pointer penalty from left-to-right pretraining that the curriculum failure below explains mechanistically. Its Swiss-Prot plan is the natural consumer of the instrumentation and cost decomposition in this note.
- [Swiss-Prot temperature-1 results and append-3 supervision scaling](../2026-09-15/swissprot-temperature-1-and-supervision-scaling.md): Establishes the temperature-1 evaluation protocol that this note's default-temperature change brings the synthetic scripts in line with.
- [Swiss-Prot partition-schedule results](../2026-09-10/swissprot-partition-schedule-results.md): Defines the random three-partition baseline that the protein version of the tombstone schedule would have to match.

# Open Questions

- Is the residual ~0.82 nats of position loss at the irreducible entropy floor,
  or is there extractable signal? The exact-posterior diagnostic (train against
  the true next-action distribution instead of a sampled one-hot) would settle
  whether the limit is gradient variance or representation, and bounds the
  payoff of every variance-reduction idea. Harder than first estimated: the
  deletion process has memory, so the posterior requires marginalising over
  deletion trajectories rather than tabulating 1289 content states.
- Does the action-space overhead (+18% observed versus +12.5% predicted for the
  50/50 mixture) really come from ruling out unused delete candidates? If so it
  grows with live sequence length and matters at protein scale.
- Is self-correction real? The model currently retires junk *it generated
  because the training distribution told it to* -- a planned insert-then-retire,
  not error detection. Testing the capability needs a probe asking whether
  deletions preferentially remove tokens that made the state unreachable;
  `reachable()` already supports this.
- How should the off-manifold hazard be redefined now that deletions are
  learned? Reachability became budget-dependent, and there is no useful
  budget-free version (the empty subsequence always embeds), so the
  low-variance metric that carried every insertion-only comparison is currently
  unavailable for tombstone runs.
- The 50/50 mixture had not converged at 500 epochs (0.719 to 0.828 across
  epochs 310-500). How much of the remaining 8-12 point gap to insertion-only
  closes with budget?
- Does the span-source choice (donor versus uniform) matter at protein alphabet
  size, where it has room to, and does it change *what* deletions the model
  learns for miniaturisation?

# Sources

- Implementation and tests listed in the frontmatter. 76 tests pass in the
  `arid_indigo` conda environment as of 2026-09-18, including a parity test
  asserting one-shot pointer distributions match cached sampling on deletion
  trajectories, and a multi-worker test asserting the curriculum epoch counter
  crosses the process boundary.
- Run outputs under `/home/users/diamant/repos/ARID/outputs/` as listed in the
  frontmatter. Note that `indigo_500_epochs_*` are 5000-epoch runs despite the
  directory names, and that all runs predating 2026-09-15 logged metrics at
  temperature 0.8; checkpoint-level numbers in this note were re-measured at
  T=1.0 for comparability.
- Entropy figures from `examples/tombstone_entropy_probe.py`, tabulated against
  a short context and therefore upper bounds on the true conditional entropy;
  the cross-arm comparisons are the reliable part.
- InDIGO (Gu, Liu & Cho, 2019), <https://arxiv.org/abs/1902.01370> -- the
  bootstrapping/position-module warning.
- sigma-GPT (2024), <https://arxiv.org/html/2404.09562v1> -- order as an input,
  and the curriculum that works when there is no position head.

## Command provenance

```bash
# Environment used throughout.
P=/home/groups/btrippe/diamant/miniforge/envs/arid_indigo/bin/python
export HF_HOME=/scratch/users/diamant/models

# Insertion-only controls (scaling comparison; d96/L4 is the default size).
$P -m examples.indigo_insertion_synthetic \
  --absolute-position-embeddings gpt2 --relative-position-embeddings indigo \
  --pointer-context gap_window --max-epochs 500 --lr 3e-4 --num-workers 4 \
  --output-dir outputs/indigo_gap_window
$P -m examples.indigo_insertion_synthetic \
  --absolute-position-embeddings gpt2 --relative-position-embeddings indigo \
  --pointer-context gap_window --min-delete-length 1 --max-delete-length 4 \
  --d-model 192 --nhead 8 --num-layers 8 --dim-feedforward 768 \
  --max-epochs 500 --lr 3e-4 --num-workers 4 --output-dir outputs/indigo_d192_L8

# Order-entropy ladder. min=max=16 is the exact left-to-right control.
$P -m examples.indigo_insertion_synthetic \
  --absolute-position-embeddings gpt2 --relative-position-embeddings indigo \
  --pointer-context gap_window --min-delete-length 16 --max-delete-length 16 \
  --max-epochs 500 --lr 3e-4 --num-workers 4 --output-dir outputs/indigo_mdl_1616

# Curriculum arms (negative result).
for curr in 50 100 200; do
  $P -m examples.indigo_insertion_synthetic \
    --absolute-position-embeddings gpt2 --relative-position-embeddings indigo \
    --pointer-context gap_window --min-delete-length 1 --max-delete-length 4 \
    --curriculum-epochs $curr --max-epochs 500 --lr 3e-4 --num-workers 4 \
    --output-dir outputs/indigo_curr${curr}_span_1_4
done

# Tombstones. Defaults are span_counts=(0,1) at 50/50, max_span=1, donor,
# random placement, d192/L8, 500 epochs, T=1.0.
$P -m examples.indigo_tombstone_synthetic \
  --num-workers 4 --output-dir outputs/indigo_tombstone_mix
# The always-1-junk arm this replaced:
$P -m examples.indigo_tombstone_synthetic \
  --junk-span-counts 1 --num-workers 4 --output-dir outputs/indigo_tombstone_1x1

# Pre-flight entropy probe; trains nothing.
PYTHONPATH=$PWD $P examples/tombstone_entropy_probe.py --trials 2000

# Tests.
$P -m pytest tests/ -q --ignore=tests/test_oracle.py --ignore=tests/test_reinforce.py
```

Scripts must be run either as modules from the repository root (`python -m
examples.foo`) or directly, now that shared task code lives in `arid.island_task`
and `arid.island_plots` rather than being imported across example scripts.
`tests/test_oracle.py` and `tests/test_reinforce.py` fail at import because
`arid/oracle.py` and `arid/reinforce.py` do not exist; this predates the work
here.
