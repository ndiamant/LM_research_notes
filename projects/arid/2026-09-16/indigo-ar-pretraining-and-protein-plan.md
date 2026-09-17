---
title: INDIGO autoregressive pretraining status and Swiss-Prot plan
date: 2026-09-16
project: arid
agent: codex
status: draft
sources:
  - /home/users/diamant/repos/ARID/examples/indigo_pretrain_finetune_synthetic.py
  - /home/users/diamant/repos/ARID/tests/test_indigo_pretrain_finetune.py
  - /home/users/diamant/repos/ARID/outputs/indigo_pretrain_finetune_synthetic
  - /home/users/diamant/repos/ARID/outputs/indigo_pretrain_finetune_0_100_168
  - /home/users/diamant/repos/ARID/outputs/indigo_pretrain_finetune_50_100_168
  - /home/users/diamant/repos/ARID/outputs/indigo_pretrain_finetune_0_100_161
  - /home/users/diamant/repos/ARID/outputs/indigo_pretrain_finetune_50_100_161
  - /scratch/users/diamant/arid_runs/partition_3_random_25_epochs/config.yaml
tags:
  - indigo
  - autoregressive
  - pretraining
  - synthetic
  - swissprot
---

# Summary

The synthetic experiments provide credible early evidence that left-to-right
autoregressive pretraining transfers to INDIGO's token-value prediction, while
the position/pointer head can initially learn more slowly than from random
initialization. Transfer is strongest when the insertion schedule stays close
to left-to-right autoregression: it is nearly perfect for fixed deletion length
16 and remains strong early for lengths 8--16. For the harder 1--16 schedule,
pretraining substantially improves the value loss but initially damages the
pointer loss; both models largely converge by 100 epochs.

The present evidence is sufficient to move on rather than optimize the
synthetic pointer architecture prematurely. The next experiment should compare
randomly initialized and AR-pretrained INDIGO models on Swiss-Prot, with the
same random three-partition schedule and approximately the same transformer
capacity as the existing `partition_3_random_25_epochs` ReverseEdit run.
ProGen2 should follow only after an architecture-matched in-house GPT-2
experiment establishes whether protein-scale AR transfer helps.

# Key Points

- `examples/indigo_pretrain_finetune_synthetic.py` pretrains a vanilla GPT-2
  next-token model and transfers its token embeddings, absolute position
  embeddings, transformer blocks, final normalization, and tied language-model
  head into a script-local INDIGO subclass.
- New INDIGO relation biases are initialized to zero. The value adapter begins
  as the pretrained AR function plus zero-initialized corrections from the
  chosen gap's left and right context. The pointer head is newly initialized.
- Absolute positional embeddings should remain transferred. Reinitializing
  WPE significantly hurt early performance, consistent with WPE and backbone
  weights being co-adapted.
- Adding learned absolute gap-position embeddings to the pointer made little
  difference and slightly worsened the final result. There is no current
  evidence for retaining that added complexity.
- The most plausible unresolved explanation for slower pretrained pointer
  learning is optimization conflict: the pretrained backbone initially
  preserves a useful value predictor while the new pointer objective asks it
  to reorganize. The earlier embedding-geometry explanation remains
  speculative.
- A pointer-specific transformer layer could preserve one-shot log
  probabilities and cached generation if it is causal and has its own KV
  cache, but it would add a capacity confound before the basic protein transfer
  question has been tested.
- Existing insertion and tombstone checkpoints remain a separate concern. The
  transfer architecture is intentionally script-local rather than replacing
  the base INDIGO module used by the older scripts.
- Current code hygiene issue: the experiment tests still exercise
  `pointer_gap_position_embeddings` and `reinitialize_wpe`, but those arguments
  are absent from the current script. On 2026-09-16 the focused test file had
  3 passing and 2 failing tests, both due to these missing arguments. The main
  pretrain/finetune path remains present; either restore the ablation options or
  remove their stale tests before treating the branch as clean.

# Details

## Synthetic architecture

The vanilla and INDIGO models use matched GPT-2 backbones. The INDIGO model
retains causal attention and adds zero-initialized INDIGO relative-attention
biases. The transferred language-model head continues to predict token values.
For each trajectory step, the value adapter combines the current causal hidden
state with selected left- and right-gap states:

```text
adapted value state = current hidden
                    + W_left(left-gap hidden)
                    + W_right(right-gap hidden)
```

Both correction matrices start at zero, so transfer initially reproduces the
AR value function. The width-one `gap_window` pointer represents a candidate
gap by its left and right endpoints. Its candidate dimension is a scoring
dimension rather than an additional token in the transformer sequence, which
allows all insertion positions to be scored in one pass and preserves cached
generation.

Setting `--pretrain-epochs 0` skips AR optimization and supplies the same
architecture with random weights, providing a true random-initialization
control. This is important because one epoch was already a substantial
pretraining treatment: the earlier one-epoch run reached approximately 79.4%
AR validation token accuracy.

## Synthetic results

### Strongly any-order schedule: deletion lengths 1--4

The initial comparison used 100 versus 1 AR pretraining epochs and 500 INDIGO
finetuning epochs:

```bash
python examples/indigo_pretrain_finetune_synthetic.py \
  --pretrain-epochs 100 --finetune-epochs 500

python examples/indigo_pretrain_finetune_synthetic.py \
  --pretrain-epochs 1 --finetune-epochs 500 \
  --output-dir outputs/indigo_pretrain_finetune_synthetic_no_pretrain
```

At finetuning epoch 100, the 100-epoch pretrained model had validity 0.207,
validation loss 0.697, and off-manifold mass 0.0552, versus 0.512, 0.663, and
0.0348 for the one-epoch model. The latter reached 50% validity around epoch
100, compared with epoch 140 for strong pretraining. By epoch 500 the gap was
negligible: validity was 0.818 versus 0.826 and loss 0.635 versus 0.633.

This result does not imply that AR pretraining is generally harmful. It shows
that a strong next-token solution can slow adaptation to a distant arbitrary-
gap factorization, while providing little asymptotic advantage on this very
small task.

### Near-autoregressive schedules

Fixed deletion length 16 transferred nearly perfectly. This is an important
positive control: it validates the weight-transfer implementation and shows
that pretraining is useful when INDIGO's generation order matches the AR
factorization.

For deletion lengths 8--16, the following comparison showed strong early value
transfer:

```bash
python examples/indigo_pretrain_finetune_synthetic.py \
  --pretrain-epochs 0 --finetune-epochs 100 \
  --output-dir outputs/indigo_pretrain_finetune_0_100_168 \
  --min-delete-length 8 --max-delete-length 16

python examples/indigo_pretrain_finetune_synthetic.py \
  --pretrain-epochs 50 --finetune-epochs 100 \
  --output-dir outputs/indigo_pretrain_finetune_50_100_168 \
  --min-delete-length 8 --max-delete-length 16
```

| Epoch | Scratch value loss | Pretrained value loss | Scratch pointer loss | Pretrained pointer loss |
|---:|---:|---:|---:|---:|
| 1 | 1.030 | **0.367** | **0.621** | 1.317 |
| 2 | 0.752 | **0.360** | **0.317** | 0.327 |
| 3 | 0.556 | **0.287** | 0.220 | **0.219** |

The pretrained total loss was better by epoch 2 (0.343 versus 0.528), but the
benefit mostly disappeared by epoch 10. Both runs finished near 98% validity
(0.982 scratch and 0.984 pretrained).

For the broader 1--16 schedule:

```bash
python examples/indigo_pretrain_finetune_synthetic.py \
  --pretrain-epochs 0 --finetune-epochs 100 \
  --output-dir outputs/indigo_pretrain_finetune_0_100_161 \
  --min-delete-length 1 --max-delete-length 16

python examples/indigo_pretrain_finetune_synthetic.py \
  --pretrain-epochs 50 --finetune-epochs 100 \
  --output-dir outputs/indigo_pretrain_finetune_50_100_161 \
  --min-delete-length 1 --max-delete-length 16
```

| Epoch | Initialization | Total loss | Pointer loss | Value loss | Pointer acc. | Value acc. |
|---:|---|---:|---:|---:|---:|---:|
| 1 | Scratch | **0.942** | **0.881** | 1.007 | **0.816** | 0.717 |
| 1 | AR-50 | 1.083 | 1.628 | **0.503** | 0.432 | **0.761** |
| 2 | Scratch | 0.685 | **0.619** | 0.755 | not recorded | not recorded |
| 2 | AR-50 | **0.614** | 0.731 | **0.489** | not recorded | not recorded |

By epoch 100 the models were close: scratch loss 0.3206, validity 0.7422, and
off-manifold mass 0.0159; pretrained loss 0.3210, validity 0.7617, and
off-manifold mass 0.0161. The clean conclusion is therefore early positive
transfer to value prediction accompanied by an early pointer penalty.

## Ablations and interpretation

Reinitializing the absolute positional embeddings did not cure the pointer
penalty. Relative to full transfer, fresh WPE worsened epoch-1 total loss from
1.083 to 1.147, pointer loss from 1.628 to 1.693, and value loss from 0.503 to
0.567. Epoch-10 validity fell from 0.270 to 0.182, although by epoch 100 it was
again similar (0.762 versus 0.754). The default should therefore retain WPE.

Learned gap-position embeddings added a head-only term to each pointer score,
so they preserved both one-shot trajectory log probabilities and cached
generation. They had little effect: epoch-1 pointer loss was 1.645 versus 1.628
without them; modest middle-epoch improvements did not persist, and final
pointer loss and validity were slightly worse (0.384 and 0.746 versus 0.367
and 0.762). The comparison may also have a small RNG-stream confound because
constructing the optional embedding consumes random draws.

The training schedules themselves make pointer prediction increasingly close
to a boundary shortcut. The target insertion was the rightmost gap in about
26.45%, 39.87%, 51.85%, and 100% of insertion examples for schedules 1--4,
1--16, 8--16, and 16--16 respectively. Boundary rates were approximately
32.8%, 43.62%, 53.17%, and 100%. This explains why schedule distance matters,
but the unsuccessful explicit gap-position feature suggests it is not the
whole story.

If the pointer behavior becomes important on proteins, the cheapest diagnostic
sequence is:

1. Record pointer and value losses before the first finetuning update. A worse
   pretrained pointer at step zero suggests a readout/accessibility problem;
   divergence only after updates suggests gradient competition.
2. Measure cosine similarity between pointer and value gradients on the shared
   backbone.
3. Try a short pointer-only warmup or a pointer-specific residual MLP before a
   new transformer layer.
4. Add a causal pointer-adaptation layer only if these diagnostics show that
   additional contextual processing is needed.

## Swiss-Prot experiment plan

The primary existing comparison is:

```bash
python examples/swissprot_proteins.py \
  --config /scratch/users/diamant/arid_runs/partition_3_random_25_epochs/config.yaml
```

That run uses the full Swiss-Prot split, batch size 64, 25 epochs, learning
rate `3e-4`, weight decay `1e-2`, and seed 7. Its ReverseEdit model has
`d_model=384`, six encoder layers, four decoder layers, six attention heads,
and feed-forward width 1536, for approximately 20.4M parameters. The schedule
uses exactly three chunks in random generation order.

Create a separate `examples/swissprot_indigo_pretrain_finetune.py` rather than
changing the existing Swiss-Prot or base INDIGO paths. Use an 11-layer GPT-2
backbone with width 384, six heads, and feed-forward width 1536, approximately
19.6M parameters. The first controlled comparison should contain:

1. The existing random-three-partition ReverseEdit baseline.
2. Protein INDIGO with random initialization.
3. The identical protein INDIGO initialized from an architecture-matched
   left-to-right Swiss-Prot AR model.

Use the same data split, tokenizer, batch size, seed, partition sampler, and 25
INDIGO epochs. Convert each sampled random three-chunk trajectory into tokenwise
insertions: preserve the sampled chunk order, then insert residues within each
chunk left-to-right at the correct gap, followed by DONE. Train on all causal
trajectory decisions in one pass.

This comparison is schedule- and capacity-matched, but not generation-step-
matched. ReverseEdit generates three spans; tokenwise INDIGO makes roughly one
edit per residue. Report optimizer updates, supervised residue decisions, wall
clock, and sampling cost alongside epoch counts. Track pointer loss/accuracy,
residue value loss/accuracy, termination, length, composition, uniqueness, and
the existing ESM-C PLL/PPPL metrics.

The current `examples/swissprot_autoregressive.py` uses
`nn.TransformerEncoder`, so its checkpoint is not directly transferable to the
GPT-2 INDIGO backbone. The protein experiment needs an architecture-matched
GPT-2 AR pretrainer. It will also need a protein-specific value vocabulary
rather than the base synthetic model's hard-coded zero/one content-token
assumptions; initially keeping this in the new experiment module minimizes
checkpoint risk.

Once random-init versus in-house AR transfer is understood, repeat the
pretrained condition with ProGen2. ProGen2 should not be the first comparison
because it simultaneously changes pretraining scale, tokenizer, positional
conventions, and potentially backbone details.

# Related Notes

- [Swiss-Prot temperature-1 results and append-3 supervision scaling](../2026-09-15/swissprot-temperature-1-and-supervision-scaling.md): Establishes the temperature-1 evaluation protocol, the importance of dense residue supervision, and the motivation for tokenwise INDIGO.
- [Swiss-Prot partition-schedule results](../2026-09-10/swissprot-partition-schedule-results.md): Defines the random three-partition baseline and reports its position-prediction and ESM-C results.
- [Swiss-Prot setup and matched baseline runs](../2026-09-06/swissprot-setup-and-baseline-runs.md): Records the dataset, environment, capacity-matched AR baseline, and evaluation assets.

# Open Questions

- Does AR pretraining improve protein INDIGO value prediction enough to improve
  ESM-C PPPL, rather than merely accelerate the first few epochs?
- Does the pretrained pointer penalty persist on diverse protein sequences, or
  is it peculiar to the tiny binary-island task?
- Is any pointer penalty already present before the first finetuning update, or
  is it induced by competition between value and pointer gradients?
- How should tokenwise INDIGO's much larger number of sampling steps be handled
  in a fair quality/compute comparison with three-span ReverseEdit generation?
- Should the stale synthetic ablation tests be removed, or should their
  optional script-local implementations be restored for reproducibility?
- After an in-house GPT-2 comparison, how much additional transfer comes from
  ProGen2 and how much adaptation is required for its tokenizer and positional
  scheme?

# Sources

- Synthetic experiment implementation and focused tests listed in the
  frontmatter; test status was checked with the `arid_indigo` conda environment
  on 2026-09-16.
- Metrics and interpretations from the synthetic runs and commands documented
  above, as reported and reviewed during the 2026-09 ARID experiment session.
- Swiss-Prot random-partition configuration at
  `/scratch/users/diamant/arid_runs/partition_3_random_25_epochs/config.yaml`.
- The three related ARID research notes linked above.
