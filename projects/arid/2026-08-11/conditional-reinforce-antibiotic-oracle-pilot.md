---
title: Conditional REINFORCE antibiotic-oracle pilot
date: 2026-08-11
project: arid
agent: Codex
status: complete
sources:
  - /scratch/users/diamant/arid_moses_reinforce/pilot_g16_lr1e-5/run_metadata.json
  - /scratch/users/diamant/arid_moses_reinforce/pilot_g16_lr1e-5/metrics.json
  - /scratch/users/diamant/arid_moses_reinforce/pilot_g16_lr1e-5/finetuned_trajectory_scores.json
  - /scratch/users/diamant/arid_moses_reinforce/pilot_g16_lr1e-5/finetuned_trajectories.pdf
  - user assessment of the completed run and trajectory visualizations
tags:
  - smiles
  - reinforcement-learning
  - reinforce
  - conditional-generation
  - antibiotics
  - oracle
---

# Summary

Conditional REINFORCE successfully fine-tuned the pretrained ARID MOSES edit
policy against the portable Morgan-fingerprint antibiotic oracle. The selected
run used 128 molecule-conditional rollouts per optimizer step, forced the first
edit, began at reverse step 5, and allowed up to eight edits. Its best fixed-set
validation mean oracle log-probability improvement was **0.4383** at epoch 40.
At the final epoch, mean validation improvement remained **0.4073**, 67.76% of
validation trajectories improved, and trajectories averaged 1.95 edits.

The user judged the aggregate metrics and edit visualizations promising. The
result establishes that grouped leave-one-out REINFORCE can shift the policy
toward higher oracle scores while retaining approximately 99.4% final validity
and making compact edits. Evaluation on fresh molecules and fresh rollout seeds
is still needed to measure generalization and oracle exploitation.

# Key Points

- Best `val/mean_reward_improvement`: `0.4383128881` at zero-indexed epoch 40.
- Final `val/mean_reward_improvement`: `0.4072883725`, corresponding to a
  geometric mean final-to-starting oracle-probability ratio of about
  `exp(0.4073) = 1.50`.
- Final `val/positive_reward_fraction`: `0.6776123047`; its maximum over training
  was `0.6906738281`.
- Final `val/mean_num_edits`: `1.9494628906`.
- Final validation validity was 99.44%, changed fraction was 97.13%, and natural
  termination before the edit cap was 99.54%.
- The final validation mean starting and ending oracle log probabilities were
  `-5.6241` and `-5.2168`, respectively.
- Training used one starting molecule per optimizer step. Each molecule produced
  a group of 128 trajectories for the leave-one-out baseline.
- The output path contains `g16` in its name, but `run_metadata.json` confirms
  that the actual group size was 128.

# Details

## Command and configuration

The run was launched from `/home/users/diamant/repos/ARID` with:

```bash
python examples/moses_reinforce.py \
  --policy-run-dir /scratch/users/diamant/arid_moses_training/h100_full_d384_l3_b512_w4 \
  --oracle-checkpoint /scratch/users/diamant/morgan_FP_abx_classifiers_v2/best.ckpt \
  --output-dir /scratch/users/diamant/arid_moses_reinforce/pilot_g16_lr1e-5 \
  --train-samples 256 \
  --val-samples 64 \
  --group-size 128 \
  --start-reverse-step 5 \
  --max-edits 8 \
  --temperature 1.0 \
  --lr 1e-5 \
  --max-epochs 50 \
  --plot-samples 12 \
  --wandb-project arid_moses_reinforce \
  --device cuda
```

The reward was the final oracle log probability minus the starting oracle log
probability. Invalid final SMILES received `-12`; valid oracle log probabilities
were floored at `-10`. The first action was forced to be an edit, after which the
policy could terminate normally. The reverse-step counter covered indices 5–12
when all eight edit opportunities were used.

Each optimizer step used one starting molecule and 128 conditional trajectories.
With 256 starting molecules and 50 epochs, the run made 12,800 optimizer updates
and sampled approximately 1.64 million training trajectories. Validation reused
64 fixed MOSES test molecules and fixed rollout seeds, making epoch-to-epoch
comparisons lower variance but not independent.

## Validation behavior

| Metric | Final epoch | Best or maximum |
| --- | ---: | ---: |
| Mean reward improvement | 0.4073 | 0.4383 |
| Positive reward fraction | 67.76% | 69.07% |
| Mean edits | 1.949 | 2.028 |
| Validity | 99.44% | 99.74% |
| Changed fraction | 97.13% | 98.05% |
| Natural termination | 99.54% | 99.77% |

The simultaneous positive mean reward and 67.8% positive-reward fraction show
that the gain was not produced only by a small number of extreme trajectories.
An average of roughly two edits is also consistent with compact conditional
optimization rather than systematically exhausting the eight-edit allowance.

There is a notable train-validation gap: final training mean reward improvement
was `1.8738`, versus `0.4073` on validation. This may reflect adaptation to the
256 training molecules and should be treated as evidence to prioritize a fresh
held-out evaluation rather than simply extending this run.

## Saved outputs

- Best policy: `/scratch/users/diamant/arid_moses_reinforce/pilot_g16_lr1e-5/model.pt`.
- Configuration and best epoch: `run_metadata.json` in the same directory.
- Full logged histories: `metrics.json`.
- Twelve trajectories with start/end oracle log-probability titles:
  `finetuned_trajectories.pdf`.
- Machine-readable scores for those examples: `finetuned_trajectory_scores.json`.

The plotted examples include both improvements and regressions, as expected for
individual stochastic samples. They are qualitative diagnostics rather than a
replacement for the aggregate 8,192-trajectory validation evaluation per epoch.

# Related Notes

- [Portable Morgan fingerprint antibiotic-activity oracle](morgan-antibiotic-oracle.md):
  Documents the selected reward model, its validation enrichment, and its
  applicability-domain limitations.
- [Conditional edit validity and reverse-timestep behavior](conditional-edit-validity-and-timestep-behavior.md):
  Motivates forcing the first edit, incrementing the timestep, and expecting
  naturally short conditional trajectories.
- [MOSES H100 full training reaches perfect validation sample validity](../2026-08-05/moses-h100-full-training-validity.md):
  Records the pretrained edit-policy checkpoint used to initialize this run.

# Open Questions

- Do the gains persist on fresh MOSES starting molecules and fresh rollout seeds?
- How much of the train-validation gap comes from memorizing favorable edits for
  the 256 training molecules?
- Are the high-scoring edits still favorable under independent antibiotic
  predictors or chemistry-based applicability-domain checks?
- Does the policy preserve molecular similarity and other desirable properties
  while increasing the Morgan-oracle score?
- Would checkpoint selection remain stable with a larger independently sampled
  validation set?

# Sources

- `/scratch/users/diamant/arid_moses_reinforce/pilot_g16_lr1e-5/run_metadata.json`
  for the exact configuration, reward definition, best epoch, and best metric.
- `/scratch/users/diamant/arid_moses_reinforce/pilot_g16_lr1e-5/metrics.json` for
  the complete training and validation metric histories.
- `/scratch/users/diamant/arid_moses_reinforce/pilot_g16_lr1e-5/finetuned_trajectories.pdf`
  and `finetuned_trajectory_scores.json` for the final qualitative examples.
- User report that the run worked well and that its aggregate metrics and
  trajectory visualizations appeared promising.
