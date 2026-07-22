---
title: Larger-span MOSES trajectory sweep
date: 2026-07-22
project: arid
agent: codex
status: draft
sources:
  - /home/users/diamant/repos/ARID/examples/moses_trajectory_diagnostics.py
  - /home/users/diamant/repos/ARID/arid/corruption.py
  - /scratch/users/diamant/arid_moses_step4/trajectory_sweep_extended/trajectory_report.md
  - /scratch/users/diamant/arid_moses_step4/trajectory_sweep_extended/trajectory_summary.csv
  - /scratch/users/diamant/arid_moses_step4/trajectory_sweep_extended/forward_trajectories.pdf
tags:
  - smiles
  - moses
  - corruption
  - larger-spans
---

# Summary

An extended MOSES noising sweep tested larger forward deletion spans, which correspond to longer reverse insertions, and larger forward insertion spans, which correspond to longer reverse deletions. After making empty content absorbing in the forward process, the best clean first larger-span schedule is `del12_ins2_tail40_cap80`: it allows 12-token reverse insertions, terminates on every sampled trajectory, keeps p100 max corrupted content length at 50, and has median/p95 trajectory lengths of about 10/17 steps.

# Key Points

- The sweep used 2,048 training molecules, seed 0, and the same three molecules for every schedule in `forward_trajectories.pdf`.
- Empty content is absorbing in the forward process; earlier non-absorbing metrics from this note were superseded by this rerun.
- All nine extended schedules reached empty content for 100% of sampled trajectories.
- No selected forward insertion was rejected by the configured corrupted-length cap.
- Increasing reverse insertion capacity from 8 to 12 tokens is cheap when forward insertions stay short.
- Larger forward insertions produce more reverse-delete targets, but drive up corrupted length and outer-step count.
- `del12_ins4_tail48_cap96` is a good balanced candidate if the deletion head needs more reverse-delete examples.
- 12-token forward insertions look too aggressive for an early run.

# Details

| schedule | p50 steps | p95 steps | p95 max corrupted len | p100 max corrupted len | reverse insert targets | reverse delete targets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| deletion_only_del8_cap64 | 9.0 | 12.0 | 42.0 | 47.0 | 0.899 | 0.000 |
| del8_ins2_tail32_cap80 | 13.0 | 21.0 | 43.0 | 51.0 | 0.694 | 0.237 |
| del12_ins2_tail40_cap80 | 10.0 | 17.0 | 43.0 | 50.0 | 0.666 | 0.246 |
| del8_ins4_tail40_cap80 | 18.0 | 30.0 | 47.0 | 66.0 | 0.633 | 0.316 |
| del12_ins4_tail48_cap96 | 13.0 | 24.0 | 47.0 | 62.0 | 0.601 | 0.330 |
| del8_ins8_tail64_cap128 | 27.0 | 47.0 | 64.0 | 108.0 | 0.638 | 0.327 |
| del12_ins8_tail72_cap128 | 18.0 | 36.0 | 60.0 | 97.0 | 0.590 | 0.361 |
| del16_ins8_tail96_cap160 | 14.0 | 28.0 | 57.0 | 110.0 | 0.573 | 0.364 |
| del12_ins12_tail96_cap160 | 31.0 | 62.0 | 93.0 | 160.0 | 0.587 | 0.384 |

The practical tradeoff is asymmetric. Larger reverse insertions, caused by longer forward deletions, reduce outer iterations without increasing corrupted length much. Larger reverse deletions, caused by longer forward insertions, are architecturally cheaper per target but they create longer corrupted states and longer trajectories. Absorbing empty removes post-empty noise excursions, shortening trajectories and lowering reverse-delete target fractions. This supports trying group-sized reverse insertions before training with very large forward insertions.

# Related Notes

- [MOSES Step 4 trajectory sweep](moses-step4-trajectory-sweep.md): Establishes the baseline 4-token-span sweep used for comparison.

# Open Questions

- Whether `del12_ins2_tail40_cap80` improves molecule validity or reconstruction over the prior `balanced_del4_ins2_cap64` setting.
- Whether 12-token random SMILES spans actually correspond to meaningful chemical groups often enough, or whether later span sampling should become branch/subgraph aware.
- Whether the reverse insertion decoder needs more capacity for 12-token spans.

# Sources

- ARID diagnostics and corruption implementation listed in the frontmatter.
- Generated artifacts under `/scratch/users/diamant/arid_moses_step4/trajectory_sweep_extended/`.
