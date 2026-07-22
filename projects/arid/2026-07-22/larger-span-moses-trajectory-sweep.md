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

An extended MOSES noising sweep tested larger forward deletion spans, which correspond to longer reverse insertions, and larger forward insertion spans, which correspond to longer reverse deletions. The best first larger-span schedule is `del12_ins2_tail40_cap80`: it allows 12-token reverse insertions, terminates on every sampled trajectory, keeps p100 max corrupted content length at 51, and has median/p95 trajectory lengths of about 23/31 steps.

# Key Points

- The sweep used 2,048 training molecules, seed 0, and the same three molecules for every schedule in `forward_trajectories.pdf`.
- All nine extended schedules reached empty content for 100% of sampled trajectories.
- No selected forward insertion was rejected by the configured corrupted-length cap.
- Increasing reverse insertion capacity from 8 to 12 tokens is cheap when forward insertions stay short.
- Larger forward insertions produce more reverse-delete targets, but drive up corrupted length and outer-step count.
- 12-token forward insertions look too aggressive for an early run.

# Details

| schedule | p50 steps | p95 steps | p95 max corrupted len | p100 max corrupted len | reverse insert targets | reverse delete targets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| deletion_only_del8_cap64 | 9.0 | 12.0 | 42.0 | 47.0 | 0.899 | 0.000 |
| del8_ins2_tail32_cap80 | 22.0 | 29.0 | 44.0 | 53.0 | 0.633 | 0.323 |
| del12_ins2_tail40_cap80 | 23.0 | 30.7 | 43.0 | 51.0 | 0.592 | 0.365 |
| del8_ins4_tail40_cap80 | 31.0 | 40.0 | 47.0 | 62.0 | 0.609 | 0.360 |
| del12_ins4_tail48_cap96 | 32.0 | 42.0 | 46.0 | 59.0 | 0.580 | 0.390 |
| del8_ins8_tail64_cap128 | 42.0 | 53.0 | 65.7 | 109.0 | 0.633 | 0.343 |
| del12_ins8_tail72_cap128 | 42.0 | 53.0 | 61.0 | 95.0 | 0.598 | 0.379 |
| del16_ins8_tail96_cap160 | 37.0 | 47.0 | 58.0 | 103.0 | 0.586 | 0.387 |
| del12_ins12_tail96_cap160 | 64.0 | 76.0 | 96.0 | 151.0 | 0.606 | 0.378 |

The practical tradeoff is asymmetric. Larger reverse insertions, caused by longer forward deletions, reduce outer iterations without increasing corrupted length much. Larger reverse deletions, caused by longer forward insertions, are architecturally cheaper per target but they create longer corrupted states and longer trajectories. This supports trying group-sized reverse insertions before training with very large forward insertions.

# Related Notes

- [MOSES Step 4 trajectory sweep](moses-step4-trajectory-sweep.md): Establishes the baseline 4-token-span sweep used for comparison.

# Open Questions

- Whether `del12_ins2_tail40_cap80` improves molecule validity or reconstruction over the prior `balanced_del4_ins2_cap64` setting.
- Whether 12-token random SMILES spans actually correspond to meaningful chemical groups often enough, or whether later span sampling should become branch/subgraph aware.
- Whether the reverse insertion decoder needs more capacity for 12-token spans.

# Sources

- ARID diagnostics and corruption implementation listed in the frontmatter.
- Generated artifacts under `/scratch/users/diamant/arid_moses_step4/trajectory_sweep_extended/`.
