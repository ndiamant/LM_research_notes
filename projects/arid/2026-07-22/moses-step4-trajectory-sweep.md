---
title: MOSES Step 4 trajectory sweep
date: 2026-07-22
project: arid
agent: codex
status: draft
sources:
  - /home/users/diamant/repos/ARID/SMILES_TRAINING_PLAN.md
  - /home/users/diamant/repos/ARID/arid/smiles.py
  - /home/users/diamant/repos/ARID/arid/corruption.py
  - /home/users/diamant/repos/ARID/examples/moses_trajectory_diagnostics.py
  - /scratch/users/diamant/arid_moses_step4/trajectory_sweep/trajectory_report.md
  - /scratch/users/diamant/arid_moses_step4/trajectory_sweep/trajectory_summary.csv
tags:
  - smiles
  - moses
  - corruption
  - trajectory-diagnostics
---

# Summary

ARID step 4 now has direct MOSES CSV-gzip loading, a copied DeepChem SMILES regex tokenizer, a training-only vocabulary, a minimal SMILES dataset, and a reusable mixed insertion/deletion corruption sampler. A 2,048-molecule trajectory sweep over official MOSES train examples suggests `balanced_del4_ins2_cap64` is the best first mixed noising schedule: it terminates reliably, keeps max corrupted lengths well under the cap, and creates enough reverse-delete targets to exercise the deletion head.

# Key Points

- Official split sizes observed: train 1,584,663; test 176,074; scaffold test 176,225.
- Training-only vocabulary size is 28 total IDs: 4 special tokens plus 24 lexical SMILES tokens.
- Tokenization round-tripped on the scanned training split, and unknown-token coverage was zero on deterministic validation, test, and scaffold-test.
- A max clean token length of 64 removes 0 molecules from train, test, and scaffold-test; observed full-split max token lengths were 54, 51, and 51 respectively.
- All six swept schedules reached empty content for 100% of the sampled trajectories.
- No sampled forward insertion was rejected by the corrupted-length cap in this sweep.

# Details

The sweep used 2,048 retained training molecules sampled with seed 0. Inserted noise tokens were sampled from the training lexical unigram distribution, excluding PAD, BOS, EOS, and UNK.

| schedule | p50 steps | p95 steps | p95 max corrupted len | p100 max corrupted len | reverse insert targets | reverse delete targets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| deletion_only_del4_cap64 | 15.0 | 19.0 | 42.0 | 47.0 | 0.936 | 0.000 |
| low_insert_del4_ins1_cap64 | 18.0 | 24.0 | 42.0 | 48.0 | 0.820 | 0.128 |
| balanced_del4_ins2_cap64 | 25.0 | 33.0 | 43.0 | 50.0 | 0.729 | 0.233 |
| balanced_del4_ins3_cap80 | 31.0 | 40.0 | 46.0 | 61.0 | 0.703 | 0.266 |
| high_insert_del4_ins4_cap96 | 48.0 | 60.0 | 56.0 | 71.0 | 0.650 | 0.330 |
| slow_delete_high_insert_cap96 | 83.0 | 101.7 | 71.0 | 96.0 | 0.724 | 0.264 |

The recommended first mixed setting is `balanced_del4_ins2_cap64`. It gives a meaningful reverse-delete signal while staying close enough to the deletion-only baseline to keep training and debugging practical. `high_insert_del4_ins4_cap96` is useful as a later stress test. `slow_delete_high_insert_cap96` is probably too expensive for a first run because p95 trajectories exceed 100 forward steps and the maximum sampled corrupted length reaches the cap.

# Related Notes

- None.

# Open Questions

- Whether the 23% reverse-delete target rate from `balanced_del4_ins2_cap64` is sufficient for stable deletion-head learning.
- Whether timestep sampling should be adjusted after observing actual mode-loss and deletion-loss gradients.
- Whether the corrupted-length cap should remain at content length 64 or move to 80 for the first bounded mixed training run.

# Sources

- ARID plan and implementation files listed in the frontmatter.
- Generated diagnostic artifacts under `/scratch/users/diamant/arid_moses_step4/trajectory_sweep/`.
