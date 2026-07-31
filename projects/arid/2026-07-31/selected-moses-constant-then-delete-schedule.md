---
title: Selected MOSES constant-then-delete schedule
date: 2026-07-31
project: arid
agent: codex
status: draft
sources:
  - /home/users/diamant/repos/ARID/arid/schedules.py
  - /home/users/diamant/repos/ARID/examples/moses_trajectory_diagnostics.py
  - /scratch/users/diamant/arid_moses_step4/trajectory_sweep_constant_then_delete/trajectory_report.md
  - /scratch/users/diamant/arid_moses_step4/trajectory_sweep_constant_then_delete/trajectory_summary.csv
  - /scratch/users/diamant/arid_moses_step4/trajectory_sweep_constant_then_delete/forward_trajectories.pdf
tags:
  - smiles
  - moses
  - schedule-decision
  - corruption
---

# Summary

Use `del12_ins4_mix10_p50_cap96` as the next MOSES forward noising schedule. This is the simple constant-then-delete schedule: for the first 10 forward steps, choose insertion with probability 0.5 and deletion with probability 0.5; after step 10, switch to delete-only until empty. Empty content is absorbing.

# Key Points

- Selected schedule: `del12_ins4_mix10_p50_cap96`.
- Forward deletion max span: 12 tokens, so reverse insertion can train on up to 12-token spans.
- Forward insertion max span: 4 tokens, so reverse deletion can delete up to 4-token noisy spans.
- Mixed window: 10 forward steps with constant 50% insertion probability.
- Corrupted content cap: 96 tokens, excluding BOS/EOS.
- Measured on 2,048 sampled MOSES training molecules: 100% reached empty, median/p95 forward steps 14/18, p95/p100 max corrupted content length 50/70, and reverse-delete targets for 33.1% of states.

# Details

This schedule was chosen over `del12_ins2_mix10_p50_cap80` because `ins4` better matches the desired reverse deletion head capacity while remaining short and bounded. The `ins2` variant is slightly cleaner in corrupted length, with p95/p100 max content length 46/59 and reverse-delete targets for 35.2% of states, but it only trains 1-2 token reverse deletions.

The 10-step mixed window was preferred over longer windows because it gives enough reverse-delete signal without extending trajectories. In the same sweep, `del12_ins4_mix12_p50_cap96` increased reverse-delete targets to 36.2% with p50/p95 steps 15/20, and `del12_ins4_mix16_p50_cap96` increased reverse-delete targets to 40.2% with p50/p95 steps 18/23. The 24-step mixed window had p95 steps 30 and looked too long for the first simple run.

Command provenance:

```bash
cd /home/users/diamant/repos/ARID
/home/groups/btrippe/diamant/miniforge/envs/g2pt/bin/python -m py_compile arid/schedules.py examples/moses_trajectory_diagnostics.py
/home/groups/btrippe/diamant/miniforge/envs/g2pt/bin/python examples/moses_trajectory_diagnostics.py --schedule-set constant --sample-size 2048 --plot-molecules 3 --output-dir /scratch/users/diamant/arid_moses_step4/trajectory_sweep_constant_then_delete
```

The AGENTS.md Python path `/Users/ndiamant/miniforge3/envs/simp_diff/bin/python` was not present on the Sherlock host, so the available `g2pt` environment was used for this diagnostic.

# Related Notes

- [Larger-span MOSES trajectory sweep](../2026-07-22/larger-span-moses-trajectory-sweep.md): Establishes the larger-span, absorbing-empty sweep that motivated simpler constant-then-delete schedules.
- [MOSES Step 4 trajectory sweep](../2026-07-22/moses-step4-trajectory-sweep.md): Establishes the original MOSES Step 4 diagnostic pipeline and baseline schedule comparisons.

# Open Questions

- Whether 33.1% reverse-delete targets is enough for stable deletion-head learning once training begins.
- Whether 12-token reverse insertion spans require a larger insertion decoder than the current first-pass model.
- Whether later span sampling should become branch- or subgraph-aware if random 12-token SMILES spans do not align well with chemical groups.

# Sources

- ARID code paths and generated artifacts listed in the frontmatter.
