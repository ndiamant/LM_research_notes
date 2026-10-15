---
title: PRISM speedup from shared DNA encoding and skipped wildtype bases (commit 66d6585) on H100
date: 2026-10-14
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - SMFNet commit 66d6585 ("faster PRISM"), src/smf_net/evaluation/ism.py
  - /home/users/diamant/repos/SMFNet/scratch/prism_bench/bench.py
  - /home/users/diamant/repos/SMFNet/scratch/prism_bench/check.py
tags: [prism, ism, performance, benchmark, h100]
---

# Summary

SMFNet commit `66d6585` makes `single_read_ism` (PRISM) about **1.7–1.8x faster** on one H100. Outputs are effectively unchanged: deltas correlate at 0.99997 with the pre-commit code.

# Key Points

- The commit changes two things in `src/smf_net/evaluation/ism.py`:
  1. **Loop order swap / shared DNA encoding.** The old code looped over unmasking orders and, for each order, over mutant batches, so the DNA encoder ran again for every order. The new code loops over mutant batches, encodes each batch once (`_dna_embedding`), and scores that embedding on every trajectory (`_score_mutations_across_trajectories`).
  2. **Wildtype base skipped.** The reference base at each position now reuses the reference score instead of being scored again. This cuts the mutants by 25%.
- Skipping the wildtype base alone accounts for about 1.33x. The remaining ~1.3x is *inferred* to come from the shared encoding; the two changes were not timed separately.
- The speedup is flat across window length and rises slightly with more orders (more trajectories share each encoding).

# Details

Setup: Sherlock node `sh04-03n07`, job 47845095, one idle NVIDIA H100 80GB HBM3. Model: `/scratch/users/diamant/smf_models/SMFNet_full_run_2026-08-11/epoch=377-step=171990.ckpt`, fp32, compiled. Input 2048 bp, output 1024, `sample_steps=32`. Inputs were random DNA and a random binary read with ~60% missing calls. Batch size was 128 and `sampling_mode="cosine"`. Each time is the best of 2 runs after a warmup.

| Settings | Old | New | Speedup |
|---|---|---|---|
| 100 bp window, 8 orders | 5.20 s | 3.04 s | 1.71x |
| 100 bp window, 32 orders | 20.86 s | 11.67 s | 1.79x |
| 300 bp window, 8 orders | 15.65 s | 9.22 s | 1.70x |
| 100 bp window, 8 orders, `context_mode="masked"`, 100 bp likelihood window | 5.23 s | 3.07 s | 1.70x |

Numerical agreement (100 bp, 8 orders): the largest |old − new| delta is 1.1e-2, against a median |delta| of 0.30 and a max of 3.4. Pearson r is 0.99997. Repeat runs of the new code are bit-identical. The small differences are *probably* floating-point effects of different batch shapes, plus wildtype entries now being exactly zero. As a result, PRISM outputs saved before this commit won't reproduce bit-for-bit (for example in `scripts/prism_order_example.py`).

Commands (run from the SMFNet repo root, directly on the compute node):

```bash
git show 66d6585^:src/smf_net/evaluation/ism.py > scratch/prism_bench/ism_old.py
$GROUP_HOME/uv-envs/SMFNet-cu124/bin/python scratch/prism_bench/bench.py   # timing table
$GROUP_HOME/uv-envs/SMFNet-cu124/bin/python scratch/prism_bench/check.py   # numerical agreement
```

The `scratch/prism_bench/` scripts are untracked in the SMFNet repo.

# Related Notes

- None found. Earlier SMFNet notes in this project don't cover PRISM/ISM runtime.

# Open Questions

- How much speedup does each change give on its own (shared encoding vs. skipped wildtype)?
- Does the speedup hold in bf16, at larger batch sizes, or with `sampling_mode="one_pos_at_a_time"`, which has many more trajectory steps?

# Sources

- SMFNet commit `66d6585`, `src/smf_net/evaluation/ism.py` (`_score_mutations_across_trajectories`, `_score_read_across_orders`)
- Benchmark scripts: `/home/users/diamant/repos/SMFNet/scratch/prism_bench/{bench.py,check.py,ism_old.py}`
