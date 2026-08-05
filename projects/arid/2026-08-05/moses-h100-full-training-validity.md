---
title: MOSES H100 full training reaches perfect validation sample validity
date: 2026-08-05
project: arid
agent: codex
status: draft
sources:
  - /home/users/diamant/repos/ARID/examples/moses_smiles.py
  - /home/users/diamant/repos/ARID/examples/moses_sample_trajectories.py
  - /scratch/users/diamant/arid_moses_training/h100_full_d192_l3_b512_w8/run_metadata.json
  - /scratch/users/diamant/arid_moses_training/h100_full_d192_l3_b512_w8/metrics.json
  - /scratch/users/diamant/arid_moses_sampling/h100_full_d192_l3_b512_w8/sample_trajectories.pdf
  - /scratch/users/diamant/arid_moses_sampling/h100_full_d192_l3_b512_w8/samples.csv
  - /scratch/users/diamant/arid_moses_sampling/h100_full_d192_l3_b512_w8/sample_trajectories.jsonl
tags:
  - smiles
  - moses
  - training-result
  - h100
  - sampling
---

# Summary

A full-dataset MOSES ARID training run using the selected `del12_ins4_mix10_p50_cap96` schedule trained successfully on one H100. The final validation sampling callback reached `val/sample_validity = 1.0`, meaning 32/32 sampled SMILES were RDKit-valid on the last sampled validation epoch.

# Key Points

- Successful run directory: `/scratch/users/diamant/arid_moses_training/h100_full_d192_l3_b512_w8`.
- W&B project: `arid-smiles`; local W&B run file is under `/scratch/users/diamant/arid_moses_training/h100_full_d192_l3_b512_w8/wandb/run-20260804_165326-fg6zi5v9/`.
- Model size: `d_model=192`, `nhead=8`, 3 encoder layers, 3 decoder layers, feedforward width 768.
- Optimization: `lr=1e-3`, `dropout=0`, `wd=0`, batch size 512, 8 workers, 100 epochs.
- Data: 1,500,000 MOSES train examples and 79,000 validation examples after the script's `max_clean_length=64` filter.
- Final metrics parsed from `metrics.json`: `val/sample_validity=1.0`, `val/loss=5.9104`, `train/loss_epoch=5.9260`, `val/mode_acc=0.7859`, `val/span_token_acc=0.5047`.
- The final five sampled-validity values were `[0.90625, 0.90625, 0.90625, 0.96875, 1.0]`.

# Details

The run used the constant-then-delete schedule selected in the prior note: 10 mixed forward steps with 50% insertion probability, then delete-only until empty. Forward deletion max length was 12 tokens, forward insertion max length was 4 tokens, corrupted content was capped at 96 tokens, and reverse trajectories were capped at 96 steps.

Training command:

```bash
cd /home/users/diamant/repos/ARID
python examples/moses_smiles.py \
  --accelerator gpu \
  --devices 1 \
  --wandb-project arid-smiles \
  --train-samples 1500000 \
  --val-samples 79000 \
  --batch-size 512 \
  --num-workers 8 \
  --max-epochs 100 \
  --dropout 0 \
  --wd 0 \
  --d-model 192 \
  --nhead 8 \
  --num-encoder-layers 3 \
  --num-decoder-layers 3 \
  --dim-feedforward 768 \
  --output-dir /scratch/users/diamant/arid_moses_training/h100_full_d192_l3_b512_w8 \
  --vocab-path /scratch/users/diamant/arid_moses_step4/trajectory_sweep_constant_then_delete/vocab.json \
  --lr 1e-3
```

The durable environment used elsewhere in this ARID work is:

```bash
/home/groups/btrippe/diamant/miniforge/envs/g2pt/bin/python
```

Sample reverse trajectories were visualized after training with `examples/moses_sample_trajectories.py`. The visualization starts from empty `[BOS, EOS]`, then plots each model reverse edit as a row. Green boxes mark model insertions, red boxes mark model deletions, and each page title reports final SMILES, RDKit validity, number of reverse edits, and final token length.

Trajectory visualization command:

```bash
cd /home/users/diamant/repos/ARID
/home/groups/btrippe/diamant/miniforge/envs/g2pt/bin/python examples/moses_sample_trajectories.py \
  --run-dir /scratch/users/diamant/arid_moses_training/h100_full_d192_l3_b512_w8 \
  --output-dir /scratch/users/diamant/arid_moses_sampling/h100_full_d192_l3_b512_w8 \
  --num-samples 24 \
  --temperature 0.8 \
  --device cuda \
  --seed 0
```

That trajectory visualization sample had 22/24 RDKit-valid final SMILES at temperature 0.8. The corresponding artifacts are:

- `/scratch/users/diamant/arid_moses_sampling/h100_full_d192_l3_b512_w8/sample_trajectories.pdf`
- `/scratch/users/diamant/arid_moses_sampling/h100_full_d192_l3_b512_w8/samples.csv`
- `/scratch/users/diamant/arid_moses_sampling/h100_full_d192_l3_b512_w8/sample_trajectories.jsonl`

The sampled trajectories showed both pure insertion paths and paths that used reverse deletions. In the 24 sampled trajectories, several valid samples used 3-6 deletion edits, which is useful evidence that the deletion head is participating in generation rather than being unused. This is observational, not yet a quantitative evaluation of deletion-head necessity.

# Related Notes

- [Selected MOSES constant-then-delete schedule](../2026-07-31/selected-moses-constant-then-delete-schedule.md): Records the exact noising schedule used for this successful training run.
- [Larger-span MOSES trajectory sweep](../2026-07-22/larger-span-moses-trajectory-sweep.md): Motivates 12-token reverse insertion capacity and 4-token reverse deletion capacity.
- [MOSES Step 4 trajectory sweep](../2026-07-22/moses-step4-trajectory-sweep.md): Establishes the original MOSES trajectory diagnostic pipeline.

# Open Questions

- Whether 32/32 sampled validity at validation time remains stable under larger sample sizes, different seeds, lower or higher temperatures, and scaffold-split evaluation.
- Whether reverse deletions materially improve distributional quality or primarily repair local syntax during sampling.
- Whether sampling should expose probability traces for mode, position, and span decisions so trajectory plots can distinguish confident edits from uncertain repairs.

# Sources

- User-provided successful training command and final validation outcome from the 2026-08-05 ARID conversation.
- ARID code paths and generated artifacts listed in the frontmatter.
