---
title: Morgan fingerprint antibiotic-activity oracle baseline
date: 2026-08-11
project: arid
agent: Codex
status: complete
sources:
  - user-provided training command and result
  - /scratch/users/diamant/morgan_FP_abx_classifiers/run_metadata.json
  - /scratch/users/diamant/morgan_FP_abx_classifiers/metrics.json
tags:
  - gneprop
  - antibiotics
  - morgan-fingerprints
  - oracle
  - auprc
---

# Summary

A one-hidden-layer MLP trained on 256-bit Morgan fingerprints from the GNEtolC dataset reached a best validation AUPRC of **0.15263** on predefined split 0. The selected checkpoint is `/scratch/users/diamant/morgan_FP_abx_classifiers/best-v11.ckpt`. The validation set has 374 actives among 16,364 molecules (2.2855% prevalence), so this AUPRC is approximately 6.68 times the random-ranking baseline.

# Key Points

- The run uses the unweighted binary cross-entropy loss.
- The best-checkpoint metric was `val/auprc`, rather than validation loss.
- The maximum recorded validation AUPRC was `0.1526295393705368`, at zero-indexed epoch 8.
- Output artifacts are under `/scratch/users/diamant/morgan_FP_abx_classifiers/`.

# Details

The model was trained from the ARID checkout with:

```bash
cd /home/users/diamant/repos/ARID
/home/groups/btrippe/diamant/miniforge/envs/g2pt/bin/python examples/gneprop_morgan_mlp.py \
  --wandb-project simple_morgan_abx_classifier \
  --output-dir /scratch/users/diamant/morgan_FP_abx_classifiers \
  --num-layers 1 \
  --checkpoint-metric val/auprc \
  --batch-size 64 \
  --lr 5e-4
```

The saved run metadata records 256 Morgan-fingerprint bits, hidden width 256, dropout 0.0, split 0, seed 0, and 20 training epochs. The input data were `/scratch/users/diamant/data/gneprop_data/GNEtolC.csv` and `/scratch/users/diamant/data/gneprop_data/dataset_100k_v1.pkl`.

The split contains 75,943 training examples and 16,364 validation examples. The saved validation-AUPRC trajectory peaks at `0.1526295393705368`; its corresponding validation loss is `0.11139035224914551`. This is an empirical result for this fixed split, not an estimate of generalization to another split or assay.

# Related Notes

None.

# Open Questions

- Compare this unweighted baseline with `--loss-weighting train-prevalence` and assess both AUPRC and calibration.
- Evaluate the remaining predefined splits to measure sensitivity to the scaffold-cluster partition.

# Sources

- User-provided training command, checkpoint path, and reported approximate AUPRC in the 2026-08-11 ARID conversation.
- `/scratch/users/diamant/morgan_FP_abx_classifiers/run_metadata.json` for the saved configuration and selected checkpoint path.
- `/scratch/users/diamant/morgan_FP_abx_classifiers/metrics.json` for the validation loss and AUPRC histories.
