---
title: Portable Morgan fingerprint antibiotic-activity oracle
date: 2026-08-11
project: arid
agent: Codex
status: complete
sources:
  - /scratch/users/diamant/morgan_FP_abx_classifiers_v2/run_metadata.json
  - /scratch/users/diamant/morgan_FP_abx_classifiers_v2/metrics.json
  - /scratch/users/diamant/morgan_FP_abx_classifiers_v2/best.ckpt
  - /home/users/diamant/repos/ARID/arid/oracle.py
  - /home/users/diamant/repos/ARID/examples/gneprop_morgan_mlp.py
tags:
  - gneprop
  - antibiotics
  - morgan-fingerprints
  - oracle
  - auprc
  - reinforcement-learning
---

# Summary

A portable replacement for the Morgan-fingerprint antibiotic-activity oracle was trained at `/scratch/users/diamant/morgan_FP_abx_classifiers_v2/best.ckpt`. It reached validation AUPRC **0.15262954** on predefined GNEtolC split 0, exactly reproducing the earlier legacy-checkpoint result. The validation active prevalence is 2.2855%, so the selected model provides approximately 6.68-fold AUPRC enrichment over random ranking.

The checkpoint uses the importable `arid.oracle.MorganFingerprintMLPConfig` rather than a class serialized under a training script's `__main__` module. A separate-process inference check loaded it through `load_morgan_fingerprint_mlp` under PyTorch's safe weights-only loader and produced logits successfully. This v2 checkpoint is the selected oracle for conditional REINFORCE work.

# Key Points

- Selected checkpoint: `/scratch/users/diamant/morgan_FP_abx_classifiers_v2/best.ckpt`.
- Best `val/auprc`: `0.1526295393705368` at zero-indexed epoch 8.
- Validation loss at the selected epoch: `0.11139035224914551`.
- Lowest validation loss was `0.1044374480843544` at epoch 2, but checkpoint selection used AUPRC.
- Architecture: 256-bit radius-2 Morgan fingerprints and one 256-unit ReLU hidden layer, with no dropout.
- Training used unweighted binary cross-entropy, batch size 64, learning rate `5e-4`, seed 0, and 20 epochs.
- Split 0 contains 75,943 training and 16,364 validation examples; the validation set contains 374 actives.

# Details

The model was trained from `/home/users/diamant/repos/ARID` with:

```bash
/home/groups/btrippe/diamant/miniforge/envs/g2pt/bin/python \
  examples/gneprop_morgan_mlp.py \
  --wandb-project simple_morgan_abx_classifier \
  --output-dir /scratch/users/diamant/morgan_FP_abx_classifiers_v2 \
  --num-layers 1 \
  --checkpoint-metric val/auprc \
  --batch-size 64 \
  --lr 5e-4
```

The input files were `/scratch/users/diamant/data/gneprop_data/GNEtolC.csv` and `/scratch/users/diamant/data/gneprop_data/dataset_100k_v1.pkl`. The run used split 0, 256 fingerprint bits, hidden width 256, dropout 0.0, no prevalence weighting, and W&B project `simple_morgan_abx_classifier`.

The validation AUPRC trajectory was:

```text
0.05037, 0.09673, 0.12796, 0.12499, 0.14612,
0.14179, 0.14487, 0.13824, 0.15263, 0.14465,
0.13980, 0.13553, 0.11906, 0.12461, 0.12108,
0.12268, 0.12822, 0.12662, 0.12645, 0.11726
```

The newly refactored reusable implementation lives in `/home/users/diamant/repos/ARID/arid/oracle.py`. The selected checkpoint was verified with:

```python
from pathlib import Path

import torch

from arid.oracle import load_morgan_fingerprint_mlp

oracle = load_morgan_fingerprint_mlp(
    Path("/scratch/users/diamant/morgan_FP_abx_classifiers_v2/best.ckpt"),
    torch.device("cuda"),
)
```

The earlier `/scratch/users/diamant/morgan_FP_abx_classifiers/best-v11.ckpt` has identical recorded metrics but serialized its configuration under `__main__`; v2 supersedes it for programmatic reuse.

# Related Notes

- [Conditional edit validity and reverse-timestep behavior](conditional-edit-validity-and-timestep-behavior.md): Establishes that the pretrained edit policy can make valid conditional edits and motivates using this oracle as its REINFORCE reward.

# Open Questions

- How quickly does the edit policy exploit regions where this classifier is poorly calibrated or outside the GNEtolC training distribution?
- Should oracle score improvements be accompanied by similarity, uncertainty, or applicability-domain diagnostics?
- How stable is oracle ranking across the other predefined scaffold-cluster splits?

# Sources

- `/scratch/users/diamant/morgan_FP_abx_classifiers_v2/run_metadata.json` for the full training configuration and selected checkpoint path.
- `/scratch/users/diamant/morgan_FP_abx_classifiers_v2/metrics.json` for validation AUPRC and loss histories.
- `/scratch/users/diamant/morgan_FP_abx_classifiers_v2/best.ckpt` for the selected model parameters and portable configuration.
- `/home/users/diamant/repos/ARID/arid/oracle.py` and `/home/users/diamant/repos/ARID/examples/gneprop_morgan_mlp.py` for model loading and training behavior.
