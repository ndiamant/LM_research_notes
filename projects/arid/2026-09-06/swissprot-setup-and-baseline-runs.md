---
title: Swiss-Prot setup and matched baseline runs
date: 2026-09-06
project: arid
agent: codex
status: draft
sources:
  - /home/users/diamant/repos/ARID/examples/swissprot_autoregressive.py
  - /home/users/diamant/repos/ARID/examples/swissprot_proteins.py
  - /home/users/diamant/repos/ARID/examples/swissprot_evaluate.py
  - /scratch/users/diamant/arid_runs/autoregressive_repro/checkpoints/run_metadata.json
  - /scratch/users/diamant/arid_runs/fixed10/config.yaml
  - /scratch/users/diamant/arid_runs/fixed10/outputs/resolved_config.yaml
tags:
  - swissprot
  - autoregressive
  - fixed-schedule
  - esmc
  - wandb
  - sherlock
---

# Summary

The Swiss-Prot environment and data were made portable from the original `/backups/chihoim` setup to Sherlock scratch. Two matched 25-epoch, full-data baselines were started on 2026-09-06: a 19.64M-parameter left-to-right autoregressive model and a 20.33M-parameter ARID model using a fixed 10-step interleaved forward schedule. They share the 220,109-sequence training split, batch size 64, learning rate `3e-4`, weight decay `1e-2`, and seed 7.

# Key Points

- Python environment: `/home/groups/btrippe/diamant/miniforge/envs/esm3/bin/python`.
- Data: `/scratch/users/diamant/datasets/swiss_prot/splits/`, containing the expected 220,109/12,251/12,233 train/validation/test split.
- ESM-C 600M cache: `/scratch/users/diamant/models/models--EvolutionaryScale--esmc-600m-2024-12` (about 2.2 GB).
- W&B project: `arid-swissprot`.
- Autoregressive W&B run ID: `ctqnmyv0`; display name `260901_ar_baseline_repro`.
- Fixed-window ARID W&B run ID: `4s3fqolb`; display name `swissprot-fixed10-seed7`.
- The autoregressive script does not evaluate ESM-C during training. It generates 1,024 proteins after training for evaluation with the shared offline evaluator.
- Fixed-window ARID evaluates 64 generated proteins with ESM-C every five epochs; final comparisons should nevertheless use the shared offline evaluator for both runs.

# Details

## Environment and assets

The existing `esm3` environment was reused rather than modifying the Python 3.10 SMILES environment. The verified package state was:

```text
Python             3.12.11
torch              2.6.0 (CUDA 12.4)
pytorch-lightning  2.6.5
esm                3.2.0
wandb              0.29.0
numpy              2.3.1
scipy              1.16.0
matplotlib         3.10.3
PyYAML             6.0.2
```

ARID is installed editable from `/home/users/diamant/repos/ARID`. `pytorch-lightning==2.6.5` was added for the training entry points. A direct `pip install wandb` selected the W&B 0.29.0 source distribution because the CentOS 7 host exposes glibc 2.17 while its pip wheel requires manylinux 2.28. That source build failed with the available Go 1.18 compiler. The working installation came from conda-forge instead:

```bash
/home/groups/btrippe/diamant/miniforge/bin/mamba install \
  -n esm3 -c conda-forge 'wandb=0.29.0'
```

The environment currently has a broken optional `charset_normalizer` native import after mixing package managers. `requests` warns and falls back; both training and W&B tracking are functioning, so it was deliberately left unchanged during the runs.

The Swiss-Prot split and filtered FASTA were copied from `/backups/chihoim`. `splits.json` now resolves its filtered FASTA locally:

```text
/scratch/users/diamant/datasets/swiss_prot/splits/splits.json
/scratch/users/diamant/datasets/swiss_prot/splits/swissprot_max256.fasta
```

The recorded split counts are 220,109 train, 12,251 validation, and 12,233 test. The ESM-C cache was also copied from `/backups/chihoim` and moved into the common Hugging Face cache root. ESM-C 600M was successfully loaded offline from that cache with 575,036,992 parameters. Jobs should export:

```bash
export HF_HOME=/scratch/users/diamant/models
export HF_HUB_CACHE=/scratch/users/diamant/models
```

## Autoregressive baseline

The launch command captured by W&B was:

```bash
cd /home/users/diamant/repos/ARID
python examples/swissprot_autoregressive.py \
  --splits /scratch/users/diamant/datasets/swiss_prot/splits/splits.json \
  --output-dir /scratch/users/diamant/arid_runs/autoregressive_repro/checkpoints \
  --analysis-dir /scratch/users/diamant/arid_runs/autoregressive_repro/outputs \
  --gpu 0 \
  --wandb-project arid-swissprot \
  --wandb-name 260901_ar_baseline_repro \
  --wandb-tags autoregressive baseline order-diagnostic
```

The executable was the `esm3` environment's Python. The run uses all training and validation sequences, 25 epochs, batch 64, 3,440 optimizer steps per epoch, seed 7, `d_model=384`, 11 transformer layers, and 19,636,632 parameters. It ran on an H100 80 GB in interactive Slurm allocation `42031803` on `sh04-03n07`.

At approximately 14:14 Pacific on 2026-09-06, it had completed epoch 14 and was 9% through epoch 15. The latest recorded validation loss was about 2.47. W&B run ID: `ctqnmyv0`.

The script writes `generated.txt` only after training finishes. ESM-C and distributional evaluation should then run as:

```bash
export HF_HOME=/scratch/users/diamant/models
export HF_HUB_CACHE=/scratch/users/diamant/models

/home/groups/btrippe/diamant/miniforge/envs/esm3/bin/python \
  examples/swissprot_evaluate.py \
  --sequences-file /scratch/users/diamant/arid_runs/autoregressive_repro/outputs/generated.txt \
  --splits /scratch/users/diamant/datasets/swiss_prot/splits/splits.json \
  --output-dir /scratch/users/diamant/arid_runs/autoregressive_repro/outputs/evaluation \
  --reference-samples 2048 \
  --esmc-batch-size 16 \
  --gpu 0
```

## Fixed-10 ARID baseline

The complete input and resolved configurations are retained at:

```text
/scratch/users/diamant/arid_runs/fixed10/config.yaml
/scratch/users/diamant/arid_runs/fixed10/outputs/resolved_config.yaml
```

The launch command captured by W&B was:

```bash
cd /home/users/diamant/repos/ARID
export HF_HOME=/scratch/users/diamant/models
export HF_HUB_CACHE=/scratch/users/diamant/models
/home/groups/btrippe/diamant/miniforge/envs/esm3/bin/python \
  examples/swissprot_proteins.py \
  --config /scratch/users/diamant/arid_runs/fixed10/config.yaml
```

The run uses all training and validation sequences, 25 epochs, batch 64, learning rate `3e-4`, seed 7, `d_model=384`, six encoder layers, four decoder layers, and approximately 20.33M parameters. It ran on an L40S in interactive Slurm allocation `42277809` on `sh04-07n09`. W&B run ID: `4s3fqolb`.

Its forward schedule isolates fixed versus length-scaled windowing while retaining the Swiss-Prot span sizes:

```text
schedule_type                    constant
mixed_steps                     10
insert_probability              0.5
max_forward_delete_length       24
max_forward_insert_length       8
forward_delete_length_decay     1.0
forward_insert_length_decay     1.0
deletion_example_fraction       0.0
max_corrupted_length            352
max_reverse_steps               128
```

Thus steps 1–10 choose forward insertion or deletion with equal probability, and steps 11 onward are deletion-only. Span lengths are uniform. The fixed run uses no decoupled deletion examples.

At approximately 14:14 Pacific on 2026-09-06, it had completed epoch 3 and was 6% through epoch 4. Epochs took about 255 seconds. The first periodic generation and ESM-C callback was therefore still pending at the time of inspection.

Both W&B metadata records point to ARID commit `3766b7642004da105c7d1bef33e3156fe73fddab`. The working tree also contained an uncommitted change exposing `--schedule-type` and `--mixed-steps` as Swiss-Prot CLI overrides. That change did not determine this run: the fixed schedule was selected through `config.yaml`, which the base commit already supported.

# Related Notes

- [Selected MOSES constant-then-delete schedule](../2026-07-31/selected-moses-constant-then-delete-schedule.md): Establishes the fixed 10-step 50/50 schedule transferred here, with smaller MOSES spans of forward-delete 12 and forward-insert 4.
- [MOSES H100 full training reaches perfect validation sample validity](../2026-08-05/moses-h100-full-training-validity.md): Provides prior full-training context for the SMILES branch and W&B-based generation monitoring.

# Open Questions

- Whether fixed-window interleaving improves final Swiss-Prot ESM-C PPPL relative to the length-scaled and decoupled-deletion ARID runs.
- How strongly the fixed 10-step window under-supplies reverse-deletion targets for long proteins.
- Whether the two runs finish cleanly and produce their expected final artifacts.
- Final comparison should use the same 1,024-generation offline evaluator for both runs; the ARID callback's 64-sample metrics are monitoring only.

# Sources

- ARID scripts, run metadata, resolved configuration, W&B metadata, logs, and cache/data paths listed in the frontmatter and details.
- [Selected MOSES constant-then-delete schedule](../2026-07-31/selected-moses-constant-then-delete-schedule.md).
