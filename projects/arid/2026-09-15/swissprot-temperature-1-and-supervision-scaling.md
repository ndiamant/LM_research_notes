---
title: Swiss-Prot temperature-1 results and append-3 supervision scaling
date: 2026-09-15
project: arid
agent: codex
status: draft
sources:
  - /home/users/diamant/repos/ARID/examples/swissprot_evaluate.py
  - /home/users/diamant/repos/ARID/examples/swissprot_autoregressive.py
  - /home/users/diamant/repos/ARID/examples/indigo_insertion_synthetic.py
  - /home/users/diamant/repos/ARID/arid/corruption.py
  - /home/users/diamant/repos/ARID/arid/indigo.py
  - /home/users/diamant/repos/ARID/submit_swissprot_ablation_evaluations_temp1.sh
  - /home/users/diamant/repos/ARID/submit_swissprot_autoregressive_eval_temp1.sh
  - /home/users/diamant/repos/ARID/submit_partition_3_append_75_epochs.sh
  - /home/users/diamant/repos/ARID/submit_partition_3_append_75_epochs_eval_temp1.sh
  - /scratch/users/diamant/arid_runs/partition_3_append_75_epochs
  - /scratch/users/diamant/arid_runs/autoregressive_repro
tags:
  - swissprot
  - partition-schedule
  - autoregressive
  - supervision-density
  - temperature
  - esmc
  - indigo
---

# Summary

Swiss-Prot evaluation is now standardized at sampling temperature 1.0, using
1,024 generated sequences and 2,048 held-out references. The 25-epoch,
three-partition append-only model remained the best 25-epoch ARID variant at
ESM-C PPPL 10.77. Training the same approximately 20M-parameter configuration
for 75 epochs improved PPPL to 8.83, versus 7.45 for the 25-epoch left-to-right
autoregressive baseline. Length W1 remained essentially unchanged at 5.64.

This closes about 59% of the append-3-to-AR PPPL difference and supports the
hypothesis that sparse residue supervision explains much of ARID's deficit.
It does not close the gap completely: append-3 remains 1.37 PPPL, or 18%, above
AR, with worse ESM-C embedding distances.

An implementation audit corrects the recent informal supervision arithmetic.
Fixed-three partition training samples uniformly from three insertion states
and one clean `DONE` state. Expected residue supervision is therefore `L/4` per
sequence visit, not `L/3`. The 75-epoch run supplies about 18.75 AR-equivalent
residue epochs; a strict residue-supervision match to 25 AR epochs is 100
append-3 epochs. The `matched-residue-supervision` W&B tag on the 75-epoch run
records the original intent but is not mathematically exact.

# Key Points

- At temperature 1.0, the ordering is AR (7.45 PPPL), append-3 at 75 epochs
  (8.83), append-3 at 25 epochs (10.77), variable 2/4 append-only (11.22),
  decoupled deletion supervision (11.42), and three-partition random order
  (11.76).
- Extending append-3 from 25 to 75 epochs improved PPPL by 18%, embedding MMD²
  from 0.854 to 0.611, and embedding Fréchet from 1.868 to 1.458.
- Length calibration did not cause the quality improvement: length W1 was
  5.68 at 25 epochs and 5.64 at 75 epochs.
- The 75-epoch model has 100% reported termination, uniqueness, and novelty in
  the 1,024-sample evaluation.
- The 25-epoch random-order penalty is 1.00 PPPL relative to append-only. This
  suggests arbitrary chunk order is harder but is not the dominant source of
  the original gap.
- Neither an approximately 84M model at 50 epochs nor cleanup deletions in the
  forward trajectory helped sample quality. Extra training of the 20M model
  and denser residue supervision have been more effective than capacity alone.
- The historical targets supplied by the user are PPPL 5.11 ideally and below
  11.29 to exceed the best prior coworker schedule. These targets may use a
  different temperature protocol, so the temperature-1 table is the clean
  comparison within the current experiment series.

# Details

## Temperature-1 offline results

All rows use the shared offline evaluator and are ordered by generated ESM-C
PPPL. Lower is better for PPPL, length W1, MMD², and Fréchet.

| Model | Epochs | ESM-C PPPL | Length W1 | Termination | MMD² | Fréchet |
|---|---:|---:|---:|---:|---:|---:|
| Autoregressive baseline | 25 | **7.455** | **4.026** | 0.999* | **0.441** | **1.215** |
| Three-partition append-only | 75 | **8.826** | 5.641 | 1.000 | 0.611 | 1.458 |
| Three-partition append-only | 25 | 10.765 | 5.684 | 1.000 | 0.854 | 1.868 |
| Variable 2/4-partition append-only | 25 | 11.223 | 6.523 | 1.000 | 0.878 | 1.948 |
| Append-only + 20% decoupled deletion supervision | 25 | 11.416 | 15.234 | 0.819 | 0.908 | 2.005 |
| Three-partition random order | 25 | 11.760 | 6.423 | 1.000 | 0.972 | 2.069 |
| Fixed-10 interleaved | 500 | 12.216 | 7.448 | 0.997 | 1.018 | 2.191 |
| Three-partition append + coupled forward deletion | 25 | 12.227 | 9.387 | 1.000 | 0.996 | 2.164 |
| Fixed-10 interleaved | 50 | 13.848 | 11.263 | 1.000 | 1.074 | 2.516 |
| Fixed-10 interleaved, approximately 84M | 50 | 13.870 | 6.222 | 1.000 | 1.083 | 2.470 |
| Fixed-10 interleaved | 25 | 14.378 | 11.622 | 1.000 | 1.089 | 2.586 |

`*` The AR generation summary records 1,023/1,024 sequences terminating
(0.9990). The file-based common evaluator reports 1.0 because a plain sequence
file does not retain termination flags.

The ESM-C pseudo-log-likelihood improved from -2.376 at 25 append-3 epochs to
-2.178 at 75 epochs; AR is -2.009. Expressed in PPPL, the append-3 excess over
AR fell from 3.311 to 1.372, closing 58.6% of the initial difference. The
improvement is also reflected in embedding distances, while k-mer-2 JS stayed
near 0.0037 and length statistics barely changed.

The final ten validation span-residue losses fluctuate around 2.44 and training
span-residue loss around 2.12. This apparent plateau coincides with the cosine
learning rate approaching zero, so it does not by itself establish an
architectural ceiling. A fresh 100-epoch schedule would be the strict
residue-supervision-matched comparison; merely resuming the final 75-epoch
checkpoint at its near-zero learning rate would not be equivalent.

## Supervision accounting

For `K=3`, each trajectory contains:

- three reverse-insertion states, each with expected chunk length `L/3`;
- one clean state whose target is `DONE` and contains no residue target.

The low-discrepancy state sampler is uniform over all four states. Therefore:

```text
E[residue targets per visit]
  = 3 * P(insertion state) * E[chunk length]
  = 3 * (1/4) * (L/3)
  = L/4
```

The approximate AR-equivalent residue epochs are consequently:

| Append-3 epochs | AR-equivalent residue epochs |
|---:|---:|
| 25 | 6.25 |
| 75 | 18.75 |
| 100 | 25.00 |
| 150 | 37.50 |

This accounting strengthens, rather than weakens, the main empirical result:
the 75-epoch append-3 model reached PPPL 8.83 despite still receiving about 25%
less expected residue supervision than the 25-epoch AR model. It remains an
approximation because AR predicts every residue plus EOS on each visit, whereas
ARID also spends capacity and gradient budget on mode, span-stop, and position
decisions.

## Implications for dense any-order autoregression

The developing INDIGO synthetic model generates one token per insertion and
supervises an entire randomly ordered insertion trajectory in one causal pass.
It therefore offers approximately full token supervision per sequence visit
while retaining arbitrary-order insertion. The current results make this a
high-value next direction: append-3 improved sharply when residue exposure was
increased, while randomizing three chunk positions imposed a comparatively
moderate penalty.

This conclusion remains a hypothesis until INDIGO is trained on Swiss-Prot.
Important remaining differences from left-to-right AR include the much larger
space of partial contexts, a position target for every residue, and the current
factorization `p(token | state) p(position | token, state)`. For proteins,
predicting position first and then conditioning residue identity on that gap
may be more sample-efficient because it exposes the relevant left and right
context before residue prediction.

## Command and Slurm provenance

The runtime setup was:

```bash
source /home/groups/btrippe/diamant/miniforge/etc/profile.d/mamba.sh
mamba activate esm3
export HF_HOME=/scratch/users/diamant/models
cd /home/users/diamant/repos/ARID
```

The temperature-1 ablation array was submitted with:

```bash
sbatch submit_swissprot_ablation_evaluations_temp1.sh
```

Job `43424089` completed all eight tasks successfully. The AR checkpoint needed
a separate generation step and was submitted with:

```bash
sbatch submit_swissprot_autoregressive_eval_temp1.sh
```

Job `43426914` completed in 2:10. Its output is under
`/scratch/users/diamant/arid_runs/autoregressive_repro/eval_temp_1.0/`.

The 75-epoch append-3 run used an exact copy of the 25-epoch model, schedule,
loss, data, optimizer, and seed, changing only epoch count and run metadata:

```bash
sbatch submit_partition_3_append_75_epochs.sh
```

Training job `43431533` completed successfully in 6:50:41. The checkpoint and
training metrics are under
`/scratch/users/diamant/arid_runs/partition_3_append_75_epochs/`. Its evaluation
was submitted with:

```bash
sbatch submit_partition_3_append_75_epochs_eval_temp1.sh
```

Evaluation job `43579581` completed successfully in 3:16. The primary result is
`/scratch/users/diamant/arid_runs/partition_3_append_75_epochs/eval_temp_1.0/metrics.json`.

The underlying evaluation command was:

```bash
python examples/swissprot_evaluate.py \
  --checkpoint /scratch/users/diamant/arid_runs/partition_3_append_75_epochs/checkpoints/model.pt \
  --splits /scratch/users/diamant/datasets/swiss_prot/splits/splits.json \
  --output-dir /scratch/users/diamant/arid_runs/partition_3_append_75_epochs/eval_temp_1.0 \
  --temperature 1.0 \
  --gpu 0
```

# Related Notes

- [Swiss-Prot partition-schedule results](../2026-09-10/swissprot-partition-schedule-results.md): Defines the partition schedules, reports the temperature-0.8 experiments, and documents the original deletion ablations.
- [Swiss-Prot setup and matched baseline runs](../2026-09-06/swissprot-setup-and-baseline-runs.md): Records the environment, dataset split, ESM-C cache, and initial AR/fixed-schedule training setup reused here.

# Open Questions

- Does a fresh 100-epoch append-3 run, which strictly matches expected residue
  supervision to 25 AR epochs, close the remaining PPPL gap?
- Is the final gap due primarily to the span decoder architecture, span-stop
  prediction, or residual differences in conditioning and optimization?
- Can Swiss-Prot INDIGO approach AR quality with full token supervision while
  retaining arbitrary-order insertion?
- For dense insertion modeling, is position-first factorization more efficient
  than the current token-first factorization?
- Are the rankings stable across training seeds and independent sets of 1,024
  generated samples?

# Sources

- ARID implementation and Slurm scripts listed in the frontmatter.
- Temperature-1 evaluation artifacts under
  `/scratch/users/diamant/arid_runs/*/eval_temp_1.0/metrics.json`.
- The 75-epoch config, checkpoint, training metrics, and evaluation under
  `/scratch/users/diamant/arid_runs/partition_3_append_75_epochs/`.
- Slurm accounting records for jobs `43424089`, `43426914`, `43431533`, and
  `43579581`.
- User-provided historical PPPL targets and experiment interpretation from the
  2026-09 ARID conversation.
