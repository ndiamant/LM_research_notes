---
title: Swiss-Prot partition-schedule results
date: 2026-09-10
project: arid
agent: codex
status: draft
sources:
  - /home/users/diamant/repos/ARID/arid/schedules.py
  - /home/users/diamant/repos/ARID/arid/corruption.py
  - /home/users/diamant/repos/ARID/examples/swissprot_proteins.py
  - /home/users/diamant/repos/ARID/examples/swissprot_evaluate.py
  - /home/users/diamant/repos/ARID/submit_partition_ablations.sh
  - /home/users/diamant/repos/ARID/submit_partition_ablation_evaluations.sh
  - /scratch/users/diamant/arid_runs/partition_3_append_25_epochs
  - /scratch/users/diamant/arid_runs/partition_3_random_25_epochs
  - /scratch/users/diamant/arid_runs/partition_2_4_append_matched_25_epochs
  - /scratch/users/diamant/arid_runs/partition_3_append_del_frac0p20_25_epochs
  - /scratch/users/diamant/arid_runs/fixed_10_50_epochs_84m
tags:
  - swissprot
  - partition-schedule
  - arbitrary-order
  - deletion-supervision
  - esmc
  - slurm
---

# Summary

Replacing repeated short forward deletions with a uniformly sampled three-way partition substantially improved 25-epoch Swiss-Prot ARID generation. At the common sampling temperature 0.8, fixed-three append-only partition generation reached ESM-C PPPL 7.89, compared with the previously reported suffix-only result of 10.11 and the 25-epoch fixed-window result of 12.13. A plain autoregressive baseline remains better at PPPL 4.98.

Randomizing the order of the same three chunks imposed only a moderate penalty, to PPPL 8.68. Varying append-only trajectories between two and four chunks while exactly matching expected residue supervision reached PPPL 8.35. Adding 20% separately constructed deletion targets reached PPPL 8.47 and learned useful deletion heads, but reduced unconditional termination to 86%; the current reverse-step semantics do not accommodate cleanup deletions cleanly.

# Key Points

- All main partition comparisons used the full 220,109-sequence training split, 25 epochs, batch size 64, seed 7, the same approximately 20M-parameter model, and offline evaluation on 1,024 generations against 2,048 validation sequences.
- A clean sequence of length `L` is split at uniformly sampled internal cut points. With `K` chunks and uniform sampling over the `K` insertion states plus DONE, expected residue supervision per protein visit is `L / (K + 1)`.
- Fixed `K=3` therefore provides `L/4` expected residue targets per visit, much denser supervision than selecting one short span from the earlier suffix-only trajectory.
- Random three-chunk generation retained most of the append-only gain: PPPL worsened by only 0.79 while length remained similar.
- The variable-count probabilities `P(K=2)=0.375` and `P(K=4)=0.625` exactly preserve `E[1/(K+1)]=1/4`; their PPPL penalty was only 0.46.
- The deletion experiment used decoupled supervision, not interleaved forward noise. Twenty percent of sampled partition states were replaced by a state with 1–8 unigram-junk residues inserted and a target that deletes that junk.
- Sampling temperature strongly traded off sample quality and distribution matching. For the fixed-three append-only model, temperature 0.8 gave PPPL 7.89 and length W1 18.87, while temperature 1.0 gave PPPL 10.77 and length W1 5.68.
- The working interpretation is that this temperature tradeoff is likely generic across the models: higher temperature improves aggregate length/distribution matching while reducing average sample quality.

# Details

## Partition schedule and supervision

For fixed-three append-only training, two cut points are selected uniformly from the clean sequence's internal boundaries. The chunks are deleted right-to-left in the forward trajectory, so the learned reverse process appends them left-to-right. Each trajectory has four states for training purposes: three reverse-insertion targets and one clean-state DONE target.

For random-order training, the same chunks are generated in a freshly sampled permutation. Step 0 still inserts into the empty sequence, while subsequent steps must choose prepend, append, or an internal chunk boundary. The aggregate validation insertion-position accuracy was 82.5% with a gain of 1.87 nats over uniform. Because one of three insertion steps has the trivial empty-sequence position, the implied average accuracy over the two nontrivial steps is approximately 74%, assuming equal representation.

The supervision-matched variable-count experiment used:

```yaml
partition_chunk_counts: [2, 4]
partition_chunk_count_probabilities: [0.375, 0.625]
partition_generation_order: left_to_right
```

This satisfies:

```text
0.375 / (2 + 1) + 0.625 / (4 + 1) = 1 / 4
```

It changes chunk sizes and makes DONE timing variable without changing expected residue supervision.

## Offline results at temperature 0.8

All rows below use the shared `examples/swissprot_evaluate.py` protocol. Lower is better for PPPL, length W1, k-mer JS, and embedding Frechet.

| Model | PPPL | PLL | Length W1 | Median length | Termination | k-mer-2 JS | Frechet |
|---|---:|---:|---:|---:|---:|---:|---:|
| Autoregressive baseline | **4.98** | -1.605 | 6.83 | 156 | 1.000 | **0.00545** | **0.688** |
| Fixed-three append-only partition | **7.89** | **-2.065** | 18.87 | 175 | 1.000 | **0.00822** | **1.357** |
| Variable 2/4 append-only partition | 8.35 | -2.122 | 19.75 | 178 | 1.000 | 0.00880 | 1.455 |
| Fixed-three + 20% deletion supervision | 8.47 | -2.137 | **12.22** | **160** | 0.860 | 0.01091 | 1.521 |
| Fixed-three random-order partition | 8.68 | -2.161 | 21.09 | 174.5 | 1.000 | 0.01209 | 1.564 |
| Fixed-10 interleaved, 500 epochs | 9.21 | -2.220 | 4.70 | 150 | 0.997 | 0.00930 | 1.681 |
| Previously reported suffix-only ARID | 10.11 | -2.313 | 7.6 | not recorded here | 1.000 | 0.0085 | 1.846 |
| Fixed-10 interleaved, 84M model, 50 epochs | 11.65 | -2.456 | 4.90 | 150 | 1.000 | 0.00957 | 2.164 |
| Fixed-10 interleaved, 50 epochs | 11.72 | -2.462 | 5.68 | 142 | 1.000 | 0.00981 | 2.221 |
| Fixed-10 interleaved, 25 epochs | 12.13 | -2.496 | 4.99 | 143 | 1.000 | 0.01023 | 2.279 |

The 500-epoch fixed-window row is not compute-matched; it is included to show that the 25-epoch partition models also exceeded that much longer run on PPPL.

Within the earlier fixed-window family, extending the approximately 20M model from 25 to 50 epochs improved PPPL from 12.13 to 11.72. Increasing capacity to approximately 84M parameters at 50 epochs changed PPPL only slightly further, to 11.65, whereas training the 20M model for 500 epochs reached 9.21. Those results motivated testing denser per-visit residue supervision rather than capacity alone.

The random-order penalty relative to fixed-three append-only was 0.096 nats per residue, or 0.79 PPPL. Length median was essentially unchanged, so the PPPL difference is not explained by a new length shift. The random-order span-residue gain over the unigram baseline was 0.201 nats versus 0.268 for append-only, suggesting that arbitrary-order conditioning is harder but tractable with three large chunks.

## Temperature sweep

The fixed-three append-only checkpoint was sampled over temperatures 0.6–1.0. ESM-C was run for temperatures 0.8 and 1.0; the other rows are distribution-only evaluations.

| Temperature | Length W1 | Median | Mean | k-mer-2 JS | ESM-C PPPL |
|---:|---:|---:|---:|---:|---:|
| 0.6 | 33.30 | 184 | 183.99 | 0.06674 | not scored |
| 0.7 | 28.52 | 187 | 179.16 | 0.02991 | not scored |
| 0.8 | 18.87 | 175 | 169.19 | 0.00822 | **7.89** |
| 0.9 | 11.39 | 163 | 160.99 | **0.00344** | not scored |
| 1.0 | **5.68** | **154** | **153.12** | 0.00376 | 10.77 |

The validation reference median is 151 residues. Higher temperature makes the span STOP token relatively more competitive, shortening generations and improving aggregate distribution matching, while noisier residue choices worsen ESM-C sample quality.

## Decoupled deletion supervision

The deletion-enabled run retained a forward trajectory containing only deletions of clean chunks; those operations create the reverse-insertion targets. After sampling a trajectory state, the collator replaced its normal target with probability 0.2 by inserting a fresh 1–8-residue unigram-junk span and asking the model to delete it. Junk never enters the trajectory and therefore cannot contaminate later insertion targets.

For fixed `K=3`, expected target proportions are 60% reverse insertion, 20% DONE, and 20% reverse deletion. Because deletion rows replace normal rows, the experiment has fixed optimizer compute but 20% less insertion-residue exposure than the append-only control.

The deletion heads learned a nontrivial task:

- deletion-start accuracy 30.6%, with a 0.482-nat gain over a uniform legal-position baseline;
- deletion-length accuracy 38.3%, versus 12.5% chance over eight lengths;
- PPPL worsened only 0.58 relative to append-only, while length W1 improved from 18.87 to 12.22.

The 86% termination rate prevents treating the length improvement as an unqualified gain. With `max_reverse_steps: 4`, three intended insertions, possible cleanup deletions, and DONE no longer fit reliably. More fundamentally, deleting injected junk restores the underlying partition state and should not necessarily advance partition progress, but the current sampler increments the reverse-step index after every edit. This progress/edit distinction should be fixed before drawing strong conclusions about deletion-enabled unconditional generation.

## Command and Slurm provenance

The runtime environment and model cache were:

```bash
source /home/groups/btrippe/diamant/miniforge/etc/profile.d/mamba.sh
mamba activate esm3
export HF_HOME=/scratch/users/diamant/models
cd /home/users/diamant/repos/ARID
```

The fixed-three append-only and random-order runs used their scratch configs directly:

```bash
python examples/swissprot_proteins.py --config /scratch/users/diamant/arid_runs/partition_3_append_25_epochs/config.yaml
python examples/swissprot_proteins.py --config /scratch/users/diamant/arid_runs/partition_3_random_25_epochs/config.yaml
```

The variable-count and deletion ablations were submitted through `submit_partition_ablations.sh`:

```bash
sbatch submit_partition_ablations.sh
```

- Training array job: `42760632`; both tasks completed on `akundaje` with exit code 0, in 2:22:57 and 2:20:23.
- Evaluation array job: `42783746`; both tasks completed on `akundaje` with exit code 0, in 2:55 each.
- Both scripts requested one GPU satisfying `GPU_MEM:48GB|GPU_MEM:80GB`, with eligibility on `hns,akundaje`.
- Slurm logs live under `/scratch/users/diamant/arid_runs/slurm_logs/`.

The evaluation command embodied by `submit_partition_ablation_evaluations.sh` was:

```bash
python examples/swissprot_evaluate.py \
  --checkpoint "$run_dir/checkpoints/model.pt" \
  --splits /scratch/users/diamant/datasets/swiss_prot/splits/splits.json \
  --output-dir "$run_dir/evaluation" \
  --gpu 0
```

Focused validation before submission comprised 16 passing partition/corruption unit tests plus a real ESM3-environment forward/backward smoke test. The smoke batch contained 11 deletion targets among 64 examples and produced a finite loss.

At the time of this note, the SwissProt partition integration and deletion support in `arid/corruption.py`, `arid/schedules.py`, and `examples/swissprot_proteins.py`, as well as the two Slurm scripts, were uncommitted changes on ARID branch `nd/chunk-partition-and-AOE-model` at base commit `d3e2522`. The runs therefore depend on the working-tree versions of those files and should be committed before treating the code provenance as durable.

# Related Notes

- [Swiss-Prot setup and matched baseline runs](../2026-09-06/swissprot-setup-and-baseline-runs.md): Records the Sherlock environment, Swiss-Prot split, ESM-C cache, and original matched baseline setup reused here.
- [Conditional edit validity and reverse-timestep behavior](../2026-08-11/conditional-edit-validity-and-timestep-behavior.md): Earlier evidence that reverse-step conditioning can dominate DONE behavior during conditional editing.
- [Selected MOSES constant-then-delete schedule](../2026-07-31/selected-moses-constant-then-delete-schedule.md): Documents the fixed schedule that motivated the earlier Swiss-Prot fixed-window comparison.

# Open Questions

- Does combining variable chunk count with random generation order preserve the current PPPL gains?
- How should partition progress be represented so cleanup deletions do not consume an insertion-progress step?
- After fixing progress semantics, does deletion supervision retain its improved length calibration and useful deletion accuracy with full termination?
- Would a 32-epoch deletion run, approximately matching insertion-residue exposure rather than optimizer compute, recover the append-only span loss?
- How much does the reverse-step embedding contribute when chunk count varies, and should it be removed for models intended for RL editing?
- Are the schedule differences stable over multiple training and sampling seeds?

# Sources

- ARID implementation: `/home/users/diamant/repos/ARID/arid/schedules.py`, `arid/corruption.py`, `examples/swissprot_proteins.py`, and `examples/swissprot_evaluate.py`.
- Run configs, resolved configs, checkpoints, training metrics, and offline evaluations under `/scratch/users/diamant/arid_runs/{partition_3_append_25_epochs,partition_3_random_25_epochs,partition_2_4_append_matched_25_epochs,partition_3_append_del_frac0p20_25_epochs}/`.
- Earlier baseline evaluation artifacts under `/scratch/users/diamant/arid_runs/{autoregressive_repro,fixed10,fixed_10_50_epochs,fixed_10_50_epochs_84m,fixed_10_500_epochs}/`.
- Slurm submission scripts `/home/users/diamant/repos/ARID/submit_partition_ablations.sh` and `/home/users/diamant/repos/ARID/submit_partition_ablation_evaluations.sh`; accounting records for jobs `42760632` and `42783746`.
- User-provided interpretation that the temperature tradeoff likely generalizes across models.
