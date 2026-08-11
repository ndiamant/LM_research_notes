---
title: Conditional edit validity and reverse-timestep behavior
date: 2026-08-11
project: arid
agent: Codex
status: complete
sources:
  - /home/users/diamant/repos/ARID/examples/moses_conditional_validity.py
  - /scratch/users/diamant/arid_moses_conditional_validity/h100_full_d384_l3_b512_w4/summary.json
  - user-provided commands and unforced-run results from the 2026-08-11 ARID conversation
tags:
  - smiles
  - conditional-generation
  - validity
  - reinforcement-learning
  - reverse-timestep
---

# Summary

The trained `h100_full_d384_l3_b512_w4` ARID model can make highly valid edits when initialized from an already-valid MOSES molecule. Without intervention, however, the mode head treats reverse step 0 differently from all other tested steps: step 0 often edits, while steps 1–16 almost always immediately predict `DONE`. Forcing only the first action to be a nonempty edit raises the final changed-molecule fraction to 85.1–91.4% across all tested starting steps while retaining 98.0–99.2% overall final validity.

This supports using a forced first edit for the planned conditional REINFORCE experiment while continuing to increment the reverse-step counter normally. The model then usually stops naturally after 1–2 edits. Exact returns to the starting token sequence occur in roughly 9–15% of forced trajectories and can serve as meaningful no-improvement samples.

# Key Points

- Model: `/scratch/users/diamant/arid_moses_training/h100_full_d384_l3_b512_w4`.
- Evaluation used 64 valid MOSES test molecules and 128 trajectories per molecule, giving 8,192 trajectories per condition at temperature 0.8.
- Without forcing, reverse step 0 produced a 44.1% changed-output rate at an eight-edit cap; steps 1–16 produced only 1.5–2.7% changed outputs.
- With a forced first edit, all trajectories made at least one edit and averaged 1.21–1.40 edits before stopping.
- Forced-run changed-output validity was 97.7–99.1%, depending on the starting reverse step.
- At least 99.9% of forced trajectories terminated naturally before the eight-edit cap.
- `changed` compares final and initial token sequences. Because the first edit was guaranteed, every unchanged final output in the forced run had edited and then returned exactly to its starting sequence.

# Details

## Unforced first edit

The unforced diagnostic was run from `/home/users/diamant/repos/ARID` with:

```bash
/home/groups/btrippe/diamant/miniforge/envs/g2pt/bin/python \
  examples/moses_conditional_validity.py \
  --num-starting-molecules 64 \
  --samples-per-molecule 128 \
  --max-edits 4 8 \
  --start-reverse-steps 0 1 2 4 8 16 \
  --device cuda
```

The eight-edit-cap results are shown below. The later forced run reused the output directory and overwrote the unforced `summary.json`; these values were captured from the earlier summary during the 2026-08-11 ARID conversation.

| Start reverse step | Final validity | Changed fraction | Validity among changed | Mean edits given at least one edit |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 98.77% | 44.14% | 97.21% | 1.422 |
| 1 | 99.96% | 1.49% | 97.54% | 1.122 |
| 2 | 100.00% | 2.01% | 100.00% | 1.042 |
| 4 | 100.00% | 1.87% | 100.00% | 1.052 |
| 8 | 100.00% | 2.67% | 100.00% | 1.082 |
| 16 | 99.99% | 2.08% | 99.41% | 1.064 |

Changing the maximum from four to eight edits had little effect. The model stopped naturally rather than reaching the cap. The behavior suggests that reverse step 0 acts as a special edit regime because `DONE` is absent there during diffusion training; a clean-looking input at later steps strongly elicits `DONE`.

## Forced first edit

The forced diagnostic used:

```bash
/home/groups/btrippe/diamant/miniforge/envs/g2pt/bin/python \
  examples/moses_conditional_validity.py \
  --num-starting-molecules 64 \
  --samples-per-molecule 128 \
  --max-edits 8 \
  --start-reverse-steps 0 1 2 4 8 16 \
  --force-first-edit \
  --device cuda
```

`--force-first-edit` masks `DONE` on the first decision and prevents an insertion from immediately generating its span-stop token. Subsequent decisions use the normal mode distribution, and the reverse-step index increments after every edit.

| Start reverse step | Final validity | Changed fraction | Validity among changed | Mean edits | Natural termination |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 98.02% | 85.10% | 97.68% | 1.399 | 99.93% |
| 1 | 98.88% | 90.28% | 98.76% | 1.264 | 99.94% |
| 2 | 99.07% | 91.00% | 98.98% | 1.251 | 99.93% |
| 4 | 99.00% | 91.39% | 98.90% | 1.246 | 99.93% |
| 8 | 99.21% | 88.54% | 99.10% | 1.304 | 99.90% |
| 16 | 98.96% | 90.99% | 98.86% | 1.213 | 99.96% |

The changed fraction remains below one because later edits can exactly undo the forced first edit. This is useful for conditional REINFORCE: a returned starting molecule receives zero activity improvement and behaves as a natural no-op comparison within its molecule-specific rollout group.

# Related Notes

- [MOSES H100 full training reaches perfect validation sample validity](../2026-08-05/moses-h100-full-training-validity.md): Documents the training run and unconditional validity evidence for the checkpoint evaluated here.

# Open Questions

- Does forced-first-edit behavior remain stable for other seeds, temperatures, and starting-molecule sets?
- How does validity change when REINFORCE shifts probability toward longer trajectories or particular edit types?
- Should conditional rollout outputs use unique run directories so forced and unforced summaries remain independently auditable?

# Sources

- `/home/users/diamant/repos/ARID/examples/moses_conditional_validity.py` for the diagnostic and forced-edit semantics.
- `/scratch/users/diamant/arid_moses_conditional_validity/h100_full_d384_l3_b512_w4/summary.json` for the completed forced-run configuration and metrics.
- Earlier unforced summary output captured in the 2026-08-11 ARID conversation before the shared output directory was overwritten.
