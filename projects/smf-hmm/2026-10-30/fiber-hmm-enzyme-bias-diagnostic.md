---
title: Fiber-seq HMM absorbs m6A context bias into its states and its DNA prior
date: 2026-10-30
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - Run /scratch/users/diamant/2026-10-22_fiber_fixed_blocks_global_kl3 (train.sbatch, slurm-49057551.log, checkpoints/best, checkpoints/best.model_spec.json; wandb run dhpxj0qy)
  - /scratch/users/diamant/kernel_denoise_proto/hmm_bias_check.py and hmm_bias_check_60x8.txt
  - smf_hmm repo, evaluation/datasets/smf/context_bias/adj_table_3_3.npz (branch nd/experiments, commit d20f916)
  - smf_hmm repo, smf_hmm/emission_modules.py (OrderedGroupSMFEmissions), smf_hmm/hmm.py, smf_hmm/checkpointing.py
  - User discussion, 2026-10-30
tags: [smf-hmm, fiberseq, m6a, gm12878, enzyme-bias, sequence-context, emissions, cnn-prior, diagnostics]
---

# Summary

The user's fixed-blocks HMM on GM12878 Fiber-seq (`2026-10-22_fiber_fixed_blocks_global_kl3`)
"didn't work very well". Its emissions are one constant P(protected | state) at every site, but
the [context-bias note](fiberseq-m6a-sequence-context-bias.md) found about a 100-fold spread in
the m6A rate across 3+3 contexts. **The model compensated by repurposing states, and its
DNA-conditioned CNN prior learned the enzyme/caller's sequence preference.**

We ran the checkpoint on 60 validation windows (492K called sites) and binned sites by context
efficiency:

- Prior "open" occupancy, from DNA alone before any calls, rises from **0.000** in the least
  efficient decile to **0.262** in the most efficient.
- Prior TF_short_strong falls from 0.105 to 0.024 across the same deciles.

Biology would predict the opposite for A-tracts, since poly(dA:dT) tends to exclude nucleosomes
(inferred, not tested here). So short-TF calls, and anything the prior learned about TF motifs,
are contaminated by A/T-rich sequence. The proposed fix is context-offset emissions with a fixed
external table. It is not implemented yet: it changes core HMM code, and the design is to be
agreed first.

# Key Points

- **Learned emissions** (wandb run summary; the value is P(call = protected | group)):

  | group | P(protected) | implied P(m6A) |
  | --- | --- | --- |
  | open | 0.091 | 0.91 |
  | TF_short_weak | 0.550, exactly the `closed_band` floor | 0.45 |
  | TF_short_strong | 0.922 | 0.08 |
  | TF_long_weak | 0.841 | 0.16 |
  | TF_long_strong / other_long / histone | 0.986 / 0.990 / 0.994 | ≤ 0.014 |

- **Interpretation (inferred):**
  - "Open" at 0.91 m6A is above any context's rate. The adjusted table tops out near 0.68 even
    in the most accessible bin, so this state is "runs of m6A at efficient contexts".
  - TF_short_weak, pinned at 0.45 m6A, is ordinary open DNA: a second open state.
  - TF_short_strong (0.08 m6A) covers open DNA at low-efficiency contexts.
- **Occupancy by efficiency decile** (`hmm_bias_check_60x8.txt`):

  | decile | median e | observed m6A | post. open | post. TF_short_strong | post. TF_short_weak | prior open | prior TF_short_strong |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | 0 | 0.022 | 0.005 | 0.000 | 0.113 | 0.015 | 0.000 | 0.105 |
  | 2 | 0.165 | 0.042 | 0.014 | 0.117 | 0.068 | 0.014 | 0.113 |
  | 4 | 0.342 | 0.107 | 0.058 | 0.088 | 0.087 | 0.059 | 0.086 |
  | 6 | 0.486 | 0.193 | 0.137 | 0.056 | 0.093 | 0.137 | 0.055 |
  | 9 | 0.614 | 0.285 | 0.257 | 0.024 | 0.055 | 0.262 | 0.024 |

  Histone posterior falls from 0.61 to 0.39 across the deciles; other_long stays about 0.21–0.24.
- **Posterior averages ≈ prior averages per decile.** This is expected in expectation for a
  calibrated model (law of total expectation), so it does **not** show that the prior dominates.
  The finding is that the DNA-only prior carries the gradient at all. The CNN prior has 3 conv
  layers of kernel 9, a receptive field of about 25 bp, which is enough to learn 7-mer bias.

# Details

## Diagnostic (`hmm_bias_check.py`)

- **Windows:** 60 random windows on the validation chromosomes chr8 and chr20 of `active.h5`
  (seed 0), 8 random molecules each, full 2114 bp. DNA comes from `smf_hmm.data.load_dna`.
- **Inference:** for each molecule, `log_posterior_state_marginals` on the densified calls
  (h5 1/0/−1 as the model's protected/m6A/missing). The prior is
  `log_prior_state_marginals(L, dna)`. States are summed to groups via the spec's
  `state_to_group`.
- **Efficiency:** `e = sigmoid(f[ctx] + g[49])` from `adj_table_3_3.npz`. That is the
  accessibility-adjusted 3+3 context offset at the most accessible active bin, with the context
  A-centred and reverse-complemented at T sites. Sites are binned into deciles of `e` over the
  called sites.
- **A 10-window × 4-molecule pilot** (40K sites) gave the same gradients.

**Running it.** The GPU in allocation 50179010 was held in `Exclusive_Process` mode by the
user's concurrent SMFNet training. The checkpoint stores a `cuda:0` sharding, so orbax would
not restore it on CPU. The script monkeypatches
`orbax.checkpoint.type_handlers._deserialize_sharding_from_json_string` to return a CPU
`SingleDeviceSharding`; this is a diagnostic-only workaround, and `smf_hmm` code was not
changed. It ran pinned to two of the allocation's cores to leave the training alone, in about
5 minutes.

```bash
cd $SCRATCH/kernel_denoise_proto
JAX_PLATFORMS=cpu XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2" \
  taskset -c 42,49 /home/groups/btrippe/diamant/miniforge/envs/hmm/bin/python hmm_bias_check.py 60 8
```

The arguments are the number of windows and molecules per window. `taskset` cores must be
among the allocation's (`taskset -cp $$`).

## Training command (from `train.sbatch`, Slurm job 49057551)

The important settings:

- `transition=fixed_blocks`, with `prior.kernel_size=9`, `prior.normalization=global`.
- `emission.open_protected_probability=0.3`, `open_band=[0.05,0.55]`,
  `closed_band=[0.55,0.995]`.
- Ordered groups: TF_short_weak 0.75 → TF_short_strong 0.95, TF_long_weak 0.75 →
  TF_long_strong 0.95, other_long 0.9, histone 0.97.
- `objective.prior_occupancy_kl_weight=3.0`, with targets open 0.31, each TF group 0.015,
  other_long 0.08 and histone 0.55.
- Chromosomes: train on the rest, val chr8 + chr20, test chr1 + chr3 + chr6; 100 epochs,
  lr 1e-3. The final loss was 61.96.

## Proposed fix (not implemented)

1. **Context-offset emissions.** Use `logit P(m6A | state, site) = logit_state + f[ctx(site)]`,
   which is the same additive-logit form as the adjusted context fit, so its `f` table plugs in
   directly. Compute the per-position offsets from the DNA the model already receives
   (`PriorConditioning.dna`): the ±3 context, reverse-complemented at T sites.
2. **Fix `f` from the external table and learn only the state logits.** If the offsets are
   learned jointly, the DNA-conditioned prior and the emission offsets can both explain low m6A
   at A-tracts, and the model cannot tell them apart.
3. **Refit the table on HMM-train chromosomes only.** `adj_table_3_3.npz` was fit on everything
   except chr8/9/10, so it saw the HMM test chromosomes (chr1/3/6) and validation chr20.
4. **Code impact:** `hmm.py` applies a constant `(K, C)` `log_emission[:, x_safe].T` at 6 call
   sites. These need per-position emission terms, and the new emission class must be registered
   in `checkpointing.py`. The emission bands and occupancy targets need re-thinking around
   context-centred logits.
5. **Success criteria:** the prior-by-efficiency gradient above flattens, TF_short_weak comes
   off its bound, open's P(m6A) drops toward the mean-efficiency open rate, and held-out site
   log-likelihood on validation improves.

# Related Notes

- [m6A sequence-context bias in GM12878 Fiber-seq](fiberseq-m6a-sequence-context-bias.md):
  Source of the efficiency table and the about 100-fold context spread. Its adjusted model is
  the emission form proposed here.
- [Posterior predictive checks of the fixed-blocks SMF HMM](../2026-10-22/posterior-predictive-checks.md):
  Predicted exactly this failure for Fiber-seq ("enzyme bias is absorbed into states ... the
  DNA-conditioned CNN prior can learn those k-mers as binding") and planned per-k-mer emission
  offsets. This note confirms it on a trained fiber model. It also reported emissions on their
  bounds for the mESC SMF run, as TF_short_weak is here.
- [Cross-read kernel denoising of undermethylated Fiber-seq](fiberseq-cross-read-kernel-denoising.md):
  The same efficiency table, applied to own-read smoothing; it reduced sequence-tracking but
  barely changed call-level AUROC.

# Open Questions

- After the fix, how much short-TF occupancy remains at ChIP-supported motifs (CTCF, EBF1, SPI1,
  PAX5, YY1) versus A/T-rich non-motif DNA?
- Is a logit-additive offset shared by all states right, or does the protected-state background
  need its own context dependence (open vs protected offsets, as the PPC note planned)?
- How much of the depletion in A/T runs is the PacBio caller or the `ML >= 230` threshold rather
  than Hia5? This matters for whether the offset belongs in the emission or in a call-level
  model.

# Sources

- Run `/scratch/users/diamant/2026-10-22_fiber_fixed_blocks_global_kl3`: `train.sbatch`,
  `slurm-49057551.log` (wandb summary), `checkpoints/best`, `checkpoints/best.model_spec.json`.
- `/scratch/users/diamant/kernel_denoise_proto/hmm_bias_check.py`, `hmm_bias_check_60x8.txt`.
- smf_hmm `evaluation/datasets/smf/context_bias/adj_table_3_3.npz` (commit `d20f916`);
  `smf_hmm/emission_modules.py`, `smf_hmm/hmm.py`, `smf_hmm/checkpointing.py`.
- User discussion, 2026-10-30.
