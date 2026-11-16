---
title: Hypothesis "deterministic spacer-TF-spacer transitions make TF calls favourable" — refuted as stated (a TF excursion costs about 8.5 nats of prior); instead the learned per-microstate prior potentials penalize TF position 9, turning the 8 bp minimum into the preferred TF length
date: 2026-11-16
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /scratch/users/diamant/fiber_eval_work/tf_mechanism/out/structure_gelu_s0.json (gpu_diagnostics.py part structure)
  - Prior output-layer biases from checkpoints: gelu_s0 (2026-11-15_fiber_gelu_head_sptf_nochain_B_seed0), collapsed (2026-11-14_fiber_ablation_sptf_nochain_seed2), ag_sptf (2026-11-12_..._ol001_spacer10_tf0067_kl10)
  - smf_hmm transition_modules.FixedBlockTransition (new tied_blocks option, uncommitted)
tags: [smf-hmm, fiberseq, transitions, spacer, prior, tf-duration]
---

# Summary

**Hypothesis** (user): with deterministic spacer chains (open → 10 bp entry spacer → TF → 10 bp exit
spacer → open), a spacer-TF-spacer excursion might be favourable.

**Cost accounting refutes it as stated.** Relative to staying open, from the plain open state:
- entering the spacer costs −6.63 nats;
- an 8 bp duration costs about −2.4;
- the open self-transitions skipped are worth about +0.5.

Net, about −8.5 nats of prior. An 8 bp TF with about 3 protected calls gains only about 2.3 nats of
emission over open (0.76 each). So TF nubs need a position-specific prior bump to appear; consistently,
×0 makes 7× fewer adjacent nubs (split-nucleosome note).

**The real transition-side factor is duration shaping by the prior.**
- **Per-state potentials.** The prior assigns a potential to each of the 40 TF microstates (positions
  along the TF chain), so it can reshape the TF duration distribution.
- **What it learned.** The GELU model has learned a −1.31 nat penalty at TF position 9 (positions 1-8:
  +0.13 each; then −0.51, −0.27, …). Continuing past the 8 bp minimum is expensive, so the
  baseline-implied duration distribution moves from 8.9% at exactly 8 bp (transitions alone) to 39%.
- **Comparison.** The collapsed run learned the same shape (+1.21 at position 8, −0.60 at 9; 38% at 8 bp).
  The older ReLU spacer + TF/3 run, with chaining allowed, did not (11%; its mean duration rises to
  24 bp).
- **Effect.** This is why the prior-added TF calls are overwhelmingly exactly 8 bp.

**Fix under test.** A new `tied_blocks` option of `FixedBlockTransitionConfig` gives every microstate of
the named blocks the potential of the block's first state, as spacers already share the open state's.
The prior can then say "TF here" but not "an 8 bp TF here". It is being trained with TF and other_long
tied (fix trials, overview note).

# Key Points

- **Prior output bias per TF microstate** (baseline potential before the sequence-dependent part):

  | model | positions 1-8 | 9 | 10 | 11 | 12 | mean 9-40 |
  |---|---|---|---|---|---|---|
  | GELU seed 0 | +0.13 | −1.31 | −0.51 | −0.27 | −0.15 | −0.08 |
  | collapsed (ReLU, seed 2) | +0.10 to +1.21 | −0.60 | −0.83 | −0.45 | −0.17 | −0.09 |
  | ReLU spacer + TF/3 (chaining) | +0.29 | −0.60 | −0.28 | −0.13 | −0.09 | +0.03 |

- **TF duration distribution implied by transitions × baseline potentials** (relative to open's potential):

  | | P(8 bp) | mean | P(≤ 10 bp) |
  |---|---|---|---|
  | transitions alone | 0.089 | 17.7 bp | 0.25 |
  | GELU seed 0 | 0.393 | 14.2 bp | 0.55 |
  | collapsed | 0.378 | 11.8 bp | 0.65 |
  | ReLU spacer + TF/3 | 0.108 | 24.1 bp | 0.20 |

- **Observed confident off-motif TF segments** (GELU ×1): 53% exactly 8 bp, against about 3% at ×0
  (split-nucleosome note figure).

# Details

- **Potentials.** Each potential is the output-layer bias minus the open state's bias. Spacer states
  share open's potential; the sequence-dependent part W·h varies the potentials along the sequence.
  The table shows only the baseline.
- **Implied pmf.** p(d) ∝ p_transition(d) · exp(Σ_{j ≤ d} (b_TF,j − b_open)).
- **Tied blocks (code).**
  - `FixedBlockTransitionConfig.tied_blocks`; `prior_potential_index` maps each tied block's states to
    its first state.
  - Transition probabilities and stationary occupancies are unchanged.
  - Tests: tests/test_transition_modules.py (3 new); the full suite passes (483).

```bash
# training (tied TF and other_long; GELU head; design B; seed 0)
sbatch /scratch/users/diamant/2026-11-16_fiber_gelu_B_tied_alphagenome_bf16/train.sbatch   # 52509099
```

# Related Notes

- [TF nub mechanism overview](fiber-tf-nub-mechanism-overview.md): All hypotheses and the fix-trial
  results.
- [Split nucleosomes and edges](fiber-tf-nubs-split-nucleosomes-and-edges.md): Where the 8 bp nubs sit.
- [TF chaining note](../2026-11-13/fiber-no-other-long-tf-chaining.md): The no-chaining change; the
  duration shaping appears in no-chaining models.

# Open Questions

- Why does no-chaining training learn the 9th-position penalty when the chaining design did not?
  Possibly the remnants beside nucleosomes are short, and without chains short TFs are the only fit
  (inferred).
- Should histone microstates be tied too? Their per-state potentials can only shift placement since
  the length is fixed.

# Sources

- Checkpoints and outputs in the frontmatter.
