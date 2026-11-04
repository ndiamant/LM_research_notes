---
title: ChIP evaluation of single-TF models — multi-TF wins on OCT4 and KLF4 through its weak TF states and its prior
date: 2026-11-03
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /scratch/users/diamant/mesc_validation/covariation_windows/{chr8,chr16}-512-valid90-reads10/analysis_2026-11-03 (individual_correlations.tsv, pair_correlations.tsv, unique_motif_summaries.tsv.gz, results.pdf)
  - New posteriors in .../posteriors/{single_tf_147_merged_cgc,single_tf_tapered,fixed_blocks_tapered,single_tf}.h5; inference script /scratch/users/diamant/mesc_validation/covariation_windows/infer_2026-11-03.sh
  - chr8 .../analysis_2026-11-03/oct4_disagreement_sites.pdf (+ oct4_sites_png/, oct4_disagreement_site_keys.txt)
  - smf_hmm repo, branch nd/experiments, commit 7cd90ea (posterior inference fix); plot_posterior_examples.py --site-key and summary pages (uncommitted at time of writing)
  - User discussion, 2026-11-03
tags: [smf-hmm, tf-chip, benchmark, oct4, klf4, ctcf, esrrb, single-tf, multi-tf, weak-tf, cnn-prior, emissions, posterior-inference, bug]
---

# Summary

We ran the motif-anchored ChIP covariation benchmark on four fixed-block models, on the chr8
and chr16 covariation windows. Results were compared with the existing geometric, phase-type,
`neural_global` and HiddenFoot posteriors.

**Main result.** The three single-TF models, which have a sequence-conditional TF emission
head, tie the multi-TF fixed-block models on CTCF and ESRRB. They are clearly worse on **OCT4**
(Spearman with ChIP 0.30–0.41 vs 0.50–0.57) and **KLF4** (0.08–0.13 vs 0.19–0.24). The new
merged-CGC / 147 bp run is no better than plain `single_tf`.

**Where the OCT4 gap comes from:**
- **Bound sites only.** At low-ChIP motifs every model gives P(TF) ≈ 0.05. At top-quartile ChIP
  motifs, multi-TF gives 0.14–0.18 and single-TF stays at about 0.08.
- **The weak TF groups.** The multi-TF gain is almost all in `TF_short_weak` plus
  `TF_long_weak`. Both learned P(x = 1) = 0.20, which is the closed band's floor; open is 0.01
  and the strong groups 0.99.
- **The sequence prior.** Before any reads are seen, prior-only P(TF) at OCT4 motifs correlates
  with ChIP at 0.39–0.43 for multi-TF, vs 0.20–0.27 for single-TF. The reads then raise every
  model by about the same amount.

**Bug fixed along the way** (commit `7cd90ea`). `smf_hmm.cli.infer_posterior` built posteriors
from prior terms only, so the sequence-conditional emission head was never applied in ChIP
inference. Training was unaffected. Earlier ChIP numbers for single-TF models would have been
computed without the head; `neural_global` has no head and is unaffected.

# Key Points

- **Individual sites:** Spearman between motif-averaged P(TF) and pooled ChIP-nexus, chr8 /
  chr16.

| model | CTCF | ESRRB | KLF4 | OCT4 |
|---|---|---|---|---|
| sites | 3780 / 1754 | 1674 / 859 | 6655 / 3329 | 810 / 520 |
| geometric | 0.60 / 0.63 | 0.25 / 0.30 | 0.11 / 0.08 | 0.26 / 0.30 |
| phase_type | 0.64 / 0.64 | 0.25 / 0.30 | 0.24 / 0.20 | 0.33 / 0.38 |
| hiddenfoot | 0.30 / 0.25 | 0.08 / 0.05 | −0.18 / −0.17 | −0.05 / −0.06 |
| neural_global (multi-TF, `2026-09-29_fixed_blocks_global_kl3`) | **0.67 / 0.69** | 0.25 / 0.30 | **0.24 / 0.21** | **0.50 / 0.57** |
| fixed_blocks_tapered (multi-TF, tapered histone) | 0.66 / 0.68 | **0.28 / 0.32** | 0.23 / 0.19 | **0.50 / 0.56** |
| single_tf (`2026-11-02_fixed_blocks_single_tf_global_kl3`) | 0.64 / 0.67 | 0.27 / 0.30 | 0.13 / 0.11 | 0.34 / 0.41 |
| single_tf_tapered (`2026-09-29_single_tf_tapered_histone_global_kl3`) | 0.65 / 0.66 | 0.21 / 0.25 | 0.12 / 0.08 | 0.30 / 0.33 |
| single_tf_147_merged_cgc (`2026-11-03_...`) | 0.65 / 0.68 | 0.24 / 0.30 | 0.11 / 0.08 | 0.33 / 0.38 |

- **Histone taper:** no consistent effect. It's neutral between the two multi-TF models and
  slightly negative among the single-TF models.
- **Motif pairs** (`nonoverlapping_chip` scope, all distances):
  - Every HMM model stays near 0 on same-TF pairs (−0.1 to 0.2), where ChIP is 0.6–0.75. None
    of the new models changes this.
  - HiddenFoot scores 0.3–0.55 on every pair, including CTCF–ESRRB, where ChIP is −0.18 (chr8)
    and −0.29 (chr16). This suggests a window-wide shared signal such as accessibility, not TF
    covariation (inferred).
- **OCT4, mean motif P(TF) by ChIP quartile** (chr8; chr16 is similar):

| ChIP quartile | neural_global | fixed_blocks_tapered | single_tf | single_tf_147_merged_cgc |
|---|---|---|---|---|
| Q1 (low) | 0.059 | 0.057 | 0.052 | 0.048 |
| Q4 (high) | 0.147 | 0.141 | 0.083 | 0.081 |

- **`neural_global` OCT4 P(TF) by group**, Q1 → Q4 (chr8; chr16 in brackets):
  - TF_short_weak: 0.013 → 0.045 (0.014 → 0.063)
  - TF_long_weak: 0.020 → 0.063 (0.018 → 0.073)
  - TF_short_strong: 0.012 → 0.017 (0.012 → 0.017)
  - TF_long_strong: 0.014 → 0.022 (0.013 → 0.024)
- **Prior-only vs posterior** Spearman with ChIP at OCT4 motifs, chr8 / chr16:

| model | prior only | posterior |
|---|---|---|
| neural_global | 0.39 / 0.43 | 0.50 / 0.57 |
| single_tf | 0.26 / 0.27 | 0.34 / 0.41 |
| single_tf_147_merged_cgc | 0.22 / 0.20 | 0.33 / 0.38 |

  Correlation with motif score shows the same pattern: 0.33–0.37 for multi-TF vs 0.17–0.25 for
  single-TF.

# Details

**Bug in posterior inference** (fixed in `7cd90ea`).
- `smf_hmm/posterior_inference.py` used `model.prior_terms` and
  `log_posterior_state_marginals_with_prior_terms(x, initial, bias)`. With no `log_emission`,
  the HMM falls back to `emission_module.log_emission_matrix()`, which is the zero-offset
  emissions.
- It now uses `model.region_terms`, as `smf_hmm/train.py` does.
- A new `--assay-type` option (`InferenceConfig.assay_type`) merges CGC calls for models trained
  on `cg_gc_accessibility_merged_cgc`.
- Two new tests:
  - inference must match `log_posterior_state_marginals` on a model with a non-zero head, and
    must differ from the zero-offset posterior;
  - merged-assay inference must equal inference on hand-merged observations.
- Still unfixed: the same omission is in `smf_hmm/state_call_inference.py` (the older
  state-call evaluation) and in the user's untracked `smf_hmm/posterior_checks/heldout.py` and
  `replicates.py`.

**Interpretation (inferred, not tested):**
- **Weak TF states are cheap at motifs.** A weak TF state at P(x = 1) = 0.2 costs only about
  log(0.8/0.99) ≈ −0.21 nats per accessible call relative to open. The multi-TF prior can
  therefore learn "OCT4 motif → weak TF" without strong read evidence.
- **Single-TF can't do this.** Its one TF state emits P(x = 1) ≈ 0.80 at zero offset, so each
  accessible call costs about −1.6 nats. Its prior learns less from motifs, and its emission
  head can't lower protection that far: across positions its range is 0.55–0.99, against the
  0.20 of the multi-TF weak groups.
- **OCT4 looks like a weak footprinter in this data.** The data's motif-minus-flank x = 1 rate
  at bound OCT4 motifs is +0.05, against +0.23 for CTCF (from the [merged-CGC run
  note](merged-cgc-147bp-histone-run.md)).
- **Caveat on what the benchmark rewards.** Part of the multi-TF advantage is the prior acting
  as a motif scanner. The benchmark rewards that, but it isn't footprint detection. In the
  example plot below, `neural_global` puts about 0.3 weak-TF probability on a bound OCT4 motif
  whose calls look mostly accessible.

**Example plots.**
- `evaluation/tf_chip/plot_posterior_examples.py` was extended (uncommitted when this was
  written):
  - `--site-key` plots chosen sites, skipping those with fewer than 10 motif-covered reads;
  - every site now gets a summary page first: a heatmap of all reads' calls sorted by the first
    model's motif P(TF), read-averaged P(TF) for all models overlaid, and per-model
    read-averaged group posteriors. The existing single-read pages follow.
- Sites chosen: the 12 chr8 OCT4 motifs in the top ChIP quartile with the largest
  `neural_global` − `single_tf_147_merged_cgc` gap in motif P(TF). Five lacked 10 motif-covered
  reads, leaving 7 sites.
- At `chr8:90670968-90670983:+` (ChIP 1489, motif score 19.0), the reads show protected calls
  over the motif. Read-averaged motif P(TF) is 0.51 for `neural_global`, almost all from its
  weak groups, vs 0.27 for `single_tf` and 0.22 for the merged-CGC model.
- At `chr8:24571331-24571346:-` (ChIP 292), the motif calls are mostly accessible.
  `neural_global` still gives 0.31, while the single-TF models give 0.03–0.05.
- The single-read pages, which show each model's highest-scoring read, look similar across
  models. The difference is in the read average.

**Proposed directions for sequence dependence** (discussed with the user; none implemented):
1. **Build on the multi-TF structure.** It already gets its sequence dependence through the
   prior choosing weak vs strong TF types. That is the "TF block types chosen by the prior"
   design proposed earlier, and it outperforms. Port the merged-CGC assay and the 147 bp
   invaded-edge histone to it.
2. **Add the emission head there, if at all.** Put it on the weak TF groups, so their
   protection can vary by motif, instead of a single TF state.
3. **Report prior-only correlations alongside posteriors**, so motif-scanning by the prior
   can't be mistaken for footprint detection.
4. **Watch the weak TF floor.** The weak groups sit at the 0.2 floor. With a closed band of 0.1
   they would move closer to open, and could absorb more prior-only calls. Keep an eye on the
   call-free TF segments noted in the merged-CGC run note.

# Commands

Everything ran in Slurm job 50835752 (akundaje partition, 1 L40S). Inference took about
6 min per chr8 model and 2 min per chr16 model. Analysis took 277 s (chr8) and 112 s (chr16),
using about 1.4 GB.

```bash
cd /home/users/diamant/repos/smf_hmm
M=/scratch/users/diamant/mesc_validation/covariation_windows
bash $M/infer_2026-11-03.sh       # 4 models x {chr8, chr16}; only the merged-CGC model gets --assay-type cg_gc_accessibility_merged_cgc
PY=/home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python
for CHR in chr8 chr16; do D=$M/$CHR-512-valid90-reads10; P=$D/posteriors
  $PY -m evaluation.tf_chip.analyze_covariation $D/dataset.h5 \
    --prediction geometric=$P/geometric.h5 --prediction phase_type=$P/phase_type.h5 \
    --prediction neural_global=$P/neural_global.h5 --prediction hiddenfoot=$P/hiddenfoot.h5 \
    --prediction fixed_blocks_tapered=$P/fixed_blocks_tapered.h5 --prediction single_tf=$P/single_tf.h5 \
    --prediction single_tf_tapered=$P/single_tf_tapered.h5 --prediction single_tf_147_merged_cgc=$P/single_tf_147_merged_cgc.h5 \
    --output-dir $D/analysis_2026-11-03
done
D=$M/chr8-512-valid90-reads10
$PY -m evaluation.tf_chip.plot_posterior_examples $D/dataset.h5 \
  --run neural_global=$D/posteriors/neural_global.h5 --run single_tf=$D/posteriors/single_tf.h5 \
  --run single_tf_147_merged_cgc=$D/posteriors/single_tf_147_merged_cgc.h5 \
  $(cat $D/analysis_2026-11-03/oct4_disagreement_site_keys.txt) \
  --output $D/analysis_2026-11-03/oct4_disagreement_sites.pdf
ml poppler && pdftoppm -r 70 -png $D/analysis_2026-11-03/oct4_disagreement_sites.pdf $D/analysis_2026-11-03/oct4_sites_png/page
```

- The quartile, group and prior-only numbers came from one-off snippets, not saved scripts:
  - the quartile and group splits were read from `unique_motif_summaries.tsv.gz`;
  - the prior-only P(TF) was computed by running each checkpoint's prior on the prepared windows'
    DNA and averaging over the motif span.
- An earlier PDF without summary pages, `oct4_disagreement_examples.pdf`, is also in the
  analysis directory.

# Related Notes

- [Merged-CGC assay and 147 bp histone run](merged-cgc-147bp-histone-run.md): The `single_tf_147_merged_cgc` model evaluated here, its call-free TF segments, and the per-TF data footprints.
- [Single-TF emission head learns a CGC artifact](../2026-11-02/single-tf-emission-head-learns-cgc-not-tf-identity.md): Why the single-TF head was expected to carry TF identity, and evidence that it carries little.
- [HiddenFoot baseline: bugs and variants](../2026-10-23/hiddenfoot-baseline-bugs-and-variants.md): Source of the `hiddenfoot.h5` posteriors and the `analysis_hiddenfoot_fixed` convention reused here.
- [Chr16 covariation, ChIP width, and window selection](../2026-09-28/chr16-covariation-chip-width-and-window-selection.md): The covariation-window benchmark design.

# Open Questions

- Does the multi-TF model with merged CGC and a 147 bp invaded-edge histone keep its OCT4/KLF4
  advantage?
- Would a sequence-conditional head on its weak TF groups add read-level evidence beyond the
  prior?
- How much of each model's ChIP correlation is prior-only motif scanning? This should be
  measured for all TFs, not only OCT4.
- Why are all HMM pair correlations near zero while ChIP pair correlations are 0.6–0.75? Is it
  the 512 bp window, or the HMM's TF-call noise?

# Sources

- Analysis outputs: `/scratch/users/diamant/mesc_validation/covariation_windows/{chr8,chr16}-512-valid90-reads10/analysis_2026-11-03/`
- Posteriors: `.../posteriors/{single_tf_147_merged_cgc,single_tf_tapered,fixed_blocks_tapered,single_tf}.h5` (new); `{geometric,phase_type,neural_global,hiddenfoot}.h5` (existing)
- Checkpoints: `/scratch/users/diamant/{2026-11-03_single_tf_147_histone_merged_cgc_global_kl3,2026-09-29_single_tf_tapered_histone_global_kl3,2026-09-29_fixed_blocks_tapered_histone_global_kl3,2026-11-02_fixed_blocks_single_tf_global_kl3,2026-09-29_fixed_blocks_global_kl3}/checkpoints`
- smf_hmm `nd/experiments`: `7cd90ea` (`smf_hmm/posterior_inference.py`, `smf_hmm/cli/infer_posterior.py`, `tests/test_posterior_inference.py`); `evaluation/tf_chip/plot_posterior_examples.py` and `tests/test_tf_chip_plot_posterior_examples.py` (uncommitted)
