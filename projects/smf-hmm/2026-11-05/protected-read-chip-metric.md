---
title: A protected-read ChIP metric — model differences vanish on reads with protection at the motif
date: 2026-11-05
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /scratch/users/diamant/mesc_validation/covariation_windows/{chr8,chr16}-512-valid90-reads10/analysis_2026-11-05 (individual_correlations.tsv, unique_motif_summaries.tsv.gz)
  - chr8 .../analysis_2026-11-04_map_targets/oct4_read_contrasts.pdf (+ oct4_read_contrasts_png/)
  - smf_hmm repo, branch nd/experiments, commits a381e5f (protected metric), 78ce9fc (read-contrast plot)
  - User discussion, 2026-11-04 and 2026-11-05
tags: [smf-hmm, tf-chip, benchmark, metrics, footprint, protected-reads, oct4, ctcf, prior-occupancy-kl, weak-tf]
---

# Summary

The user proposed scoring ChIP agreement only on reads with at least one protected call inside
the motif. Such reads carry protection evidence, so the score measures how well a model tells
a TF footprint from other protection (mostly nucleosomes). It ignores how much prior TF mass a
model puts on unprotected reads, i.e. the weak-TF floor of the
[occupancy-target note](../2026-11-04/controlled-chip-metrics-and-context-offsets.md).

On these reads **all HMM models agree almost exactly**: `neural_global`, its rerun, the
MAP-target model, `single_tf` and the merged-CGC single-TF model. So the KL target changed only
how unprotected reads are scored.
- **CTCF:** a strong signal (0.68–0.69 raw, 0.65–0.66 partial given accessibility), well above a
  read-only footprint score (0.37).
- **OCT4:** a modest signal beyond accessibility (0.13–0.14 on chr8, 0.22–0.25 on chr16).
- **ESRRB and KLF4:** about 0 beyond accessibility.

The metric is now part of `analyze_covariation` as the "protected" variant, with two read-only
reference scores.

# Key Points

- **Definition** (`analyze_covariation.py`, commit `a381e5f`):
  - `<model>_protected_mean`: mean motif P(TF) over reads with a protected call inside the
    motif.
  - Correlations use only sites with at least `MIN_PROTECTED_READS = 3` such reads (column
    `n_sites_protected`).
  - Output columns: `spearman_rho_protected` and `partial_rho_protected_accessibility`.
  - Read-only references:
    - `ref_protected_fraction`: protected reads / motif-covered reads;
    - `ref_footprint_fraction`: among protected reads, those whose calls on each side 5–30 bp
      beyond the motif are all accessible, with at least one call per side.
- **Results** (chr8 / chr16; raw Spearman with pooled ChIP-nexus, then partial given
  accessibility). Sites: CTCF 3431 / 1579, ESRRB 983 / 520, KLF4 4452 / 2134, OCT4 507 / 350.

| | CTCF | ESRRB | KLF4 | OCT4 |
|---|---|---|---|---|
| **raw, protected reads** | | | | |
| protected fraction (reads only) | 0.27 / 0.27 | −0.24 / −0.21 | −0.46 / −0.39 | −0.35 / −0.40 |
| footprint fraction (reads only) | 0.43 / 0.43 | 0.14 / 0.10 | 0.15 / 0.13 | 0.19 / 0.22 |
| motif score | 0.56 / 0.54 | 0.18 / 0.22 | 0.25 / 0.24 | 0.31 / 0.44 |
| phase_type | 0.65 / 0.63 | 0.14 / 0.15 | 0.27 / 0.20 | 0.23 / 0.29 |
| neural_global_rerun (original KL target) | 0.68 / 0.69 | 0.13 / 0.11 | 0.25 / 0.19 | 0.29 / 0.37 |
| neural_global_map_targets | 0.68 / 0.68 | 0.12 / 0.11 | 0.24 / 0.18 | 0.29 / 0.35 |
| single_tf | 0.68 / 0.68 | 0.14 / 0.12 | 0.23 / 0.18 | 0.29 / 0.37 |
| neural_global prior only | 0.67 / 0.68 | 0.07 / 0.10 | 0.13 / 0.10 | 0.39 / 0.45 |
| **partial given accessibility, protected reads** | | | | |
| protected fraction (reads only) | 0.57 / 0.58 | −0.03 / 0.03 | −0.03 / −0.07 | 0.02 / 0.04 |
| footprint fraction (reads only) | 0.37 / 0.38 | 0.08 / 0.05 | −0.01 / 0.02 | 0.09 / 0.13 |
| motif score | 0.55 / 0.54 | 0.19 / 0.21 | 0.18 / 0.19 | 0.32 / 0.44 |
| phase_type | 0.61 / 0.60 | 0.04 / 0.04 | 0.01 / 0.02 | 0.06 / 0.14 |
| neural_global_rerun | 0.65 / 0.66 | 0.04 / 0.01 | −0.01 / 0.01 | 0.13 / 0.24 |
| neural_global_map_targets | 0.65 / 0.66 | 0.03 / 0.01 | −0.00 / 0.01 | 0.14 / 0.24 |
| single_tf | 0.65 / 0.65 | 0.05 / 0.03 | −0.01 / 0.01 | 0.14 / 0.25 |
| neural_global prior only | 0.65 / 0.66 | 0.00 / 0.03 | −0.01 / 0.01 | 0.29 / 0.35 |

- **For OCT4, sequence still beats reads** on these sites: the prior alone (partial 0.29–0.35)
  and the motif score (0.32–0.44) exceed every posterior (0.13–0.25).
- **Bound CTCF sites are more often protected** beyond what accessibility explains (protected
  fraction partial 0.57). The HMMs add a further footprint-versus-nucleosome distinction.

# Details

**Why the models converge** (inferred). On a read with a protected motif call, open is
heavily disfavoured, so P(TF) is decided between TF and nucleosome or `other_long` states. All
these models share the fixed-block transitions and protected emission levels near 0.99. The
weak-TF floor only matters on reads where open competes, i.e. unprotected ones. This is
consistent with the read-category contributions in the occupancy-target note: protected-read
contributions were unchanged between the two KL targets.

**Read-contrast plots** (`evaluation/tf_chip/plot_read_contrasts.py`, commit `78ce9fc`):
- For 8 bound chr8 OCT4 sites, the PDF shows the MAP-target model's top footprint read and the
  control's top read whose motif calls are all accessible, with both models on each page.
- **Footprint reads** get P(TF) of 0.78–0.90 from both models.
- **Open reads** get, e.g., 0.38 from the control (mostly `TF_long_weak`) and 0.00 from the
  MAP-target model.
- The control also places weak-TF bumps of about 0.1–0.3 along most open stretches, not only at
  motifs.
- Caveat: the "open" criterion only checks calls inside the motif. At `chr8:90670968` the chosen
  read has a protected call just outside the motif edge, and the MAP-target model calls a short
  strong-TF footprint there (0.25).

**Limits:**
- Requiring 3 protected reads drops sites, e.g. ESRRB loses about 40%.
- Prior-only rows are identical across reads, so their "protected" value is the motif prior on
  the protected-site subset.
- `ref_accessibility` partial values are NaN or 0 by construction.

# Commands

The analysis ran in Slurm job 51132065 (sh04-06n09); `analyze_covariation` runs on CPU.

```bash
cd /home/users/diamant/repos/smf_hmm
M=/scratch/users/diamant/mesc_validation/covariation_windows
bash $M/analysis_2026-11-05.sh     # analyze_covariation with the protected variant, chr8 and chr16
PY=/home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python
D=$M/chr8-512-valid90-reads10
$PY -m evaluation.tf_chip.plot_read_contrasts $D/dataset.h5 \
  --reference map_targets=$D/posteriors/neural_global_map_targets.h5 \
  --contrast kl_target_control=$D/posteriors/neural_global_rerun.h5 \
  --tf oct4 --output $D/analysis_2026-11-04_map_targets/oct4_read_contrasts.pdf
```

A one-off snippet computed the first version of these numbers; the saved metric reproduces
them.

# Related Notes

- [Controlled ChIP metrics and occupancy targets](../2026-11-04/controlled-chip-metrics-and-context-offsets.md): The weak-TF floor and the MAP-target model compared here; the floor is what this metric ignores.
- [ChIP evaluation: single-TF vs multi-TF](../2026-11-03/chip-eval-single-tf-vs-multi-tf.md): The OCT4 footprint reads and the small binding-specific OCT4 footprint that the residual OCT4 signal here fits.

# Open Questions

- For OCT4, why do sequence (the prior, the motif score) predict ChIP better than read evidence
  on protected-read sites? Is it because ChIP-nexus intensity tracks motif strength more than
  occupancy per molecule?
- Do other ChIP-nexus TFs with strong footprints behave like CTCF? MAX::MYC, NF-Y and ZFP281
  were being evaluated overnight on 2026-11-05.

# Sources

- Analysis: `/scratch/users/diamant/mesc_validation/covariation_windows/{chr8,chr16}-512-valid90-reads10/analysis_2026-11-05/`, script `.../covariation_windows/analysis_2026-11-05.sh`
- Plots: `.../chr8-512-valid90-reads10/analysis_2026-11-04_map_targets/oct4_read_contrasts.pdf`
- smf_hmm `nd/experiments`: `a381e5f`, `78ce9fc`
