---
title: Hypothesis "the prior's extra TF bumps are other TFs" — mostly refuted; prior TF peaks and prior-added TF calls are not enriched at footprinting-motif matches at strict or relaxed thresholds, apart from individual cases
date: 2026-11-16
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - /scratch/users/diamant/fiber_eval_work/footprint_motif_hits.py (48 TFs the screen found to footprint, tf_screen/genome_screen_ranked.tsv reliable & footprint)
  - /scratch/users/diamant/fiber_eval_work/prior_tf_bumps/{bump_sources.py, bias_vs_bumps.py} -> bump_sources.tsv, relaxed_motifs.tsv
  - /scratch/users/diamant/fiber_eval_work/tf_mechanism/segments.py -> segments/tf_segments.tsv.gz, summary.txt
  - Annotated figures /scratch/users/diamant/gm12878_validation/gelu_vs_collapsed_motifs/ (gelu_motif_panels.py)
tags: [smf-hmm, fiberseq, motifs, prior, tf-calls]
---

# Summary

**Hypothesis** (user, 2026-11-15): the GELU AG prior's TF bumps away from benchmark motifs, and the extra
posterior TF calls, might be other TFs rather than spurious.

**Test.** Scan every val window for matches of the 48 motifs that the GM12878 Fiber-seq TF screen found to
footprint, and of all 100 PWMs at relaxed thresholds. Compare the share of prior TF peaks and of
off-benchmark TF calls on the motif mask with the mask's share of bases.

**Result: no enrichment at any threshold.**
- **Prior TF peaks** land on motifs at chance: 13.5% vs 13.0% area (strict), 52% vs 52% and 74% vs 72%
  (relaxed). Prior TF *mass* is mildly enriched at strict matches (17.2% vs 13.0%).
- **Prior-added TF calls avoid motifs.**
  - Off-benchmark TF calls of GELU ×1 fall on strict matches *below* chance (10.0% vs 13.0%), and
    below the constant-prior run (12.8%).
  - Among confident TF segments, the share on a strict footprinting match is 12.6% (adjacent to a
    nucleosome) and 17.2% (isolated) for ×1, against 20.0% / 27.8% for ×0.
  - So the TFs the prior adds are less motif-associated than the ones the reads call on their own.
- **Individual exceptions.** Some bumps do sit on matches, e.g. an ARNT/USF2 E-box at about 120-135 bp
  in the NRF1 chr20:23,361,937 window. Reads with no sequence information call TF there too, so it is
  protected in the data.

**Conclusion.** Other known TFs explain a minority of bumps. The bulk tracks poorly labelled sequence and
nucleosome edges instead (companion notes).

**Caveat.** Only strong PWM matches of 48 (or 100) TFs are covered. TFs outside the set, weak sites, or
footprints wider than the motif (e.g. CTCF flanks) would be missed.

# Key Points

| mask | share of bases | prior TF peaks on mask (GELU seed 0) | GELU ×1 off-benchmark TF calls | collapsed ×1 |
|---|---|---|---|---|
| 48 footprinting TFs, FPR 1e-5, ±5 bp | 13.0% | 13.5% | 10.0% | 12.8% |
| 48 footprinting TFs, FPR 1e-4, ±15 bp | 52.1% | 52.2% | 41.2% | 48.2% |
| all 100 PWMs, FPR 1e-4, ±15 bp | 71.8% | 73.6% | 66.4% | 67.9% |

- **Prior TF mass on the strict mask** (expected by area: 13.0%):

  | prior | GELU seeds 0/1 | ReLU base | ReLU spacer + TF/3 |
  |---|---|---|---|
  | share | 17.2% | 17.7% | 18.7% |

- **Scan rate.** Strong matches average 3.9 merged hits per 512 bp window. Chance alone would give about
  0.5, so the scan finds real sequence matches.

# Details

- **Scanner.** Bio.motifs log-odds PSSMs with background A/T 0.295, C/G 0.205; thresholds at FPR 1e-5
  (strict) or 1e-4 (relaxed); both strands; overlapping matches merged.
- **Peaks.** Prior TF peaks are local maxima of prior P(TF) > 0.05, at least 20 bp apart.
- **Off-benchmark TF calls.** Marginal P(TF) > 0.5 more than 20 bp from the window's benchmark motif.

```bash
python /scratch/users/diamant/fiber_eval_work/prior_tf_bumps/bump_sources.py
python /scratch/users/diamant/fiber_eval_work/prior_tf_bumps/bias_vs_bumps.py
python /scratch/users/diamant/fiber_eval_work/gelu_motif_panels.py   # annotated per-read figures
```

# Related Notes

- [TF nub mechanism overview](fiber-tf-nub-mechanism-overview.md): All hypotheses together.
- [Bias-correction hypothesis](fiber-tf-nubs-bias-correction.md): Where the bumps do sit (poorly
  labelled sequence).
- [Split nucleosomes and edges](fiber-tf-nubs-split-nucleosomes-and-edges.md): The other location the
  bumps track.

# Open Questions

- Do de novo motifs enriched at prior TF peaks (beyond the 100 PWMs) exist? Not tested.

# Sources

- Scripts and outputs in the frontmatter; motif source JASPAR PWMs listed in tf_screen/motifs.json.
