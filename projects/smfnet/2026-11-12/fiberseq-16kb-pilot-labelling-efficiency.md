---
title: 16 kb Fiber-seq pilot (16,384 bp of reads from 18,432 bp of DNA); its gains come from reading the rest of the molecule, its generated reads lack half the molecule-wide (labelling) variance, two inserted CTCF elements co-vary even 12 kb apart (orientation-independent), and a design for modelling labelling efficiency and enzyme bias
date: 2026-11-12
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - SMFNet-dev branch nd_unet_figwig, 961d8ba (configurable dna_flank, min_read_fraction); runs from the detached checkout ~/repos/SMFNet_unet_figwig_runs at 961d8ba
  - Model smf_models/GM12878_fiberseq_figwig_16kb_flank1024_alphagenome_unet_attn_muon_ema_k9_dilated (epoch24-step91275.ckpt), trained in allocation 51969713 (H100, sh04-03n07)
  - Outputs in $SCRATCH/smfnet_figwig/eval16/ and $SCRATCH/tm_inspect/coaccessibility/
  - fiberseq-16kb-pilot-labelling-efficiency/scripts/ (copied from $SCRATCH)
tags: [fiberseq, gm12878, figwig, long-context, unet, co-accessibility, ctcf, labelling-efficiency, enzyme-bias, hia5, emission, evaluation, design]
---

# Summary

A 25-epoch pilot of the 2 kb run's best U-Net at long context: 16,384 bp of reads with
1,024 bp of DNA on each side (`dna_flank=1024`), reads drawn only if they cover >= 75% of
the read window. Its val loss (0.157) is lower than the 2 kb model's (0.179) mostly because
16 kb windows are mostly easier nucleosomal sites. On the same sites (the central 1 kb, same
reads and masks) it is 4.1% better than the 2 kb model, and all of that comes from reading
the rest of the molecule: with the flanking calls masked it is 0.7% worse.

Generated reads reproduce profiles and distance-dependent co-accessibility, but only about
half of the between-molecule variance (spread of read means 0.42 of observed; long-range
co-accessibility plateau 0.083 vs 0.158), at 32, 64 or 128 sampling steps alike. That
component is technical: per-molecule labelling efficiency. Adjusting for each read's
labelling level, observed and generated co-accessibility match at long range and the
generated short-range (nucleosome-scale) coupling is 10-20% weak. Specific long-range
co-accessibility does exist between CTCF sites (partial r +0.14 given labelling level vs
+0.02 for random position pairs), and the model generates it: two inserted 3x CTCF elements
co-vary on generated molecules (partial r ~0.19 at 2-8 kb, 0.15 at 12 kb, against ~0 for an
element and an empty site), far beyond the convolutions' +/-1.5 kb reach, with no
dependence on motif orientation. Profiles are near the pseudo-replicate ceiling (98-99.6%
of the gap closed) but run lengths are not (83-85%), consistent with the labelling deficit.

A per-molecule labelling efficiency b_m, estimated from each molecule's offset from the
other molecules at the same background sites, validates well: Spearman 0.89 with the raw
molecule rate, nucleosome footprints shrinking from 158 to 118 bp across its quintiles, and
almost no information about region activity (+0.05 against window accessibility measured
on independent molecules, vs +0.64 for the raw rate). With it as a covariate the CTCF-CTCF
pair coupling holds (+0.112 vs +0.005 for random pairs).

The agreed next step models labelling explicitly with that efficiency, in an emission head
P(m6A) = q * sigma(c_context + b_molecule) + (1 - q) * eps, with the same quantities given to
the network as conditioning.

# Key Points

- **Architecture: `CNNConfig.dna_flank`** (default 512 = the old middle-half crop at
  2,048 bp, so every 2,048 bp checkpoint is unchanged; rescored val loss 0.17912 = logged).
  With AlphaGenome it must be a multiple of 64 bp. A fixed flank instead of a context that
  scales with the window: the encoder's 64 bp features see about +/-600 bp, so 1,024 bp
  clears them; 18 kb of DNA is 44% less encoder work than 32 kb.
- **Reads at 16 kb:** windows with >= 32 reads covering >= 75% of the read window: 98.2% of
  active and 96.9% of background kept-window centres (median 61 reads). Under the 2 kb (H5)
  rule only 8% of drawn reads would cover a 16 kb window.
- **Training:** 32 regions x 32 reads (77 GB of 80 at peak; 16 x 32 is 40 GB), Muon lr 0.01,
  EMA 0.999, 25 epochs x 3,651 steps (matched to the 2 kb run's 1.58T read positions),
  10 h 58 min; 116,807 train / 11,923 val regions; val loss 0.15686, still falling.
- **Why its val loss is lower:** per-site entropy of the scored sites (the loss a model
  knowing only per-site rates would get) is 0.272 nats for 16 kb windows vs 0.310 for the
  2 kb model's centres, ratio 0.878, the same as the val-loss ratio 0.881.
- **Same sites, central 1 kb** (11,923 val regions; `rescore_centre.py`):

  | | 2 kb model | 16 kb model | 16 kb, calls outside centre masked |
  |---|---|---|---|
  | all | 0.1779 | 0.1707 (-4.1%) | 0.1791 (+0.7%) |
  | active | 0.1855 | 0.1785 (-3.8%) | 0.1870 (+0.8%) |
  | most accessible 10% | 0.2558 | 0.2467 (-3.5%) | 0.2581 (+0.9%) |

- **Generated reads** (1,000 val regions x 512 reads; generated subsampled to the observed
  count for curves):

  | distance | raw obs | raw gen | adjusted obs | adjusted gen |
  |---|---|---|---|---|
  | 200 bp | 0.509 | 0.430 | 0.408 | 0.372 |
  | 1 kb | 0.228 | 0.140 | 0.065 | 0.052 |
  | 2-4 kb | 0.170 | 0.089 | -0.005 | -0.005 |
  | 8-16 kb | 0.158 | 0.083 | -0.011 | -0.012 |

  Same plateau at 32 / 64 / 128 steps (0.079 / 0.079 / 0.078 at 8-16 kb). Spread of per-read
  means, generated / observed (generated reads given observed no-call patterns): 2 kb U-Net
  0.92, 2 kb CNN 0.89, 16 kb pilot 0.42; so not the read norm, and only visible at 16 kb
  where local variation averages out.
- **CRE pairs 2-15 kb apart**, partial correlation given molecule labelling level (window
  level excluding the elements), >= 50 shared molecules:

  | pairs | n | mean r | p < 0.001 | r > 0.3 |
  |---|---|---|---|---|
  | CTCF-CTCF | 493 | +0.136 | 6.7% | 16.8% |
  | CTCF-other active | 8,000 | +0.049 | 1.5% | 4.0% |
  | other active | 44,448 | +0.029 | 1.2% | 2.9% |
  | random positions | 13,288 | +0.020 | 0.85% | 2.4% |

  At >= 100 molecules: CTCF-CTCF 18% p < 0.001 (n = 55) vs random 0.85%. 480 / 52,941 real
  pairs at FDR < 0.05.
- **metrics.py** on the 1,000 val regions (no positional baselines): profile Pearson 0.756
  (ceiling 0.770, gap closed 0.981), Spearman 0.779 (0.783, 0.996); run-length Wasserstein,
  accessible 0.515 bp (ceiling 0.168, null 2.25, 0.833), protected 4.35 bp (2.28, 16.2,
  0.852). Not directly comparable with the 2 kb model's 0.96-0.97 (other windows and read
  rule), but run lengths depend strongly on labelling (light molecules merge footprints), so
  they are a natural check for the efficiency conditioning.
- **Per-molecule efficiency b_m** (19.0M molecules; median -0.07, IQR -0.48 to +0.28):
  Spearman 0.89 with the raw molecule rate; protected-stretch peak 158 / 148 / 138 / 128 /
  118 bp across b_m quintiles; window accessibility vs mean b_m Spearman +0.05 on independent
  molecules (+0.24 on the same molecules, sampling noise; raw background-site rate +0.64).
- **CRE pairs with b_m** (subsampled; partial r given window level and b_m): CTCF-CTCF
  +0.112 (4.5% p < 0.001), CTCF-other +0.044, other active +0.027, random positions +0.005.
- **Inserted CTCF pairs** (Till's suggestion; see Details): partial r 0.19-0.21 / 0.18-0.19 /
  0.17-0.19 / 0.14-0.15 at 2 / 5 / 8 / 12 kb; convergent minus divergent within
  +/-0.016 (about 1 SE) at every distance.

# Details

## Observed co-accessibility by distance, and the receptive field

The run's U-Net reaches +/-1.5 kb with its convolutions (gradient support with randomised
weights); only the 2 bottleneck attention layers (64 bp tokens) see further. Observed
co-accessibility (16 kb windows, 75% reads, 100 bp bins) with each molecule's mean removed is
0.41 at 200 bp, 0.06-0.07 at 1 kb and <= 0.01 past ~2 kb, while the raw correlation stays at
0.14-0.17 out to 16 kb (chr21 and chr8 agree; active and background alike).

![Observed co-accessibility by distance](fiberseq-16kb-pilot-labelling-efficiency/observed_coaccessibility_by_distance.png)

*Distance-dependent co-accessibility is within the conv reach; beyond it there is a flat,
molecule-wide component (blue).*

Element pairs on observed data (genome-wide, cCREs of candidates.bed): CTCF-CTCF pairs
co-vary (+0.093 with the molecule mean subtracted vs -0.01 to -0.04 for other groups), not
through haplotype (+0.091 within haplotype) and flat with distance (+0.089 / +0.092 / +0.096
at 2-5 / 5-10 / 10-15 kb). At matched accessibility much of it goes (less-open element
0.20-0.35: CTCF +0.108, mixed CTCF-other +0.058, other +0.060), leaving a CTCF-specific
+0.04-0.05; the per-pair partial-correlation test above confirms CTCF-CTCF pairs stand out.

## The pilot's generated reads

![Generated vs observed co-accessibility by distance](fiberseq-16kb-pilot-labelling-efficiency/generated_vs_observed_by_distance.png)

*Raw curves (solid) differ by the molecule-wide component; adjusted for each read's
labelling level (dashed) they match at long range and differ slightly at short range. The
pair panel (right) has too few CTCF pairs in these regions to judge.*

The adjustment regresses each bin's accessibility on the read's mean over the window (slope
per bin, so open bins can respond more to labelling) and correlates the residuals; it
also removes any biological molecule-wide variation.

![Whole-window co-accessibility, raw and adjusted](fiberseq-16kb-pilot-labelling-efficiency/full_window_raw_vs_adjusted_chr20.png)

*chr20:41.77 Mb: the observed raw map is red throughout (labelling); adjusted, observed and
all 512 generated reads show the same near-diagonal structure.*

![Example region](fiberseq-16kb-pilot-labelling-efficiency/example_region_chr8.png)

*chr8:144.09 Mb: protection profiles match; observed reads sorted by protection show a
labelling gradient the generated reads lack; the 512-read zoom shows the ~190 bp nucleosome
bands seen in the observed zoom.*

Plot conventions settled with the user: whole window at 100 bp plus a 3 kb zoom at 20 bp
around the most accessible active element; diverging blue-grey-red clipped at +/-0.6
(preferred to the Mola paper's sequential 0-0.6 scale, which hides negatives); a panel with
all generated reads; protection on the read heatmaps and profiles. Mola (Wang et al., bioRxiv
2026) bins at 20 bp over 3 kb views, shows 10x generated molecules, and scores maps by mean
row correlation against a split-half Spearman-Brown ceiling.

## In-silico insertion of two CTCF elements (Till's suggestion)

Till (colleague): does the model learn interactions beyond 5 kb? Insert two convergent CTCF
elements, three motifs each, as the strongest case. Mola ran a similar screen (two synthetic
DHSs of three homotypic motifs at offsets 0, 50, 100 bp, 1-8 kb apart).

- Design: JASPAR MA0139.1 consensus `TGGCCACCAGGGGGCGCTA` x3 at offsets 0, 50, 100 bp per
  element; two elements centred in the 16 kb read window, 2, 5, 8 or 12 kb apart;
  convergent (forward on the left element, reverse complement on the right, pointing
  inward), divergent, tandem, each element alone, and no insertion; 50 val background
  windows (chr8/chr20, GC-matched background class); 256 generated reads each (32 steps,
  bf16), 1,050 sequences, 18 min on the H100.
- Readouts: accessibility within +/-100 bp of each element; the two elements' raw and
  partial correlation across reads given each read's labelling level (its accessibility over
  the window excluding +/-300 bp around both positions); read-level-adjusted maps averaged
  over backgrounds.

| d | element access | pair partial r | element vs empty site | convergent - divergent |
|---|---|---|---|---|
| 2 kb | 0.30 | 0.19-0.21 | 0.00-0.02 | +0.014 +/- 0.012 |
| 5 kb | 0.29 | 0.18-0.19 | 0.01-0.02 | -0.003 +/- 0.012 |
| 8 kb | 0.28 | 0.17-0.19 | 0.00-0.01 | +0.016 +/- 0.011 |
| 12 kb | 0.27 | 0.14-0.15 | 0.00-0.02 | +0.010 +/- 0.010 |

Background accessibility at the same positions is about 0.12. A partner raises an element's
accessibility by only 0.002-0.009.

![Inserted CTCF pairs](fiberseq-16kb-pilot-labelling-efficiency/ctcf_insertion.png)

*Two inserted CTCF elements co-vary on generated molecules at every distance, with an
off-diagonal spot in the averaged maps, regardless of orientation; one element and an empty
site do not.*

Interpretation (inferred): the coupling at 12 kb has to pass through the bottleneck
attention, so the model does represent long-range dependence between strong elements; its
size is close to the observed CTCF-CTCF partial correlation (+0.136, real sites are weaker
than consensus triplets). Orientation regulates loop extrusion and contacts, with no
strong expectation for co-accessibility, so its absence may be biology or may mean the
model treats CTCF elements as generic strong openers. Untested control: the same insertion
with a strong non-CTCF opener (e.g. Mola's IRF cluster) and mixed pairs; observed non-CTCF
active pairs barely exceed the null (+0.029), so a model that learned the biology should
couple them much more weakly.

## Per-molecule labelling efficiency: fit and validation

Model: logit P(m6A at site s on molecule m) = a_s + b_m. a_s comes from the pileup of all
molecules at the site (smoothed logit of its m6A fraction, sites with >= 10 calling
molecules), so it absorbs the region's activity and the context bias at that site. b_m is
the molecule's offset, fitted on its background sites only (outside +/-500 bp of IDR peaks
and active cCREs) over its whole aligned span, by solving sum_k n_k sigma(a_k + b) = m6A
calls over a 96-bin histogram of a_s. Two passes per chromosome over the lane files, 4
processes, about an hour genome-wide; 99.98% of molecules fitted.

![Efficiency validation](fiberseq-16kb-pilot-labelling-efficiency/efficiency_validation.png)

*(a) b_m is unimodal with a low-labelling shoulder near -2 (the shoulder of the per-read
m6A-rate note); (b) it tracks the raw molecule rate (Spearman 0.89); (c) nucleosome
footprints shrink from 158 to 118 bp and merged > 250 bp stretches fall from 10.5% to 0.7%
across its quintiles, as the linker note found with raw rates; (d) it carries almost no
region activity (Spearman +0.05 on independent molecules).*

On the same molecules, window accessibility and mean b_m correlate at +0.24, but that is
sampling: the spread of window means (sd 0.04) is what ~130 molecules per window give, and a
window that happens to get better-labelled molecules also looks more accessible. Measured on
the other half of each window's molecules it drops to +0.05 (12,303 val windows). Active and
background windows have the same mean b_m (-0.173 / -0.169). A linker-based estimate was
not needed.

CRE pairs with b_m as a covariate (all 493 CTCF-CTCF pairs, 2,000 CTCF-other, 4,000 other,
3,000 random; ~2 min):

| pairs | raw r | given window level | given b_m | given both | p < 0.001 (both) |
|---|---|---|---|---|---|
| CTCF-CTCF | +0.309 | +0.136 | +0.119 | +0.112 | 4.5% |
| CTCF-other | +0.259 | +0.053 | +0.060 | +0.044 | 1.3% |
| other active | +0.251 | +0.031 | +0.057 | +0.027 | 1.0% |
| random positions | +0.313 | +0.017 | +0.037 | +0.005 | 0.8% |

Both covariates give the cleanest null. 69 of 6,493 sampled pairs pass FDR < 0.05; the
strongest have r 0.5-0.67 over 50-105 shared molecules, CTCF and not, e.g.
chr12:56,118,170-56,126,332 (CTCF-CTCF, 8.2 kb, r +0.62) and chrX:47,556,265-47,561,107
(4.8 kb, r +0.67).

## Training-design discussion (recorded in FIGWIG_TODO.md, not yet built)

- Batch size: no lower lr or accumulation; a 16 kb step holds 2x the read positions of a
  2 kb step, and Muon normalises update size.
- Encoding DNA once for more reads per region gains little (encoder ~6% of a step; 8 x 128
  reads 0.34 vs 0.38 s/step) and 16 kb regions have only ~61 qualifying reads.
- Curriculum / mixing: weights are length- and flank-independent (the 2 kb model loads
  strictly into a 16 kb config), so warm starts work; preferred: mixed-length batches
  (short element-centred + long), possibly ramped. CRE loss weighting (calls in IDR peaks:
  15.1% of the 2 kb model's centre, 4.2% of a 16 kb window, so w ~ 4) with `val/loss_cre` /
  `val/loss_bg` tracking; region weights for multi-CRE windows (>= 3 active cCREs in 32.8%
  of 16 kb windows; >= 2 CTCF-bound in 2.5%). Jitter (the advisor's point): it is 66 bp now;
  for 1 kb training widen it to +/-512 bp. Segment training (U-Net only on distant CRE
  segments, attention across them) is a later option, ~2.7x cheaper, with a train/inference
  mismatch.

## Labelling efficiency and enzyme bias (design agreed)

Per-molecule labelling efficiency is technical and wide (see the related notes), and the
m6A context bias is strong (+/-3 bp, 64x odds 5th-95th percentile, the same in open and
closed chromatin and on both strands).

1. **Per-molecule efficiency b_m without region information:**
   logit P(m6A at site s on molecule m) = a_s + b_m, with a_s from the pileup of all molecules
   (absorbs region activity and context bias) and b_m fitted on the molecule's background
   sites only (outside +/-500 bp of IDR peaks and active cCREs) over its whole span. Done
   genome-wide and validated (section above).
2. **Model:** the network outputs q = P(accessible); P(m6A) = q * e + (1 - q) * eps with
   e = sigma(c_context + b_m), c from a frozen +/-3 bp table computed in the model from the DNA
   (RC- and jitter-safe). The user's idea, adopted: the same quantities also as conditioning,
   so the network can read its input read correctly (an m6A at an inefficient context on a
   lightly labelled read is strong evidence of accessibility): b_m as per-read FiLM, e (or c)
   as an input channel. Interpretability: attributions on q are free of enzyme bias, and
   generation can be run at a fixed efficiency. Check: with all calls masked q must not
   depend on b_m. Caveat: non-local labelling effects (unmethylated short linkers on light
   reads) are not per-site and would leak into q.
3. Filtering: drop molecules below ~2%; no 8% ceiling.

## Commands

All from the H100 allocation 51969713 with the `~/repos/SMFNet_unet_figwig` pixi env and
`PYTHONPATH=~/repos/SMFNet_unet_figwig_runs/src` (see each script's header):

```bash
bash $SCRATCH/smfnet_figwig/smoke_16kb.sh 32 bs32            # batch-size smoke tests
nohup bash $SCRATCH/smfnet_figwig/train_k9_dilated_fiber_16kb.sh > ...out 2>&1
nohup bash $SCRATCH/smfnet_figwig/eval_16kb_val.sh <epoch24 ckpt> pilot25 > ...out 2>&1
python rescore_centre.py --long <16 kb ckpt> --short <2 kb ckpt> --out ... --mask-outside-centre
bash $SCRATCH/smfnet_figwig/eval_16kb_steps.sh               # 32 / 64 / 128 steps, 300 regions
python compare_samples.py --samples <samples.h5> --out <dir> --elements candidates.bed
python element_pairs_partial.py                               # CRE pairs, ~30 min on CPU
python ctcf_insertion.py --checkpoint <16 kb ckpt> --out $SCRATCH/smfnet_figwig/ctcf_pairs/run50
python analyse_ctcf_insertion.py $SCRATCH/smfnet_figwig/ctcf_pairs/run50
python molecule_efficiency.py <figwig dir> <out dir> chr1 ... # 4 processes in parallel
python validate_efficiency.py; python split_half_region_check.py; python plot_efficiency_validation.py
python element_pairs_efficiency.py                            # subsampled pairs with b_m, ~2 min
```

The hg38 FASTA the Fiber-seq jobs used
(`/scratch/groups/btrippe/brian/projects/SMF/fiberseq/reference/hg38.fa`) was removed on
2026-11-11; UCSC hg38 was re-downloaded to `/scratch/groups/btrippe/ndiamant/reference/hg38.fa`
(its index is byte-identical to the old one). `$SCRATCH/smfnet_paper_data/reference/hg38.fa`
is now a dangling symlink.

# Related Notes

- [GM12878 Fiber-seq through the figwig loader, U-Net vs CNN](../2026-11-11/fiberseq-figwig-unet-vs-cnn.md):
  the 2 kb model this pilot is compared with.
- [Per-molecule m6A rate in GM12878 Fiber-seq](fiberseq-per-read-m6a-rate.md): the rate
  distribution and filter options behind the labelling discussion.
- [Fiber read m6A rate evaluation (smf-hmm)](../../smf-hmm/2026-11-12/fiber-read-m6a-rate-evaluation.md):
  model states track labelling; drop reads below 2%, no 8% ceiling.
- [Fiber-seq m6A sequence-context bias (smf-hmm)](../../smf-hmm/2026-10-30/fiberseq-m6a-sequence-context-bias.md):
  the +/-3 bp table for the emission head.
- [Fiber HMM enzyme-bias diagnostic (smf-hmm)](../../smf-hmm/2026-10-30/fiber-hmm-enzyme-bias-diagnostic.md):
  how a DNA-conditioned model absorbs enzyme bias without context-offset emissions.
- [Fiber other_long footprints and linkers (smf-hmm)](../../smf-hmm/2026-11-11/fiber-other-long-footprint-reads-and-linkers.md):
  the non-local labelling effect (light reads merge footprints).
- [Span inpainting](../2026-11-10/inpainting.md): the architecture's configuration.

# Open Questions

- Are the strongest co-accessible pairs (FDR < 0.05) known co-regulated elements, and does
  the model reproduce them on their own sequences?
- Does the emission head with efficiency conditioning restore the between-molecule spread
  and keep q invariant to b_m at full mask?
- Does the 16 kb model reproduce the CTCF-CTCF partial correlations (needs regions chosen
  to contain such pairs)? Is the 10-20% short-range deficit a sampler or model limit
  (remasking test)?
- Is the inserted-pair coupling CTCF-specific? (Insert a strong non-CTCF opener, and mixed
  pairs.) Does orientation matter for any readout once labelling is modelled?

# Sources

- Code: SMFNet-dev `nd_unet_figwig` 961d8ba; scripts in
  [fiberseq-16kb-pilot-labelling-efficiency/scripts/](fiberseq-16kb-pilot-labelling-efficiency/scripts/).
- Mola: Wang R, Jin J, Pott S, Zhou J. "Predicting single-molecule chromatin distribution
  from DNA sequence." bioRxiv 2026, doi 10.64898/2026.10.29.762992 (local PDF in
  ~/repos/SMFNet_public).
- IDR ATAC peaks ENCODE [ENCFF614SMH](https://www.encodeproject.org/files/ENCFF614SMH/);
  UCSC hg38 `hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`.
