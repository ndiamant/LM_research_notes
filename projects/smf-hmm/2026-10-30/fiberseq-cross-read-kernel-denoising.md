---
title: Cross-read kernel denoising of undermethylated Fiber-seq, and why own-read smoothing beats it
date: 2026-10-30
project: smf-hmm
agent: Claude Code
status: draft
sources:
  - smf_hmm repo, evaluation/datasets/smf/kernel_denoise/ (README.md, make_widget.py, _template.html, fiberseq_kernel_denoise.html), branch nd/experiments, commits f20556f, 623ae24 and b3580ca
  - /scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/{active.h5,manifest.tsv}
  - /scratch/users/diamant/kernel_denoise_proto/ (proto.py, proto2.py, proto3.py, sweep.js, extract.js, core_test.js, cmp.py, dom_test.js, prof.js)
  - smf_hmm repo, evaluation/datasets/smf/context_bias/adj_table_3_3.npz (commit d20f916), for the bias correction
  - smf_hmm repo, evaluation/datasets/smf/simple_HMM_tests/ (page style and call-raster rendering reused)
  - User discussion, 2026-10-30
tags: [smf-hmm, fiberseq, m6a, gm12878, denoising, non-local-means, smoothing, thinning-evaluation, widget]
---

# Summary

The user proposed denoising undermethylated GM12878 Fiber-seq calls with non-local means.
Cut each molecule into small tiles (about 10 bp) and replace each tile with a weighted average
of the same tile in the other molecules at the locus. The weights are a kernel similarity,
normalised to sum to one, and include the molecule itself. We prototyped it in numpy, then
built a self-contained interactive page where every parameter is live. The page includes a
binomial-thinning evaluation.

**Result: the cross-read denoiser beats the column mean easily, but never beats smoothing a
molecule's own calls.** At 50% thinning over 24 windows, the mean AUROC for recovering removed
m6A is 0.794 for the column mean, 0.878 for the best kernel setting and 0.890 for own-read
Gaussian smoothing (bw 15 bp). No setting of 96 wins on any of the 24 windows. The kernel also
shifts the window-wide m6A rate, which own-read smoothing does not (see Details).

# Key Points

- **The kernel has to respect one-sided noise.** Undermethylation makes a 1 reliable and a 0
  ambiguous, and plain Hamming counts the dominant 0–0 agreement as evidence. We added a
  *likelihood* kernel: the per-site pointwise mutual information of two calls under a shared
  accessibility state. With P(m6A|acc) 0.45, P(m6A|prot) 0.03 and P(acc) 0.3, the scores are
  m6A–m6A +0.925, mismatch −0.330 and no–no +0.051. It beats Hamming at every setting tried.
- **Similarity has to come from a wider context than the tile.** A 10 bp tile holds about 5 A/T
  sites, so the weights use ±C bp around it. The best C was 25–50 bp.
- **Self-weight from the kernel collapses.** A molecule's similarity to itself always
  dominates, giving AUROC 0.68 in the prototype. Setting the self-weight equal to the best
  other molecule's weight, as in non-local means, works.
- **Averaging the tile's m6A rate beats averaging site by site:** 0.880 vs 0.859 in the
  prototype, at τ 1, C 50, W 10.
- **The thinning test favors own-read smoothing (inferred).** Sites are thinned independently,
  so every removed m6A keeps unthinned same-molecule neighbours a few bp away. m6A comes in
  runs (linkers, NDRs), and that local evidence is worth more than other molecules' calls at the
  same site. AUROC is rank-only, so blurring a footprint's edges costs little. A fairer test
  would thin in 20–50 bp blocks per molecule, so own-read smoothing has nothing local to use.
  This has not been run yet.
- **The kernel weights bias the mean m6A rate.** Over 24 windows the raw rate is 0.144. The
  likelihood kernel raises it to 0.169 (×1.18), or ×1.24 with self excluded. Hamming lowers it
  to 0.120 (×0.83). Neighbours are chosen partly by how methylated they are: the likelihood
  kernel rewards shared m6A, Hamming shared no-m6A. Own-read RBF preserves the mean (0.168 vs
  0.166 on one window).

# Details

## Data

`h5_ont_polarity/active.h5` encodes 1 = unmethylated (protected), 0 = m6A and −1 = no call;
the attribute `smf_encoding` says so. The page converts this to the `simple_HMM_tests`
convention: 0 no call, 1 m6A, 2 no m6A.

`make_widget.py` samples 8 active windows each from `CTCF-only,CTCF-bound`, `PLS` and `dELS`,
with `full_span >= 100` and seed 0. It keeps molecules with at least 80% of the window's sites
called, caps them at 100 and orders them by haplotype. The h5 `strand` field is the read's
alignment orientation, not a biological strand. Each CCS molecule reports both strands, with
m6A already merged onto reference A/T positions.

## Prototype numbers (numpy, 12 windows, 50% thinning)

The windows are 4 each of CTCF-only, PLS and dELS. The table gives mean AUROC.

| method | AUROC |
| --- | --- |
| column mean | 0.803 |
| own-read RBF, bw 10 / 15 / 25 / 40 bp | 0.889 / 0.889 / 0.881 / 0.857 |
| likelihood, τ 3, C 100, W 10, self = max, per site | 0.863 |
| same, self = kernel | 0.679 |
| Hamming, τ 3, C 100 | 0.844 |
| likelihood, τ 1, C 50, W 10, tile-rate average | 0.880 |
| kernel average of RBF-smoothed molecules, best | 0.888, still not above RBF alone |

Sharpening τ drives the denoiser toward self-only, which converges on own-read smoothing but
never passes it.

## Page and verification

The page `fiberseq_kernel_denoise.html` is 2.7 MB with data inlined and runs in any browser.
It has four views:

- observed and denoised rasters, where clicking a row selects a molecule;
- a region-average panel: per-site mean m6A rate with no smoothing, with own-read RBF and with
  the kernel, plus window-wide means;
- the selected molecule's traces, a molecule × tile weight heatmap and the effective neighbour
  count per tile;
- a thinning-evaluation dot plot, with an "Evaluate all regions" table.

Checks done:

- The JS denoiser matches the numpy prototype to float32 precision, about 3e-8. The only
  differences are conventions: NaN vs 0 where no molecule has a call, and the RBF cutoff at
  4·bw.
- A node run with a mock DOM exercised every handler, with about 2.7M draw calls and none NaN
  or off-canvas.
- A full update costs about 0.3 s of compute in V8, or 0.6 s with 1 bp tiles.
- No browser is available on Sherlock, so the rendered layout was not inspected by the agent.

The 24-window sweep (`sweep.js`, thinning 0.5) covered τ ∈ {0.25, 0.5, 1, 2},
C ∈ {25, 50, 100, 200} and W ∈ {5, 10, 20} for each kernel, with tile averaging and
self = max. The best was likelihood, τ 1, C 25–50, W 10, at 0.878 vs RBF 0.890 and 0/24 wins.
The page defaults to that setting.

## Commands

Environment: `PY=/home/groups/btrippe/diamant/miniforge/envs/hmm/bin/python`, on a Sherlock
GPU allocation (sh04-03n07, Slurm job 50179010). Node is at
`/share/software/user/open/nodejs/25.3.0/bin/node`.

```bash
# build the page (reads active.h5 + manifest.tsv)
cd evaluation/datasets/smf/kernel_denoise && $PY make_widget.py
# numpy prototypes (12 windows)
cd $SCRATCH/kernel_denoise_proto && $PY -W ignore proto.py 0.5 && $PY -W ignore proto2.py && $PY -W ignore proto3.py
# extract the page's compute core, then check it against numpy and sweep it
node extract.js <path>/fiberseq_kernel_denoise.html && node prof.js && $PY -W ignore cmp.py
node sweep.js
node dom_test.js <path>/fiberseq_kernel_denoise.html   # mock-DOM smoke test
```

Profile and check the core through `extract.js` → `core_mod.js`. A direct `eval` of the
script blocks V8 optimisation and inflated the timings about tenfold.

## Correcting own-read smoothing for context bias (commit `b3580ca`)

Added after the [context-bias experiment](fiberseq-m6a-sequence-context-bias.md).

- **Efficiency per site.** `make_widget.py` looks up each site's 3+3 context in hg38,
  reverse-complemented at T sites. It takes a relative efficiency `e` from the adjusted table
  `adj_table_3_3.npz`: `sigmoid(f + g[top active bin])`, scaled to mean 1 over the payload.
  The lookup matches the fit's own `ctx_id(3, 3)` encoding on 5000 chr21 sites.
- **Bias-corrected RBF.** `â = Σ w·m6A / Σ w·e` estimates accessibility under
  `P(m6A) = a · e(context)`. The ratio form avoids the variance blow-up of reweighting each call
  by 1/e where e is about 0.002. Because the thinning target is a call, it is scored as `â · e`.

Thinning at 0.5, bw 15 bp, 24 windows, mean AUROC:

| score | all 24 | 3 windows on the table's held-out chromosomes |
| --- | --- | --- |
| `e` alone | 0.731 | 0.717 |
| own-read RBF | 0.890 | 0.898 |
| RBF × e (control) | 0.916 | 0.917 |
| corrected RBF × e | 0.917 | 0.918 |

- **Nearly all of the gain is the factor of `e`, not the correction.** Removed m6A sits mostly
  at efficient sites. Corrected × e beats the RBF × e control on 24 of 24 windows, but by only
  +0.001, so this AUROC barely measures accessibility quality.
- **The correction does change the accessibility estimate.** Averaged over molecules, the
  corrected curve's correlation with Gaussian-weighted local mean efficiency falls from 0.62
  to 0.36, at every bw from 10 to 25. The remaining 0.36 may be real GC-accessibility
  covariation, or under-correction by the adjusted table (unresolved).
- **It is not leakage.** The table was fit on chromosomes containing 21 of the 24 windows, but
  the 3 windows on held-out chr8 and chr9 match the full set.
- **The page now has the RBF × e control** in its dot plot and all-regions table, so the
  comparison is not misread.

# Related Notes

- [Fiber-seq GM12878 data preparation for SMFNet](../../smfnet/2026-10-05/fiberseq-gm12878-data-preparation.md):
  Builds the `active.h5` / `background.h5` windows used here, and documents the strand
  symmetry that lets one treat A and T sites alike.
- [Posterior predictive checks of the fixed-blocks SMF HMM](../2026-10-22/posterior-predictive-checks.md):
  Argues for held-out site predictions as the informative check. The thinning evaluation here
  is the same idea applied to denoisers.
- [m6A sequence-context bias in GM12878 Fiber-seq](fiberseq-m6a-sequence-context-bias.md):
  Follow-up from the same day. It finds about a 100-fold context-dependent spread in the m6A
  rate, which the flat P(m6A|acc) used by this kernel ignores, and supplies the 3+3 table used
  for the bias correction above.

# Open Questions

- Does cross-read borrowing help under **block** thinning (20–50 bp per molecule), or on a
  resolution-sensitive target such as the protected dip at ChIP-supported CTCF motifs?
- Can the mean-rate bias be removed, for example by normalising each molecule's overall m6A
  rate before computing weights, without losing the AUROC gain over the column mean?
- Would context-dependent emissions in the likelihood kernel change the ranking against
  own-read smoothing?
- Is the bias-corrected accessibility closer to the truth? The call-level thinning AUROC cannot
  tell. Block thinning, or flat protection at ChIP-supported motifs instead of sequence-driven
  dips, could.

# Sources

- smf_hmm `evaluation/datasets/smf/kernel_denoise/` on branch `nd/experiments` (commits
  `f20556f`, `623ae24`, `b3580ca`): `README.md`, `make_widget.py`, `_template.html`,
  `fiberseq_kernel_denoise.html`.
- `/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/active.h5` and
  `manifest.tsv`.
- Prototype and verification scripts: `/scratch/users/diamant/kernel_denoise_proto/`.
- User discussion, 2026-10-30.
