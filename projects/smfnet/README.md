# SMFNet

## Notes

- 2026-11-09: [Private dev repo kundajelab/SMFNet-dev next to the public kundajelab/SMFNet, for unreleased work with Till](2026-11-09/private-dev-repo.md)
- 2026-11-09: [Kernel width vs depth in the attention U-Net, roofline and measured costs, and a k=9 one-block-per-down-level run at matched receptive field](2026-11-09/kernel-width.md)
- 2026-11-08: [Do SMFNet models memorize training regions' DNA? Train vs val loss by mask rate, with own, wrong and mutated DNA](2026-11-08/memorization.md)
- 2026-11-07: [Patched U-Net (4 bp patches, 1 bp head) matches the attention U-Net early in training at 0.6-0.8x the step and 0.7x the memory; Muon+EMA vs AdamW+EMA final results](2026-11-07/patched-unet.md)
- 2026-11-06: [Attention U-Net training step 292 → 184 ms with a fused Triton norm kernel, channels-last activations and a Toeplitz output conv (H100)](2026-11-06/unet-training-speed.md)
- 2026-11-03: [U-Net SMF body with AlphaGenome-style blocks and per-level FiLM from frozen AlphaGenome features](2026-11-03/unet-smf-body.md)
- 2026-11-03: [Cheri blocks and Muon in SMFNet, the fused Triton kernel on PyTorch 2.5, memory, and Muon learning rate](2026-11-03/cheri-muon-smfnet.md)
- 2026-10-29: [Frozen AlphaGenome encoder as SMFNet's DNA encoder (ONT mESC), design and val results vs baseline and ATAC](2026-10-29/frozen-alphagenome-encoder.md)
- 2026-10-14: [PRISM speedup from shared DNA encoding and skipped wildtype bases (commit 66d6585) on H100](2026-10-14/prism-speedup-benchmark-h100.md)
- 2026-10-08: [Region count for the sampling-step appendix figure](2026-10-08/sampling-step-region-count.md)
- 2026-10-08: [What ends observed SMF runs, and whether to bridge missing gaps in run-length metrics](2026-10-08/run-length-missingness-gap-bridging.md)
- 2026-10-07: [GCH overfitting, regularization, EMA, and sampling-length experiments](2026-10-07/regularization-ema-sampling-experiments.md)
- 2026-10-06: [Why SMFNet predicts cell-type differences in pair correlations weakly](2026-10-06/covariation-differences-diagnostics.md)
- 2026-10-06: [Differential cell-type evaluation with held-out molecules, effect-size strata, and between/within-region scores](2026-10-06/differential-evaluation-held-out-molecules.md)
- 2026-10-05: [Chr16 feasibility of differential cell-type profiles and correlation patterns](2026-10-05/chr16-differential-effects-feasibility.md)
- 2026-10-05: [Fiber-seq GM12878 data preparation, design choices and pipeline run](2026-10-05/fiberseq-gm12878-data-preparation.md)
- 2026-10-05: [Flank baseline for profile metrics, logistic regression instead of a k-mer table](2026-10-05/flank-baseline-logistic-regression.md)
- 2026-09-29: [Train-set positional baselines and noise-matched correlation-of-correlations](2026-09-29/positional-covariation-baselines-and-noise-matched-corr-of-corr.md)
