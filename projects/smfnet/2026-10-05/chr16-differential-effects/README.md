# Differential-effect feasibility on chr16

This experiment compares observed reads only. No models were sampled or evaluated.

## Results

| Cell types | Profile testable | Profile significant | CoC testable | CoC significant | Both significant |
|---|---:|---:|---:|---:|---:|
| NP-MEL | 5,254 | 4,575 | 5,095 | 4,917 | 4,358 |
| NP-ES | 5,300 | 4,421 | 5,136 | 4,940 | 4,204 |
| NP-C2C12 | 5,139 | 4,202 | 4,911 | 4,594 | 3,848 |
| MEL-ES | 5,302 | 4,727 | 5,123 | 4,928 | 4,458 |
| MEL-C2C12 | 5,165 | 4,242 | 4,923 | 4,620 | 3,895 |
| ES-C2C12 | 5,177 | 4,236 | 4,929 | 4,609 | 3,888 |

All comparisons start with the same 5,766 regions. Significant means BH q ≤ 0.05 across all 34,596 region–cell-type-pair tests, separately for each metric; untestable entries receive p=1 for correction.

**The headline tests target shape differences visible to profile Pearson and correlation-of-correlations.** Each vector is centered and normalized to unit length. The statistic is squared distance between these normalized vectors, exactly 2*(1-Pearson). Separate direct-difference tests also capture offset/scale differences; their counts are saved in direct_difference_summary.tsv.

![Results](results.png)

## Calibration and sensitivity

- profile: 37/512 (7.2%) same-cell-type split controls have nominal p < 0.05; 0 survive BH correction. Four cell types × 150 random regions, split into disjoint molecule halves.
- profile: 5,325 distinct regions are BH-significant in at least one contrast; 19,136 region–contrast tests survive the more conservative Benjamini–Yekutieli correction for arbitrary dependence.
- coc: 23/453 (5.1%) same-cell-type split controls have nominal p < 0.05; 0 survive BH correction. Four cell types × 150 random regions, split into disjoint molecule halves.
- coc: 5,223 distinct regions are BH-significant in at least one contrast; 25,215 region–contrast tests survive the more conservative Benjamini–Yekutieli correction for arbitrary dependence.

The shape-profile test has mild nominal inflation (7.2% versus 5%) in the split controls, so these are exploratory counts, not a guarantee of exact 5% FDR. Shape-CoC controls reject at 5.1%. The direct-difference tests reject at 4.9% for both metrics; neither direct nor shape controls has a BH discovery. BY and q≤0.01 counts provide sensitivity summaries but do not repair miscalibrated marginal p-values.

## Design

- Input: combined GCH SMF HDF5 datasets referenced by NP, MEL, ES, and C2C12_fixed saved configs; validation chromosome chr16, mm10. Evaluate the central [512,1536) bp of each 2,048-bp region, matching the models’ 1,024-bp output window.
- Join exact shared region IDs and positions. Profiles require ≥30 calls per site in both groups and ≥10 sites. Pairwise analysis requires separation ≤500 bp, ≥50 jointly covering molecules in both groups, ≥5 calls of each state at each member of the pair in each group’s jointly covering molecules, and ≥10 eligible pairs. The profile-level coverage filter is also applied before building pairs.
- Shape null: the population vectors have the same centered normalized direction (Pearson=1). Profile vectors contain site means; CoC vectors contain pairwise-complete within-molecule Pearson correlations. Statistic: sum of squared differences of unit vectors. Direct-difference sensitivity analysis instead tests equality of unnormalized vectors.
- Estimate each vector’s influence function separately within each cell type. Generate independent standard-normal multipliers for each molecule, sharing a molecule’s multiplier across every site and pair it contributes to. The squared length of the difference of the two perturbed vectors gives the null reference. Thus shared-molecule dependence between site pairs is retained, while different mean accessibility, read coverage, and noise levels can be retained between cell types.
- For site means, influence contributions are valid*(call-mean)/coverage. For pairwise correlation rho=(p11-mx*my)/sqrt(vx*vy), use its derivatives with respect to mx, my, and p11, multiplied by pair-valid/joint-coverage. Apply sqrt(n/(n-1)) finite-count variance corrections. For shape, propagate each influence vector through centering and unit normalization: (center(psi)-u*(center(psi) dot u))/norm(center(theta)). See bootstrap.py for formulas.
- Use 1,999 multiplier draws, reproducible region/contrast seeds, and p=(1+number of null statistics ≥ observed)/(2,000). No zero p-values. This is an asymptotic plug-in bootstrap, not an exact finite-sample test. Coverage/minority-call thresholds and null controls support its exploratory use.
- Pairwise statistics were checked against the repository’s covariation.py and influence derivatives against finite differences. Input calls and positions are validated, including unique read IDs within each input region.

- The complete assay-mask audit checked all 5,766 regions in all four input files: zero sites in the evaluated window violate the current model’s DGCH mask. See assay_mask_check.json and check_assay.py.

## Why the initial permutation screen is not the headline result

The initial screen shuffled whole reads within identical coverage masks and tested low cross-dataset Pearson agreement. It tests equality of full conditional read distributions. In a synthetic check with identical population correlation matrices but different mean accessibility, it rejected 30/30 cases for correlation-of-correlations. Pooling groups with different means changes the reference correlation structure. Same-cell-type controls alone do not reveal this failure mode.

The revised direct-difference bootstrap rejected 7/100 cases in that synthetic check at nominal 5%; the shape bootstrap rejected 8/100. Both detected profile changes in 100/100. The direct test also detected a changed dependence structure. The permutation artifacts are retained for transparency, but their counts are not interpreted as discoveries of changed pairwise correlations.

## Interpretation and limitations

- These are pooled-dataset differences. Biological and batch effects remain inseparable. Pairwise Pearson changes can themselves be driven partly by changed marginal accessibility; a significant change is not an isolated mechanistic dependence effect.
- Significant does not imply a large effect. summary.tsv includes median absolute RMS differences among shape discoveries, estimated sampling-noise RMS, and descriptive profile Pearson / correlation-of-correlations. Raw shape-run RMS columns are in normalized-vector units, not accessibility units.
- All observed reads are used for this discovery screen. Future model evaluation should select regions on one read subset and score predictions on independent reads; counts will generally change after splitting.
- Data-dependent coverage/minority-call eligibility and approximate influence-function inference warrant sensitivity checks before confirmatory use. Same-cell-type controls assess read-sampling calibration, not biological replication.
- Overlapping regions/site pairs are dependent. Molecule multipliers retain within-region dependence; BY counts supplement BH for dependence between region-level tests. A 0.0005 p-value floor limits stringent discoveries. Each metric and shape/direct variant is a separate testing family.

## Reproduction and artifacts

Run from the SMFNet repository with load_smf active:

```bash
python analyses/differential_chr16_20261005/bootstrap.py --workers 2 --output "$SCRATCH/differential_chr16_20261005/bootstrap_full"
python analyses/differential_chr16_20261005/bootstrap.py --workers 2 --limit 150 --controls --output "$SCRATCH/differential_chr16_20261005/bootstrap_controls"
python analyses/differential_chr16_20261005/bootstrap.py --workers 6 --shape --output "$SCRATCH/differential_chr16_20261005/shape_full"
python analyses/differential_chr16_20261005/bootstrap.py --workers 1 --limit 150 --controls --shape --output "$SCRATCH/differential_chr16_20261005/shape_controls"
python analyses/differential_chr16_20261005/summarize_bootstrap.py
```

Headline per-region outputs and metadata: `/scratch/users/diamant/differential_chr16_20261005/shape_full/`. Shape controls: `/scratch/users/diamant/differential_chr16_20261005/shape_controls/`. Direct tests: `/scratch/users/diamant/differential_chr16_20261005/bootstrap_full/` and `/scratch/users/diamant/differential_chr16_20261005/bootstrap_controls/`. Initial permutation screen: `/scratch/users/diamant/differential_chr16_20261005/full/` and `/scratch/users/diamant/differential_chr16_20261005/controls/`. Scratch artifacts are reproducible and temporary; this directory retains the code, report, compact summaries, and figures.

The full-distribution null of label permutation is described in the [SciPy permutation-test documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html).
