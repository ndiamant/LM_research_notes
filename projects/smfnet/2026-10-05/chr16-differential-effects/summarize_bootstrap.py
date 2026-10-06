"""Summarize the parameter-specific bootstrap, not the permutation screen."""
import json
import os
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=Path(os.environ['SCRATCH'])/'differential_chr16_20261005'
df=pd.read_csv(ROOT/'shape_full/regions.tsv',sep='\t')
summary=pd.read_csv(ROOT/'shape_full/summary.tsv',sep='\t')
controls=pd.read_csv(ROOT/'shape_controls/regions.tsv',sep='\t')
direct=pd.read_csv(ROOT/'bootstrap_full/regions.tsv',sep='\t').set_index(['region','pair'])
meta=json.loads((ROOT/'shape_full/metadata.json').read_text())
order=['NP-MEL','NP-ES','NP-C2C12_fixed','MEL-ES','MEL-C2C12_fixed','ES-C2C12_fixed']
labels=[s.replace('_fixed','').replace('-',' / ') for s in order]
rows=[]
for pair in order:
    f=df[df.pair==pair];p=f.profile_q_global<=.05;c=f.coc_q_global<=.05
    row=dict(pair=pair,regions=len(f),profile_testable=int(f.profile_p.notna().sum()),
        profile_significant=int(p.sum()),coc_testable=int(f.coc_p.notna().sum()),
        coc_significant=int(c.sum()),both_significant=int((p&c).sum()),
        profile_BY=int((f.profile_q_global_by<=.05).sum()),coc_BY=int((f.coc_q_global_by<=.05).sum()),
        profile_q01=int((f.profile_q_global<=.01).sum()),coc_q01=int((f.coc_q_global<=.01).sum()))
    raw=direct.loc[pd.MultiIndex.from_frame(f[['region','pair']])].reset_index()
    row['profile_direct_significant']=int((raw.profile_q_global<=.05).sum())
    row['coc_direct_significant']=int((raw.coc_q_global<=.05).sum())
    for m,sig in [('profile',p),('coc',c)]:
        row.update({m+'_median_r':f[m+'_r'].median(),m+'_median_r_significant':f.loc[sig,m+'_r'].median(),
            m+'_median_absolute_rms_shape_significant':raw.loc[sig.to_numpy(),m+'_rms_difference'].median(),
            m+'_median_absolute_noise_rms_shape_significant':raw.loc[sig.to_numpy(),m+'_expected_noise_rms'].median()})
    rows.append(row)
table=pd.DataFrame(rows)
table.to_csv(HERE/'summary.tsv',sep='\t',index=False)
summary.to_csv(HERE/'detailed_summary.tsv',sep='\t',index=False)
pd.read_csv(ROOT/'shape_controls/summary.tsv',sep='\t').to_csv(HERE/'controls_summary.tsv',sep='\t',index=False)
pd.read_csv(ROOT/'bootstrap_full/summary.tsv',sep='\t').to_csv(HERE/'direct_difference_summary.tsv',sep='\t',index=False)
pd.read_csv(ROOT/'bootstrap_controls/summary.tsv',sep='\t').to_csv(HERE/'direct_difference_controls.tsv',sep='\t',index=False)
(HERE/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(15,4.7),gridspec_kw={'width_ratios':[1,1,.85]})
for ax,m,title,color in zip(axes[:2],['profile','coc'],['Profile Pearson: shape differences','Correlation of correlations: shape'],['#377eb8','#984ea3']):
    s=summary[summary.metric==m].set_index('pair').loc[order];x=np.arange(6)
    ax.bar(x,s.testable/s.regions*100,color='#e5e5e5',label='Testable')
    ax.bar(x,s.bh_global_05/s.regions*100,color=color,label='BH q ≤ 0.05')
    ax.scatter(x,s.by_global_05/s.regions*100,color='black',marker='_',s=130,label='BY q ≤ 0.05',zorder=4)
    ax.set_xticks(x,labels,rotation=45,ha='right');ax.set_ylim(0,105);ax.set_title(title)
    ax.set_ylabel('% of 5,766 shared regions')
axes[0].legend(frameon=False,fontsize=9,loc='lower left')
for m,color,label in [('profile','#377eb8','Profile'),('coc','#984ea3','Correlations')]:
    axes[2].hist(controls[m+'_p'].dropna(),bins=np.linspace(0,1,11),density=True,histtype='step',lw=2,color=color,label=label)
axes[2].axhline(1,color='gray',ls='--',lw=1)
axes[2].set(xlabel='Multiplier-bootstrap p-value',ylabel='Density',title='Same-cell-type split controls',xlim=(0,1))
axes[2].legend(frameon=False,fontsize=9)
fig.suptitle('Observed dataset differences on chr16 — exploratory shape tests',fontsize=14)
fig.tight_layout();fig.savefig(HERE/'results.png',dpi=180);fig.savefig(HERE/'results.pdf')

lines=['# Differential-effect feasibility on chr16','','This experiment compares observed reads only. No models were sampled or evaluated.','','## Results','',
    '| Cell types | Profile testable | Profile significant | CoC testable | CoC significant | Both significant |',
    '|---|---:|---:|---:|---:|---:|']
for r in rows:
    lines.append(f"| {r['pair'].replace('_fixed','')} | {r['profile_testable']:,} | {r['profile_significant']:,} | {r['coc_testable']:,} | {r['coc_significant']:,} | {r['both_significant']:,} |")
lines.extend(['','All comparisons start with the same 5,766 regions. Significant means BH q ≤ 0.05 across all 34,596 region–cell-type-pair tests, separately for each metric; untestable entries receive p=1 for correction.',
    '', '**The headline tests target shape differences visible to profile Pearson and correlation-of-correlations.** Each vector is centered and normalized to unit length. The statistic is squared distance between these normalized vectors, exactly 2*(1-Pearson). Separate direct-difference tests also capture offset/scale differences; their counts are saved in direct_difference_summary.tsv.',
    '', '![Results](results.png)', '', '## Calibration and sensitivity',''])
for m in ['profile','coc']:
    n=int(controls[m+'_p'].notna().sum());raw=int((controls[m+'_p']<.05).sum());sig=int((controls[m+'_q_global']<=.05).sum())
    lines.append(f'- {m}: {raw}/{n} ({raw/n:.1%}) same-cell-type split controls have nominal p < 0.05; {sig} survive BH correction. Four cell types × 150 random regions, split into disjoint molecule halves.')
    significant=df[m+'_q_global']<=.05;strict=df[m+'_q_global_by']<=.05
    lines.append(f"- {m}: {df.loc[significant,'region'].nunique():,} distinct regions are BH-significant in at least one contrast; {int(strict.sum()):,} region–contrast tests survive the more conservative Benjamini–Yekutieli correction for arbitrary dependence.")
lines.extend(['', 'The shape-profile test has mild nominal inflation (7.2% versus 5%) in the split controls, so these are exploratory counts, not a guarantee of exact 5% FDR. Shape-CoC controls reject at 5.1%. The direct-difference tests reject at 4.9% for both metrics; neither direct nor shape controls has a BH discovery. BY and q≤0.01 counts provide sensitivity summaries but do not repair miscalibrated marginal p-values.', '', '## Design', '',
    '- Input: combined GCH SMF HDF5 datasets referenced by NP, MEL, ES, and C2C12_fixed saved configs; validation chromosome chr16, mm10. Evaluate the central [512,1536) bp of each 2,048-bp region, matching the models’ 1,024-bp output window.',
    '- Join exact shared region IDs and positions. Profiles require ≥30 calls per site in both groups and ≥10 sites. Pairwise analysis requires separation ≤500 bp, ≥50 jointly covering molecules in both groups, ≥5 calls of each state at each member of the pair in each group’s jointly covering molecules, and ≥10 eligible pairs. The profile-level coverage filter is also applied before building pairs.',
    '- Shape null: the population vectors have the same centered normalized direction (Pearson=1). Profile vectors contain site means; CoC vectors contain pairwise-complete within-molecule Pearson correlations. Statistic: sum of squared differences of unit vectors. Direct-difference sensitivity analysis instead tests equality of unnormalized vectors.',
    '- Estimate each vector’s influence function separately within each cell type. Generate independent standard-normal multipliers for each molecule, sharing a molecule’s multiplier across every site and pair it contributes to. The squared length of the difference of the two perturbed vectors gives the null reference. Thus shared-molecule dependence between site pairs is retained, while different mean accessibility, read coverage, and noise levels can be retained between cell types.',
    '- For site means, influence contributions are valid*(call-mean)/coverage. For pairwise correlation rho=(p11-mx*my)/sqrt(vx*vy), use its derivatives with respect to mx, my, and p11, multiplied by pair-valid/joint-coverage. Apply sqrt(n/(n-1)) finite-count variance corrections. For shape, propagate each influence vector through centering and unit normalization: (center(psi)-u*(center(psi) dot u))/norm(center(theta)). See bootstrap.py for formulas.',
    '- Use 1,999 multiplier draws, reproducible region/contrast seeds, and p=(1+number of null statistics ≥ observed)/(2,000). No zero p-values. This is an asymptotic plug-in bootstrap, not an exact finite-sample test. Coverage/minority-call thresholds and null controls support its exploratory use.',
    '- Pairwise statistics were checked against the repository’s covariation.py and influence derivatives against finite differences. Input calls and positions are validated, including unique read IDs within each input region.', '',
    '- The complete assay-mask audit checked all 5,766 regions in all four input files: zero sites in the evaluated window violate the current model’s DGCH mask. See assay_mask_check.json and check_assay.py.', '',
    '## Why the initial permutation screen is not the headline result', '',
    'The initial screen shuffled whole reads within identical coverage masks and tested low cross-dataset Pearson agreement. It tests equality of full conditional read distributions. In a synthetic check with identical population correlation matrices but different mean accessibility, it rejected 30/30 cases for correlation-of-correlations. Pooling groups with different means changes the reference correlation structure. Same-cell-type controls alone do not reveal this failure mode.', '',
    'The revised direct-difference bootstrap rejected 7/100 cases in that synthetic check at nominal 5%; the shape bootstrap rejected 8/100. Both detected profile changes in 100/100. The direct test also detected a changed dependence structure. The permutation artifacts are retained for transparency, but their counts are not interpreted as discoveries of changed pairwise correlations.', '',
    '## Interpretation and limitations', '',
    '- These are pooled-dataset differences. Biological and batch effects remain inseparable. Pairwise Pearson changes can themselves be driven partly by changed marginal accessibility; a significant change is not an isolated mechanistic dependence effect.',
    '- Significant does not imply a large effect. summary.tsv includes median absolute RMS differences among shape discoveries, estimated sampling-noise RMS, and descriptive profile Pearson / correlation-of-correlations. Raw shape-run RMS columns are in normalized-vector units, not accessibility units.',
    '- All observed reads are used for this discovery screen. Future model evaluation should select regions on one read subset and score predictions on independent reads; counts will generally change after splitting.',
    '- Data-dependent coverage/minority-call eligibility and approximate influence-function inference warrant sensitivity checks before confirmatory use. Same-cell-type controls assess read-sampling calibration, not biological replication.',
    '- Overlapping regions/site pairs are dependent. Molecule multipliers retain within-region dependence; BY counts supplement BH for dependence between region-level tests. A 0.0005 p-value floor limits stringent discoveries. Each metric and shape/direct variant is a separate testing family.', '',
    '## Reproduction and artifacts', '', 'Run from the SMFNet repository with load_smf active:', '', '```bash',
    'python analyses/differential_chr16_20261005/bootstrap.py --workers 2 --output "$SCRATCH/differential_chr16_20261005/bootstrap_full"',
    'python analyses/differential_chr16_20261005/bootstrap.py --workers 2 --limit 150 --controls --output "$SCRATCH/differential_chr16_20261005/bootstrap_controls"',
    'python analyses/differential_chr16_20261005/bootstrap.py --workers 6 --shape --output "$SCRATCH/differential_chr16_20261005/shape_full"',
    'python analyses/differential_chr16_20261005/bootstrap.py --workers 1 --limit 150 --controls --shape --output "$SCRATCH/differential_chr16_20261005/shape_controls"',
    'python analyses/differential_chr16_20261005/summarize_bootstrap.py', '```', '',
    f'Headline per-region outputs and metadata: `{ROOT}/shape_full/`. Shape controls: `{ROOT}/shape_controls/`. Direct tests: `{ROOT}/bootstrap_full/` and `{ROOT}/bootstrap_controls/`. Initial permutation screen: `{ROOT}/full/` and `{ROOT}/controls/`. Scratch artifacts are reproducible and temporary; this directory retains the code, report, compact summaries, and figures.', '',
    'The full-distribution null of label permutation is described in the [SciPy permutation-test documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html).'])
(HERE/'README.md').write_text('\n'.join(lines)+'\n')
print(table.to_string(index=False))
