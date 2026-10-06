"""Parameter-specific whole-molecule multiplier bootstrap feasibility screen.

The global test uses squared L2 differences of site means or pairwise Pearson
correlations. Sampling covariance is estimated separately within each group,
using the influence function of each statistic. Shared Gaussian weights per
molecule preserve dependence across sites/pairs. This is an asymptotic test,
not an exact permutation test; calibration checks are essential.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import itertools
import json
import multiprocessing
import os
from pathlib import Path
import time

import h5py
import numpy as np
import pandas as pd
import torch
import yaml
import investigate as base

FILES = {}


def init(paths):
    torch.set_num_threads(1)
    FILES.update({c:h5py.File(p,'r') for c,p in paths.items()})


def features(a, rows, cols):
    v = (a != -1).double()
    x = a.clamp_min(0).double()
    n = v.sum(0)
    mu = x.sum(0)/n
    psi_mean = (x-v*mu)/n * (n/(n-1)).sqrt()
    j = v[:,rows]*v[:,cols]
    nj = j.sum(0)
    xj, yj = x[:,rows]*j, x[:,cols]*j
    mx, my = xj.sum(0)/nj, yj.sum(0)/nj
    xy = xj*yj
    p11 = xy.sum(0)/nj
    vx, vy = mx*(1-mx), my*(1-my)
    scale = (vx*vy).sqrt()
    rho = (p11-mx*my)/scale
    dx = -my/scale-rho*(1-2*mx)/(2*vx)
    dy = -mx/scale-rho*(1-2*my)/(2*vy)
    psi_pair = ((xy-j*p11)/scale + dx*(xj-j*mx) + dy*(yj-j*my))/nj * (nj/(nj-1)).sqrt()
    return mu, psi_mean, rho, psi_pair


def test_vectors(theta_a, psi_a, theta_b, psi_b, resamples, seed, shape=False):
    if shape:
        # For centered unit vectors u,v, ||u-v||^2 = 2*(1-Pearson).
        # Propagate whole-molecule influence through centering/normalization.
        def normalize(theta, psi):
            centered=theta-theta.mean()
            norm=centered.norm()
            if norm<=1e-12:
                return None
            unit=centered/norm
            projected=psi-psi.mean(1,keepdim=True)
            projected=(projected-(projected@unit)[:,None]*unit)/norm
            return unit,projected
        aa,bb=normalize(theta_a,psi_a),normalize(theta_b,psi_b)
        if aa is None or bb is None:
            return dict(p=np.nan,r=np.nan,rms_difference=np.nan,mean_difference=np.nan,
                expected_noise_rms=np.nan,null_rms95=np.nan)
        theta_a,psi_a=aa;theta_b,psi_b=bb
    diff = theta_a-theta_b
    observed = diff.square().sum()
    psi = torch.cat([psi_a, -psi_b]).float()
    gen = torch.Generator().manual_seed(seed)
    exceed = 0
    noise = []
    for start in range(0,resamples,128):
        count = min(128,resamples-start)
        weights = torch.randn((count,len(psi)),generator=gen)
        samples = weights@psi
        null = samples.square().sum(1)
        exceed += int((null>=observed.float()-1e-10).sum())
        noise.extend(null.tolist())
    r,_ = base.pearson_rows(theta_a[None],theta_b[None])
    return dict(p=(1+exceed)/(resamples+1), r=float(r[0]),
        rms_difference=float(diff.square().mean().sqrt()), mean_difference=float(diff.mean()),
        expected_noise_rms=float((psi.square().sum()/len(diff)).sqrt()),
        null_rms95=float(np.sqrt(np.quantile(noise,.95)/len(diff))))


def compare(pa,a,pb,b,resamples=1999,site_coverage=30,joint_coverage=50,minor=5,seed=0,shape=False):
    pos,ia,ib=np.intersect1d(pa,pb,return_indices=True)
    a,b=a[:,ia],b[:,ib]
    keep=((a!=-1).sum(0)>=site_coverage)&((b!=-1).sum(0)>=site_coverage)
    pos,a,b=pos[keep],a[:,keep],b[:,keep]
    a,b=a[(a!=-1).any(1)],b[(b!=-1).any(1)]
    row=dict(n_sites=len(pos),reads_a=len(a),reads_b=len(b),n_pairs=0)
    if len(pos)<10:
        return dict(row,status='fewer_than_10_sites')
    a,b=torch.tensor(a),torch.tensor(b)
    i,j=np.triu_indices(len(pos),1)
    keep=(pos[j]-pos[i]<=500)
    i,j=torch.tensor(i[keep]),torch.tensor(j[keep])
    eligible=torch.ones(len(i),dtype=torch.bool)
    for x in [a,b]:
        valid=(x!=-1)
        joint=valid[:,i]&valid[:,j]
        n=joint.sum(0)
        sx=((x[:,i]==1)&joint).sum(0)
        sy=((x[:,j]==1)&joint).sum(0)
        eligible &= (n>=joint_coverage)&(sx>=minor)&(sy>=minor)&(n-sx>=minor)&(n-sy>=minor)
    i,j=i[eligible],j[eligible]
    row['n_pairs']=len(i)
    ma,pma,ra,pra=features(a,i,j)
    mb,pmb,rb,prb=features(b,i,j)
    for name,ta,psia,tb,psib in [('profile',ma,pma,mb,pmb),('coc',ra,pra,rb,prb)]:
        if len(ta)<10:
            continue
        vals=test_vectors(ta,psia,tb,psib,resamples,seed+(name=='coc'),shape=shape)
        row.update({name+'_'+k:v for k,v in vals.items()})
    return dict(row,status='ok')


def worker(task):
    region,resamples,site_coverage,joint_coverage,controls,shape=task
    data={c:base.load_region(f,region) for c,f in FILES.items()}
    pairs=[(c,c) for c in base.CELL_TYPES] if controls else list(itertools.combinations(base.CELL_TYPES,2))
    result=[]
    for ca,cb in pairs:
        pa,a=data[ca];pb,b=data[cb]
        seed=int.from_bytes(hashlib.sha256(f'{region}/{ca}/{cb}/bootstrap20261005'.encode()).digest()[:4],'little')
        if controls:
            idx=np.random.default_rng(seed).permutation(len(a))
            a,b=a[idx[:len(idx)//2]],a[idx[len(idx)//2:]]
        row=dict(region=region,pair=ca+'-'+cb)
        row.update(compare(pa,a,pb,b,resamples,site_coverage,joint_coverage,seed=seed,shape=shape))
        result.append(row)
    return result


def summarize(df,output):
    rows=[]
    for metric in ['profile','coc']:
        df[metric+'_q_global']=base.adjust(df[metric+'_p'])
        df[metric+'_q_global_by']=base.adjust(df[metric+'_p'],by=True)
        for pair,ix in df.groupby('pair').groups.items():
            df.loc[ix,metric+'_q_pair']=base.adjust(df.loc[ix,metric+'_p'])
            f=df.loc[ix]
            sig=f[metric+'_q_global']<=.05
            rows.append(dict(pair=pair,metric=metric,regions=len(f),testable=int(f[metric+'_p'].notna().sum()),
                nominal_p05=int((f[metric+'_p']<.05).sum()),bh_global_05=int(sig.sum()),
                bh_pair_05=int((f[metric+'_q_pair']<=.05).sum()),by_global_05=int((f[metric+'_q_global_by']<=.05).sum()),
                median_r=f[metric+'_r'].median(),median_rms_significant=f.loc[sig,metric+'_rms_difference'].median(),
                median_r_significant=f.loc[sig,metric+'_r'].median()))
    df.to_csv(output/'regions.tsv',sep='\t',index=False)
    summary=pd.DataFrame(rows)
    summary.to_csv(output/'summary.tsv',sep='\t',index=False)
    print(summary.to_string(index=False),flush=True)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--limit',type=int)
    ap.add_argument('--resamples',type=int,default=1999)
    ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--site-coverage',type=int,default=30)
    ap.add_argument('--joint-coverage',type=int,default=50)
    ap.add_argument('--controls',action='store_true')
    ap.add_argument('--shape',action='store_true',help='Test centered normalized vectors: statistic = 2*(1-Pearson).')
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    paths={};keys={}
    for c in base.CELL_TYPES:
        cfg=yaml.safe_load((Path(os.environ['SCRATCH'])/'smf_models'/c/'.hydra/config.yaml').read_text())
        paths[c]=cfg['data']['cfg']['h5_path']
        with h5py.File(paths[c],'r') as f:
            keys[c]={r for r in f['smf_mat'] if r.startswith('chr16:')}
    regions=sorted(set.intersection(*keys.values()));np.random.default_rng(20261005).shuffle(regions)
    if args.limit:regions=regions[:args.limit]
    metadata=dict(chromosome='chr16',window=[512,1536],resamples=args.resamples,
        site_coverage=args.site_coverage,joint_coverage=args.joint_coverage,min_minority_calls=5,
        min_sites=10,min_pairs=10,max_pair_distance=500,common_regions=len(set.intersection(*keys.values())),
        analyzed_regions=len(regions),controls=args.controls,paths=paths,shape=args.shape,
        null='equal centered normalized vectors (Pearson=1)' if args.shape else 'equal site-mean vectors or equal pairwise Pearson-correlation vectors, respectively',
        method='separate-group influence-function Gaussian whole-molecule multiplier bootstrap of squared L2 differences',
        multiple_testing='BH across all region-pair tests separately by metric; untestable p=1; BY sensitivity')
    (args.output/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    started=time.monotonic();results=[]
    tasks=[(r,args.resamples,args.site_coverage,args.joint_coverage,args.controls,args.shape) for r in regions]
    with ProcessPoolExecutor(max_workers=args.workers,initializer=init,initargs=(paths,),mp_context=multiprocessing.get_context('spawn')) as ex:
        with (args.output/'raw.jsonl').open('w') as stream:
            for n,batch in enumerate(ex.map(worker,tasks,chunksize=4),1):
                results.extend(batch)
                for row in batch:stream.write(json.dumps(row)+'\n')
                stream.flush()
                if n%50==0 or n==len(regions):print(f'{n}/{len(regions)} regions; {time.monotonic()-started:.1f}s',flush=True)
    summarize(pd.DataFrame(results),args.output)


if __name__=='__main__':main()
