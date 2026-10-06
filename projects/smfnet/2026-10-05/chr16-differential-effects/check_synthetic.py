"""Reproduce the changed-means/equal-correlations calibration experiment."""
import json
from pathlib import Path
import numpy as np
import torch
import bootstrap
import investigate

torch.set_num_threads(1)
positions=np.arange(24)*10
probs=np.linspace(.12,.4,12)


def draw(rng,n,flip=False):
    columns=[]
    for p in probs:
        q=p*p+.4*p*(1-p)
        state=rng.choice(4,n,p=[1-2*p+q,p-q,p-q,q])
        columns.extend([state//2,state%2])
    x=np.asarray(columns,dtype='int8').T
    return 1-x if flip else x


result={}
rng=np.random.default_rng(75)
permutation=[]
for k in range(30):
    a,b=draw(rng,400),draw(rng,400,True)
    permutation.append(investigate.compare(positions,a,positions,b,199,20,100+k,device='cpu'))
result['permutation_coc_rejections_30']=int(sum(r['coc_p']<.05 for r in permutation))
rng=np.random.default_rng(75)
multiplier=[]
for k in range(100):
    a,b=draw(rng,400),draw(rng,400,True)
    multiplier.append(bootstrap.compare(positions,a,positions,b,399,seed=k))
result['multiplier_coc_rejections_100']=sum(r['coc_p']<.05 for r in multiplier)
result['multiplier_profile_rejections_100']=sum(r['profile_p']<.05 for r in multiplier)
a=draw(rng,500);b=draw(rng,500)[:,rng.permutation(24)]
r=bootstrap.compare(positions,a,positions,b,999,seed=123)
result['changed_dependence_p']=r['coc_p']
assert r['coc_p']<=.01
i,j=np.triu_indices(24,1)
_,_,rho,psi=bootstrap.features(torch.tensor(a),torch.tensor(i),torch.tensor(j))
assert torch.max(psi.sum(0).abs())<1e-12
eps=1e-4;v=(a>=0).astype(float);x=np.maximum(a,0);w=np.ones(len(a));w[0]+=eps
n=w@v[:,i];mx=(w@(x[:,i]*v[:,j]))/n;my=(w@(x[:,j]*v[:,i]))/n;p11=(w@(x[:,i]*x[:,j]))/n
rr=(p11-mx*my)/np.sqrt(mx*(1-mx)*my*(1-my))
assert np.allclose((rr-rho.numpy())/eps,psi[0].numpy()/np.sqrt(len(a)/(len(a)-1)),atol=1e-8)
result['influence_finite_difference_check']='passed'
rng=np.random.default_rng(75)
shapes=[]
for k in range(100):
    a,b=draw(rng,400),draw(rng,400,True)
    shapes.append(bootstrap.compare(positions,a,positions,b,399,seed=k,shape=True))
result['shape_coc_rejections_100']=sum(r['coc_p']<.05 for r in shapes)
result['shape_profile_rejections_100']=sum(r['profile_p']<.05 for r in shapes)
result['description']='Independent correlated binary pairs: within-pair rho=0.4, across-pair rho=0. Group B complements all calls; population correlation matrices are unchanged but means differ.'
print(json.dumps(result,indent=2))
Path(__file__).with_name('synthetic_checks.json').write_text(json.dumps(result,indent=2)+'\n')
