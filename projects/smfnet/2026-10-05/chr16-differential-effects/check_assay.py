"""Verify every analyzed HDF5 site is eligible under the current model mask."""
import json
import os
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from pyfaidx import Fasta
import torch
import yaml
from smf_net.data import dna_str_to_tensor, get_valid_mask, parse_region_string

torch.set_num_threads(1)
root=Path(os.environ['SCRATCH'])/'differential_chr16_20261005'
metadata=json.loads((root/'shape_full/metadata.json').read_text())
regions=pd.read_csv(root/'shape_full/regions.tsv',sep='\t',usecols=['region']).region.unique()
cfg=yaml.safe_load((Path(os.environ['SCRATCH'])/'smf_models/NP/.hydra/config.yaml').read_text())
fasta=Fasta(cfg['data']['cfg']['fasta_path'],rebuild=False)
files={c:h5py.File(p,'r') for c,p in metadata['paths'].items()}
counts={c:dict(sites=0,ineligible_sites=0,affected_regions=0) for c in files}
for region in regions:
    chrom,start,end=parse_region_string(region)
    dna=dna_str_to_tensor(fasta[chrom][start-1:end].seq.upper())
    assert len(dna)==2048
    mask=get_valid_mask(dna,'gch_accessibility').numpy()
    for cell,f in files.items():
        pos=f['smf_pos'][region][:].astype(np.int64)
        pos=pos[(pos>=512)&(pos<1536)]
        n=int((~mask[pos]).sum())
        counts[cell]['sites']+=len(pos)
        counts[cell]['ineligible_sites']+=n
        counts[cell]['affected_regions']+=int(n>0)
result=dict(regions=len(regions),assay='gch_accessibility',rule='DGCH',window=[512,1536],counts=counts)
print(json.dumps(result,indent=2))
Path(__file__).with_name('assay_mask_check.json').write_text(json.dumps(result,indent=2)+'\n')
for f in files.values():f.close()
assert all(x['ineligible_sites']==0 for x in counts.values())
