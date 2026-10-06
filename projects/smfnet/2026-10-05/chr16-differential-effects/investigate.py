"""Exploratory, coverage-stratified whole-molecule permutation analysis.

This is an observed-data feasibility experiment, not model evaluation.
Outputs belong on scratch; the script and methods are kept with the code.
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

CELL_TYPES = ('NP', 'MEL', 'ES', 'C2C12_fixed')
WORKER_FILES = {}


def initialize_worker(paths):
    torch.set_num_threads(1)
    WORKER_FILES.update({c:h5py.File(p, 'r') for c,p in paths.items()})


def process_region(task):
    region, permutations, coverage, controls, device = task
    data = {c:load_region(f, region) for c,f in WORKER_FILES.items()}
    pairs = [(c,c) for c in CELL_TYPES] if controls else list(itertools.combinations(CELL_TYPES, 2))
    results = []
    for ca, cb in pairs:
        seed = int.from_bytes(hashlib.sha256(f'{region}/{ca}/{cb}/20261005'.encode()).digest()[:4], 'little')
        pa, a = data[ca]
        pb, b = data[cb]
        if controls:
            idx = np.random.default_rng(seed).permutation(len(a))
            a, b = a[idx[:len(idx)//2]], a[idx[len(idx)//2:]]
        row = dict(region=region, pair=ca+'-'+cb)
        row.update(compare(pa, a, pb, b, permutations, coverage, seed, device=device))
        results.append(row)
    return results


def pearson_rows(a, b, minimum=10):
    keep = torch.isfinite(a) & torch.isfinite(b)
    n = keep.sum(-1)
    a = torch.where(keep, a, 0)
    b = torch.where(keep, b, 0)
    ac = torch.where(keep, a - a.sum(-1, keepdim=True) / n.clamp_min(1)[:, None], 0)
    bc = torch.where(keep, b - b.sum(-1, keepdim=True) / n.clamp_min(1)[:, None], 0)
    denom = (ac.square().sum(-1) * bc.square().sum(-1)).sqrt()
    r = (ac * bc).sum(-1) / denom.clamp_min(1e-30)
    r = torch.where((n >= minimum) & (denom > 1e-15), r.clamp(-1, 1), torch.nan)
    return r, n


def pair_r(n, sx, sy, sxy):
    scale = ((n * sx - sx.square()) * (n * sy - sy.square())).clamp_min(0).sqrt()
    return torch.where((n > 0) & (scale > 0), (n * sxy - sx * sy) / scale.clamp_min(1e-30), torch.nan)


def load_region(handle, region):
    pos = handle['smf_pos'][region][:].astype(np.int64)
    x = handle['smf_mat'][region][:]
    ids = handle['read_id'][region][:]
    if len(np.unique(ids)) != len(ids):
        raise ValueError(f'Duplicate molecules within {region}')
    if x.shape != (len(ids), len(pos)) or len(np.unique(pos)) != len(pos):
        raise ValueError(f'Invalid dimensions/positions: {region}')
    if not np.isin(x, [-1, 0, 1]).all():
        raise ValueError(f'Nonbinary calls: {region}')
    keep = (pos >= 512) & (pos < 1536)
    return pos[keep], x[:, keep].astype(np.int8)


def compare(pa, a, pb, b, permutations, coverage, seed, device='cuda', batch_size=128):
    start = time.monotonic()
    pos, ia, ib = np.intersect1d(pa, pb, return_indices=True)
    a, b = a[:, ia], b[:, ib]
    keep = ((a != -1).sum(0) >= coverage) & ((b != -1).sum(0) >= coverage)
    pos, a, b = pos[keep], a[:, keep], b[:, keep]
    a, b = a[(a != -1).any(1)], b[(b != -1).any(1)]
    out = dict(n_sites=len(pos), reads_a=len(a), reads_b=len(b))
    if len(pos) < 10:
        return dict(out, status='fewer_than_10_sites')
    pooled = np.concatenate([a, b])
    valid_np = pooled != -1
    _, strata = np.unique(np.packbits(valid_np, axis=1), axis=0, return_inverse=True)
    na = np.bincount(strata[:len(a)], minlength=strata.max()+1)
    nb = np.bincount(strata[len(a):], minlength=len(na))
    exchangeable = (na > 0) & (nb > 0)
    out.update(exchangeable_read_fraction=float(exchangeable[strata].mean()),
               n_shared_masks=int(exchangeable.sum()))
    torch.manual_seed(seed)
    gen = torch.Generator(device=device).manual_seed(seed)
    v = torch.as_tensor(valid_np, device=device, dtype=torch.float32)
    x = torch.as_tensor(np.maximum(pooled, 0), device=device, dtype=torch.float32)
    rows, cols = np.triu_indices(len(pos), 1)
    dist = pos[cols] - pos[rows]
    rows, cols = rows[(dist > 0) & (dist <= 500)], cols[(dist > 0) & (dist <= 500)]
    ri, ci = torch.as_tensor(rows, device=device), torch.as_tensor(cols, device=device)
    joint = v[:, ri] * v[:, ci]
    pair_keep = (joint[:len(a)].sum(0) >= coverage) & (joint[len(a):].sum(0) >= coverage)
    ri, ci, joint = ri[pair_keep], ci[pair_keep], joint[:, pair_keep]
    s, p = len(pos), len(ri)
    features = torch.cat([v, x, joint, x[:, ri]*v[:, ci], v[:, ri]*x[:, ci], x[:, ri]*x[:, ci]], dim=1)
    total = features.sum(0, dtype=torch.float64)
    labels = torch.zeros(len(pooled), device=device)
    labels[:len(a)] = 1
    order0 = np.argsort(strata, kind='stable')
    template = labels[torch.as_tensor(order0, device=device)]
    strata_t = torch.as_tensor(strata, device=device, dtype=torch.float64)
    all_profile, all_coc = [], []
    observed_means = observed_pairs = None
    observed_n_pairs = 0
    for offset in range(0, permutations+1, batch_size):
        count = min(batch_size, permutations+1-offset)
        rand = torch.rand((count, len(pooled)), device=device, dtype=torch.float64, generator=gen)
        order = torch.argsort(rand + strata_t, dim=1)
        weights = torch.zeros((count, len(pooled)), device=device)
        weights.scatter_(1, order, template.expand(count, -1))
        if offset == 0:
            weights[0] = labels
        sums_a = (weights @ features).double()
        sums_b = total - sums_a
        ma, mb = sums_a[:, s:2*s] / sums_a[:, :s], sums_b[:, s:2*s] / sums_b[:, :s]
        pr, _ = pearson_rows(ma, mb)
        aa = sums_a[:, 2*s:].reshape(count, 4, p)
        bb = sums_b[:, 2*s:].reshape(count, 4, p)
        ra, rb = pair_r(*aa.unbind(1)), pair_r(*bb.unbind(1))
        cr, npair = pearson_rows(ra, rb)
        if offset == 0:
            observed_means = (ma[0].cpu().numpy(), mb[0].cpu().numpy())
            observed_pairs = (ra[0].cpu().numpy(), rb[0].cpu().numpy())
            observed_n_pairs = int(npair[0])
        all_profile.extend(pr.cpu().tolist())
        all_coc.extend(cr.cpu().tolist())
    out.update(n_coverage_pairs=p, n_defined_pairs=observed_n_pairs)
    for name, vals, vectors in [('profile', all_profile, observed_means), ('coc', all_coc, observed_pairs)]:
        vals = np.array(vals)
        out[f'{name}_r'] = vals[0]
        out[f'{name}_null_median_r'] = float(np.nanmedian(vals[1:])) if np.isfinite(vals[1:]).any() else np.nan
        out[f'{name}_invalid_null_fraction'] = float((~np.isfinite(vals[1:])).mean())
        # An undefined null correlation is maximally dissimilar: conservative.
        null = np.nan_to_num(vals[1:], nan=-1.0)
        out[f'{name}_p'] = (1 + np.count_nonzero(null <= vals[0] + 1e-10)) / (permutations+1) if np.isfinite(vals[0]) else np.nan
        d = vectors[0]-vectors[1]
        out[f'{name}_rms_difference'] = float(np.sqrt(np.nanmean(d*d))) if np.isfinite(d).any() else np.nan
        out[f'{name}_mean_difference'] = float(np.nanmean(d)) if np.isfinite(d).any() else np.nan
    out.update(status='ok', seconds=time.monotonic()-start)
    return out


def adjust(p, by=False):
    p = np.nan_to_num(np.asarray(p, float), nan=1.0)
    n = len(p)
    order = np.argsort(p)
    factor = np.sum(1 / np.arange(1, n+1)) if by else 1
    vals = np.minimum.accumulate((p[order]*n*factor/np.arange(1,n+1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(vals, 1)
    return out


def summarize(frame, output):
    summary = []
    for metric in ['profile', 'coc']:
        pcol = metric+'_p'
        frame[metric+'_q_global'] = adjust(frame[pcol])
        frame[metric+'_q_global_by'] = adjust(frame[pcol], by=True)
        for pair, ix in frame.groupby('pair').groups.items():
            frame.loc[ix, metric+'_q_pair'] = adjust(frame.loc[ix, pcol])
            f = frame.loc[ix]
            summary.append(dict(pair=pair, metric=metric, regions=len(f), testable=int(f[pcol].notna().sum()),
                nominal_p05=int((f[pcol]<.05).sum()), bh_pair_05=int((f[metric+'_q_pair']<=.05).sum()),
                bh_global_05=int((f[metric+'_q_global']<=.05).sum()), by_global_05=int((f[metric+'_q_global_by']<=.05).sum()),
                median_r=f[metric+'_r'].median(), median_null_r=f[metric+'_null_median_r'].median(),
                median_exchangeable_fraction=f['exchangeable_read_fraction'].median()))
    frame.to_csv(output/'regions.tsv', sep='\t', index=False)
    summary = pd.DataFrame(summary)
    summary.to_csv(output/'summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--limit', type=int)
    ap.add_argument('--permutations', type=int, default=999)
    ap.add_argument('--coverage', type=int, default=20)
    ap.add_argument('--controls', action='store_true')
    ap.add_argument('--device', default='cuda')
    ap.add_argument('--workers', type=int, default=1)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(4)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.set_float32_matmul_precision('highest')
    files, paths, keys = {}, {}, {}
    for cell in CELL_TYPES:
        cfg = yaml.safe_load((Path(os.environ['SCRATCH'])/'smf_models'/cell/'.hydra/config.yaml').read_text())
        paths[cell] = cfg['data']['cfg']['h5_path']
        files[cell] = h5py.File(paths[cell], 'r')
        keys[cell] = {k for k in files[cell]['smf_mat'] if k.startswith('chr16:')}
    regions = sorted(set.intersection(*keys.values()))
    np.random.default_rng(20261005).shuffle(regions)
    metadata = dict(chromosome='chr16', window=[512,1536], permutations=args.permutations,
        min_site_and_joint_coverage=args.coverage, min_sites=10, min_pairs=10, max_pair_distance=500,
        common_regions=len(regions), regions_per_cell={c:len(k) for c,k in keys.items()},
        paths=paths, controls=args.controls, seed=20261005,
        torch_version=torch.__version__, numpy_version=np.__version__, device=torch.cuda.get_device_name() if args.device == 'cuda' else 'cpu',
        null='equal call-vector distributions conditional on exact coverage mask',
        multiple_testing='BH separately per metric, across all region-pair tests; untestable p=1; BY sensitivity')
    if args.limit:
        regions = regions[:args.limit]
    metadata['analyzed_regions'] = len(regions)
    (args.output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    results = []
    started = time.monotonic()
    tasks = [(r,args.permutations,args.coverage,args.controls,args.device) for r in regions]
    if args.workers > 1:
        if args.device != 'cpu':
            raise ValueError('Multiple workers are only supported on CPU')
        executor = ProcessPoolExecutor(max_workers=args.workers, initializer=initialize_worker,
            initargs=(paths,), mp_context=multiprocessing.get_context('spawn'))
        batches = executor.map(process_region, tasks, chunksize=4)
    else:
        initialize_worker(paths)
        batches = map(process_region, tasks)
    with (args.output/'raw.jsonl').open('w') as stream:
        for i, batch in enumerate(batches):
            for row in batch:
                results.append(row)
                stream.write(json.dumps(row)+'\n')
            stream.flush()
            if (i+1)%25 == 0 or i+1 == len(regions):
                print(f'{i+1}/{len(regions)} regions; {time.monotonic()-started:.1f}s', flush=True)
    if args.workers > 1:
        executor.shutdown()
    summarize(pd.DataFrame(results), args.output)
    for f in files.values():
        f.close()


if __name__ == '__main__':
    main()
