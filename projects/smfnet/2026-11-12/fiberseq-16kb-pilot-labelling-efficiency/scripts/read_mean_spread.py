"""Do generated reads vary as much between molecules as observed ones?

Per region, the variance across reads of each read's mean accessibility (accessible
fraction of its called sites). Generated reads are given the call pattern (no-call
positions) of a randomly paired observed read first, so both sides carry the same
sampling noise from missing calls. Reports the median over regions of var(generated) /
var(observed), and the pooled ratio.
"""

import sys

import numpy as np

from smf_net.evaluation.sample_io import read_samples_h5

rng = np.random.default_rng(0)
for label, path in zip(sys.argv[1::2], sys.argv[2::2]):
    ratios, v_obs, v_gen = [], [], []
    for region in read_samples_h5(path):
        L = len(region.dna)
        obs = np.full((region.observed_smf.shape[0], L), -1, np.int8)
        obs[:, region.observed_positions.numpy()] = region.observed_smf.numpy()
        gen = np.full((region.predicted_smf.shape[0], L), -1, np.int8)
        gen[:, region.predicted_positions.numpy()] = region.predicted_smf.numpy()
        pattern = obs[rng.integers(0, len(obs), len(gen))] != -1
        gen = np.where(pattern, gen, -1)

        def read_means(x):
            n = (x != -1).sum(1)
            keep = n >= 50
            return ((x == 0).sum(1)[keep] / n[keep])

        mo, mg = read_means(obs), read_means(gen)
        if len(mo) < 20 or len(mg) < 20:
            continue
        v_obs.append(mo.var())
        v_gen.append(mg.var())
        ratios.append(mg.var() / mo.var())
    ratios = np.array(ratios)
    print(f"{label:28s} {len(ratios):5d} regions: var(read means) generated / observed: "
          f"median {np.median(ratios):.2f} (IQR {np.percentile(ratios, 25):.2f}-{np.percentile(ratios, 75):.2f}), "
          f"pooled {np.sum(v_gen) / np.sum(v_obs):.2f}")
