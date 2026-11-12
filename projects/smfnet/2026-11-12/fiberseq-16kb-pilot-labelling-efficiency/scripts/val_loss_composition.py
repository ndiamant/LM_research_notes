"""How hard are the sites each model is scored on? Per-site entropy, centre 1 kb vs 16 kb.

For val regions (kept-window centres on chr8/chr20) and reads as the 16 kb run draws them
(covering >= 75% of the 16,384 bp SMF window), each A/T site's accessible fraction p across
the reads that call it, and its binary entropy H(p) in nats: the cross-entropy a model that
knew only each site's rate would score. Averaged over calls (as the loss is), in the central
1,024 bp (the 2 kb model's SMF window), the 15 kb around it, and the whole window; and for
the central 1,024 bp with reads drawn by the 2 kb (H5) rule instead.
"""

import os

import numpy as np

from smf_net.figwig_data import open_lanes, read_bed_regions, read_span

S = os.environ["SCRATCH"] + "/fiberseq/GM12878"
lanes = open_lanes(S + "/figwig")
rng = np.random.default_rng(0)
regions = [(c, s) for c, s in read_bed_regions(S + "/kept_windows.bed", 2114) if c in ("chr8", "chr20")]
regions = [regions[i] for i in rng.choice(len(regions), 400, replace=False)]


def site_stats(calls):
    """Per-site (calls, accessible) counts."""
    called = (calls != -1).sum(0)
    accessible = (calls == 0).sum(0)
    return called, accessible


def entropy(called, accessible):
    ok = called >= 5
    p = np.clip(accessible[ok] / called[ok], 1e-6, 1 - 1e-6)
    h = -(p * np.log(p) + (1 - p) * np.log(1 - p))
    return (h * called[ok]).sum(), called[ok].sum(), (accessible[ok]).sum()


totals = {k: np.zeros(3) for k in ("16kb: centre 1 kb", "16kb: outer 15 kb", "16kb: whole window", "2kb rule: centre 1 kb")}
for chrom, start in regions:
    centre = start + 1057
    l = lanes[chrom]
    a16 = centre - 8192
    if a16 < 0 or a16 + 16384 > l.chrom_length:
        continue
    idx = l.overlapping(*read_span(a16, 16384, 500, 0, 0.75))
    if len(idx) < 32:
        continue
    calls = l.calls(idx, a16, 16384).numpy()
    mid = slice(8192 - 512, 8192 + 512)
    c, a = site_stats(calls)
    inner = np.zeros(16384, bool)
    inner[mid] = True
    totals["16kb: centre 1 kb"] += entropy(c[inner], a[inner])
    totals["16kb: outer 15 kb"] += entropy(c[~inner], a[~inner])
    totals["16kb: whole window"] += entropy(c, a)
    # 2 kb rule: molecules overlapping the 2,048 bp input by 500 bp, scored on its middle half
    a2 = centre - 1024
    idx2 = l.overlapping(a2, a2 + 2048, 500)
    calls2 = l.calls(idx2, a2 + 512, 1024).numpy()
    totals["2kb rule: centre 1 kb"] += entropy(*site_stats(calls2))

print(f"{len(regions)} val regions; per-site entropy H(p) averaged over calls (nats)")
for k, (h, n, acc) in totals.items():
    print(f"  {k:22s} H {h / n:.4f}   accessible fraction of calls {acc / n:.3f}   calls {int(n):,}")
