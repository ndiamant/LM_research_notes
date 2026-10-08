"""What ends observed SMF runs, and how often bridging a missing gap would be wrong.

Reproduces the tables in ../run-length-missingness-gap-bridging.md. Two input
modes:

  samples  an SMFNet samples.h5 (observed reads aligned to model positions)
  gch      a raw GCH short-read .h5 with smf_mat/ and smf_pos/ groups

Usage (SMFNet group env, CPU only):

  CUDA_VISIBLE_DEVICES="" $GROUP_HOME/uv-envs/SMFNet-cu124/bin/python \
      analyze_run_boundaries.py samples PATH/samples.h5 --regions 300
  CUDA_VISIBLE_DEVICES="" $GROUP_HOME/uv-envs/SMFNet-cu124/bin/python \
      analyze_run_boundaries.py gch PATH/file.h5 --regions 1500 --reads-per-region 80

For every maximal run of one value inside a read, each side is classified as a
real call of the other value, the window edge, a read end (missing up to the
read's first/last call), or an internal missing gap. Runs touching the window
edge are excluded, matching smf_net.evaluation.run_lengths.run_length_bp.

The bridging error rate is simulated: take fully called stretches of k+2
consecutive sites whose two end calls agree, hide the k interior calls, and
record whether any hidden call had the other value. This is the probability
that filling a gap of that span with the flanking value would erase a real
opposite run.
"""

import argparse
import itertools

import numpy as np
import pandas as pd

SPAN_BINS = [0, 10, 20, 30, 40, 60, 80, 120, 10_000]


def _samples_regions(path, limit):
    from smf_net.evaluation.metrics import align_sample_region
    from smf_net.evaluation.sample_io import read_samples_h5

    for region in itertools.islice(read_samples_h5(path), limit):
        aligned = align_sample_region(region)
        yield aligned.observed_smf.numpy(), aligned.positions.numpy()


def _gch_regions(path, limit, region_length=2114, seed=0):
    import h5py

    with h5py.File(path, "r") as h5:
        keys = list(h5["smf_mat"].keys())
        rng = np.random.default_rng(seed)
        for key in rng.choice(keys, min(limit, len(keys)), replace=False):
            reads = h5["smf_mat"][key][:].astype(np.int8)
            pos = h5["smf_pos"][key][:]
            keep = (pos >= 0) & (pos < region_length)
            yield reads[:, keep], pos[keep]


def analyze(regions, reads_per_region=None, seed=0):
    rng = np.random.default_rng(seed)
    runs, gaps, hidden, spacing, spans = [], [], [], [], []
    for reads, pos in regions:
        n = reads.shape[1]
        if n < 3:
            continue
        spacing.append(np.median(np.diff(pos)))
        for read in reads[:reads_per_region]:
            called = np.nonzero(read != -1)[0]
            if called.size < 2:
                continue
            first, last = called[0], called[-1]
            values, called_pos = read[called], pos[called]
            spans.append(called_pos[-1] - called_pos[0])
            for j in np.nonzero(np.diff(called) > 1)[0]:
                gaps.append(
                    (
                        called_pos[j + 1] - called_pos[j],
                        called[j + 1] - called[j] - 1,
                        values[j] == values[j + 1],
                    )
                )

            def side(i, left):
                if i < 0 or i >= n:
                    return "edge"
                if read[i] != -1:
                    return "call"
                return "read end" if (i < first if left else i > last) else "gap"

            j = 0
            while j < n:
                if read[j] == -1:
                    j += 1
                    continue
                k = j
                while k + 1 < n and read[k + 1] == read[j]:
                    k += 1
                runs.append((side(j - 1, True), side(k + 1, False)))
                j = k + 1

            for hide in (1, 2, 3, 4, 6):
                if called.size < hide + 2:
                    break
                start = rng.integers(0, called.size - hide - 1)
                window = called[start : start + hide + 2]
                if window[-1] - window[0] != hide + 1:
                    continue  # not consecutive sites
                stretch = values[start : start + hide + 2]
                if stretch[0] != stretch[-1]:
                    continue
                hidden.append(
                    (
                        called_pos[start + hide + 1] - called_pos[start],
                        stretch[0],
                        bool((stretch[1:-1] != stretch[0]).any()),
                    )
                )

    runs = pd.DataFrame(runs, columns=["left", "right"])
    runs = runs[(runs.left != "edge") & (runs.right != "edge")]
    bound = np.where(
        (runs.left == "call") & (runs.right == "call"),
        "complete",
        np.where(
            (runs.left == "read end") | (runs.right == "read end"),
            "read end",
            "internal gap only",
        ),
    )
    gaps = pd.DataFrame(gaps, columns=["bp", "sites", "flanks_agree"])
    hidden = pd.DataFrame(hidden, columns=["bp", "value", "hidden_opposite"])
    print(f"median site spacing {np.median(spacing)} bp, median read span {np.median(spans):.0f} bp")
    print("interior runs bounded by:", pd.Series(bound).value_counts(normalize=True).round(3).to_dict())
    print(
        "internal gaps: per read %.2f, single-site %.2f, median span %.0f bp, flanks agree %.2f"
        % (len(gaps) / len(spans), (gaps.sites == 1).mean(), gaps.bp.median(), gaps.flanks_agree.mean())
    )
    gaps["bin"] = pd.cut(gaps.bp, SPAN_BINS)
    hidden["bin"] = pd.cut(hidden.bp, SPAN_BINS)
    table = pd.DataFrame(
        {
            "share of real gaps": gaps.bin.value_counts(normalize=True).sort_index(),
            "bridge wrong, zeros": hidden[hidden.value == 0].groupby("bin", observed=True).hidden_opposite.mean(),
            "bridge wrong, ones": hidden[hidden.value == 1].groupby("bin", observed=True).hidden_opposite.mean(),
        }
    )
    print(table.round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["samples", "gch"])
    parser.add_argument("path")
    parser.add_argument("--regions", type=int, default=300)
    parser.add_argument("--reads-per-region", type=int, default=None)
    args = parser.parse_args()
    source = _samples_regions if args.mode == "samples" else _gch_regions
    analyze(source(args.path, args.regions), args.reads_per_region)
