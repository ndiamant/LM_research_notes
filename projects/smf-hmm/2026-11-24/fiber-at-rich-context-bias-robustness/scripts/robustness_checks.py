"""Robustness of model outputs to the two m6A labelling biases: molecule labelling (b_m) and sequence context.

1. **Molecule labelling.** The reads of each window are split into within-window b_m tertiles (reads with a
   finite b_m; windows with at least 6), so every tertile covers the same loci.
   - Per tertile: share of read positions in each element (open, TF, other_long, histone), for the posterior
     marginals (all windows) and MAP paths (the windows of `--segments`).
   - Reference: the kernel's mean P(bound) (model-free; it corrects for b_m), over read positions with an
     informative call within ±h bp.
   - Summary: ratio of the top to the bottom tertile per element.
2. **Sequence context.** Positions are binned into deciles of local context efficiency: the mean of
   sigmoid(beta + f(context)) over the A/T sites within ±`--context-half-width` bp (kernel reference table,
   no b_m; positions without an A/T site are dropped). Decile edges are over all window positions.
   - Per decile: element shares for the prior (weighted by the window's reads), posterior and MAP, and the
     kernel's mean P(bound).
   - Motif exclusion, as a sensitivity: all positions (`none`), outside the strong-footprint mask (`strong`,
     the TFs of `--strong-targets`), outside the 100-PWM mask (`all_pwms`: `--motifs` plus the benchmark
     MEME files); hits at `--motif-fpr`, widened by `--mask-pad` bp.
   - Summary: ratio of the highest- to the lowest-efficiency decile per element.

Intervals: 95% bootstrap over windows (the same resampled windows for every bin of a row, so ratios are paired).

Outputs in OUT_DIR: `bm_tertiles.tsv`, `context_deciles.tsv`, `summary.tsv` (the ratios).

Usage:
  python -m evaluation.tf_chip.robustness_checks OUT_DIR DATASET_DIR [DATASET_DIR ...] \\
      --model NAME[:LABEL] [--model ...] --efficiency-path EFF.h5 [--segments NAME:MAP=viterbi_segments.tsv.gz ...]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from evaluation.tf_chip.analyze_covariation import read_table
from evaluation.tf_chip.describe_model_outputs import DEFAULT_TARGETS, group_index
from evaluation.tf_chip.kernel_pbound import call_llr, kernel_pbound, kmer_codes, load_reference
from evaluation.tf_chip.plot_model_comparison import DEFAULT_KERNEL_REFERENCE, GROUPS, motif_scanners
from evaluation.tf_chip.specificity_metrics import BENCHMARK_PWMS, DEFAULT_MOTIFS, REFERENCES, meme_scanner, motif_mask

EXCLUSIONS = ["none", "strong", "all_pwms"]
WIDTH = 1 + len(GROUPS)  # positions, then one sum per group


def raw_efficiencies(efficiency_file: h5py.File, molecule_ids) -> np.ndarray:
    """b_m per read, NaN where the molecule has no estimate."""
    out = []
    for mid in molecule_ids:
        region, row = str(mid).rsplit(":", 1)
        out.append(efficiency_file["efficiency"][region][int(row)])
    return np.asarray(out, np.float64)


def within_window_tertiles(b_m: np.ndarray, min_reads: int = 6) -> np.ndarray:
    """Tertile 0-2 of each read by b_m rank within the window; -1 for NaN b_m or too few reads."""
    out = np.full(len(b_m), -1)
    ok = np.flatnonzero(np.isfinite(b_m))
    if len(ok) < min_reads:
        return out
    rank = np.empty(len(ok)); rank[np.argsort(b_m[ok], kind="stable")] = np.arange(len(ok))
    out[ok] = np.minimum((3 * rank / len(ok)).astype(int), 2)
    return out


def context_efficiency(dna: np.ndarray, beta: float, f: np.ndarray, half_width: int) -> np.ndarray:
    """Mean expected m6A rate sigmoid(beta + f(context)) over A/T sites within ±half_width; NaN without one."""
    codes = kmer_codes(dna)
    has = codes >= 0
    e = np.where(has, 1.0 / (1.0 + np.exp(-(beta + f[np.maximum(codes, 0)]))), 0.0)
    cs, cn = np.r_[0, np.cumsum(e)], np.r_[0, np.cumsum(has)]
    idx = np.arange(len(dna))
    lo, hi = np.clip(idx - half_width, 0, len(dna)), np.clip(idx + half_width + 1, 0, len(dna))
    n = cn[hi] - cn[lo]
    return np.where(n > 0, (cs[hi] - cs[lo]) / np.maximum(n, 1), np.nan)


def map_groups(segments: pd.DataFrame, n_reads: int, length: int) -> np.ndarray:
    """(reads, L) GROUPS index of one window's MAP paths from its segment rows (-1 where uncovered)."""
    out = np.full((n_reads, length), -1, np.int8)
    for r, st, s, n in zip(segments.read.to_numpy(), segments.state.to_numpy(), segments.start.to_numpy(),
                           segments.length.to_numpy()):
        out[r, s:s + n] = GROUPS.index(st)
    return out


def binned_sums(shares: np.ndarray, weight: np.ndarray, bins: np.ndarray, n_bins: int) -> np.ndarray:
    """(n_bins, WIDTH): positions (sum of weight) and per-group sums of `shares` (L, 4) weighted, per bin."""
    keep = bins >= 0
    out = np.zeros((n_bins, WIDTH))
    out[:, 0] = np.bincount(bins[keep], weight[keep], minlength=n_bins)
    for g in range(len(GROUPS)):
        out[:, 1 + g] = np.bincount(bins[keep], shares[keep, g], minlength=n_bins)
    return out


def summarize(acc: dict, n_bins: int, n_boot: int, seed: int):
    """Share tables with window-bootstrap intervals, and the ratio of the last to the first bin."""
    rng = np.random.default_rng(seed)
    rows, ratios = [], []
    for key, a in acc.items():  # a: (windows, n_bins, WIDTH)
        used = np.flatnonzero(a[:, :, 0].sum(1) > 0)
        a = a[used]
        boot = rng.integers(0, len(a), (n_boot, len(a)))
        tot, bt = a.sum(0), a[boot].sum(1)  # (n_bins, WIDTH), (n_boot, n_bins, WIDTH)
        labels = ["P(bound)"] if key[2] == "kernel" else GROUPS
        for j, g in enumerate(labels):
            share = tot[:, 1 + j] / tot[:, 0]
            bs = bt[:, :, 1 + j] / bt[:, :, 0]
            for b in range(n_bins):
                rows.append(dict(zip(("check", "model", "level", "exclusion"), key), bin=b, element=g, share=share[b],
                                 lo=np.nanpercentile(bs[:, b], 2.5), hi=np.nanpercentile(bs[:, b], 97.5),
                                 positions=tot[b, 0], windows=len(a)))
            r, br = share[-1] / share[0], bs[:, -1] / bs[:, 0]
            ratios.append(dict(zip(("check", "model", "level", "exclusion"), key), element=g, ratio_last_first=r,
                               lo=np.nanpercentile(br, 2.5), hi=np.nanpercentile(br, 97.5)))
    return pd.DataFrame(rows), pd.DataFrame(ratios)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("dataset_dirs", type=Path, nargs="+")
    parser.add_argument("--model", action="append", required=True, help="NAME[:LABEL] (posteriors/NAME.h5, NAME_prior.h5)")
    parser.add_argument("--segments", action="append", default=[], help="NAME:MAP=segments.tsv.gz")
    parser.add_argument("--efficiency-path", type=Path, required=True)
    parser.add_argument("--kernel-reference", type=Path, default=DEFAULT_KERNEL_REFERENCE)
    parser.add_argument("--kernel-eps", type=float, default=0.02)
    parser.add_argument("--kernel-half-width", type=int, default=7)
    parser.add_argument("--context-half-width", type=int, default=10)
    parser.add_argument("--n-context-bins", type=int, default=10)
    parser.add_argument("--motifs", type=Path, default=DEFAULT_MOTIFS)
    parser.add_argument("--strong-targets", type=Path, default=DEFAULT_TARGETS)
    parser.add_argument("--motif-fpr", type=float, default=1e-6)
    parser.add_argument("--mask-pad", type=int, default=10)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    models = [m.split(":", 1) if ":" in m else (m, m) for m in args.model]
    labels = dict(models)
    segments = {}
    for spec in args.segments:
        key, path = spec.split("=", 1)
        name, level = key.split(":", 1)
        t = pd.read_csv(path, sep="\t")
        segments[(name, level)] = {k: g for k, g in t.groupby(["chrom", "window_id"])}
    beta, f = load_reference(args.kernel_reference)
    strong = motif_scanners(args.motifs, args.motif_fpr, set(args.strong_targets.read_text().split()))
    every = motif_scanners(args.motifs, args.motif_fpr) + [meme_scanner(REFERENCES / f"{m}.meme", args.motif_fpr)
                                                          for m in BENCHMARK_PWMS]
    efficiency_file = h5py.File(args.efficiency_path)

    # pass 1: context efficiency and motif masks per window; decile edges over all positions
    windows = []
    for d in args.dataset_dirs:
        with h5py.File(d / "dataset.h5") as fh:
            meta = read_table(fh["window_metadata"]).set_index("window_id")
            for wid in sorted(fh["windows"]):
                dna = fh["windows"][wid]["dna"][:].astype(np.int64)
                windows.append(dict(dir=d, wid=wid, chrom=str(meta.chrom[wid]),
                                    ctx=context_efficiency(dna, beta, f, args.context_half_width),
                                    masks={"none": np.zeros(len(dna), bool), "strong": motif_mask(dna, strong, args.mask_pad),
                                           "all_pwms": motif_mask(dna, every, args.mask_pad)}))
    ctx_all = np.concatenate([w["ctx"] for w in windows])
    edges = np.nanquantile(ctx_all, np.linspace(0, 1, args.n_context_bins + 1)[1:-1])
    print(f"{len(windows)} windows; context-efficiency decile edges {np.round(edges, 3)}", flush=True)

    nb = args.n_context_bins
    acc: dict[tuple, np.ndarray] = {}

    def add(key, wi, n_bins, values):
        acc.setdefault(key, np.zeros((len(windows), n_bins, WIDTH)))[wi] += values

    by_dir: dict[Path, list[int]] = {}
    for wi, w in enumerate(windows):
        by_dir.setdefault(w["dir"], []).append(wi)
    for d, wis in by_dir.items():
        files = {(m, lvl): h5py.File(d / "posteriors" / (f"{m}_prior.h5" if lvl == "prior" else f"{m}.h5"))
                 for m, _ in models for lvl in ("prior", "posterior")}
        cols = {k: group_index(fh) for k, fh in files.items()}
        offsets = {k: dict(zip(fh["window_ids"].asstr()[:], zip(fh["window_offsets"][:-1], fh["window_offsets"][1:])))
                   for k, fh in files.items()}
        with h5py.File(d / "dataset.h5") as fh:
            for n_done, wi in enumerate(wis):
                w = windows[wi]
                win = fh["windows"][w["wid"]]
                dna = win["dna"][:].astype(np.int64); L = len(dna)
                calls = np.full((len(win["observations"]), L), -1, np.int8)
                calls[:, win["positions"][:]] = win["observations"][:]
                R = len(calls)
                b_m = raw_efficiencies(efficiency_file, win["molecule_id"].asstr()[:])
                tert = within_window_tertiles(b_m)
                k = kernel_pbound(call_llr(calls, dna, np.nan_to_num(b_m), beta, f, args.kernel_eps), args.kernel_half_width)
                kin = np.isfinite(k)
                ctx_bin = np.where(np.isfinite(w["ctx"]), np.searchsorted(edges, w["ctx"], side="right"), -1)
                # kernel reference (model independent)
                for e in EXCLUSIONS:
                    bins = np.where(w["masks"][e], -1, ctx_bin)
                    kv = np.zeros((nb, WIDTH))
                    kv[:, 0] = np.bincount(bins[bins >= 0], kin.sum(0)[bins >= 0], minlength=nb)
                    kv[:, 1] = np.bincount(bins[bins >= 0], np.where(kin, k, 0).sum(0)[bins >= 0], minlength=nb)
                    add(("context", "kernel", "kernel", e), wi, nb, kv)
                tv = np.zeros((3, WIDTH))
                for t in range(3):
                    sel = tert == t
                    tv[t, 0], tv[t, 1] = kin[sel].sum(), np.where(kin[sel], k[sel], 0).sum()
                if (tert >= 0).any():
                    add(("bm", "kernel", "kernel", "none"), wi, 3, tv)
                for m, label in models:
                    levels = {}
                    for lvl in ("prior", "posterior"):
                        if w["wid"] not in offsets[(m, lvl)]:
                            continue
                        a, b = offsets[(m, lvl)][w["wid"]]
                        post = files[(m, lvl)]["posterior"][a:b]
                        levels[lvl] = np.stack([post[..., c].sum(-1) for c in cols[(m, lvl)]], -1)  # (R, L, 4)
                    seg = segments.get((m, "MAP"), {}).get((w["chrom"], w["wid"]))
                    if seg is not None:
                        z = map_groups(seg, R, L)
                        levels["MAP"] = np.stack([(z == g) for g in range(len(GROUPS))], -1).astype(np.float64)
                    for lvl, g in levels.items():
                        per_pos = g.sum(0) if lvl != "prior" else g[0] * R  # prior rows are identical across reads
                        for e in EXCLUSIONS:
                            bins = np.where(w["masks"][e], -1, ctx_bin)
                            add(("context", label, lvl, e), wi, nb, binned_sums(per_pos, np.full(L, float(R)), bins, nb))
                        if lvl != "prior" and (tert >= 0).any():
                            tv = np.zeros((3, WIDTH))
                            for t in range(3):
                                sel = tert == t
                                tv[t, 0], tv[t, 1:] = sel.sum() * L, g[sel].sum((0, 1))
                            add(("bm", label, lvl, "none"), wi, 3, tv)
                if n_done % 500 == 0:
                    print(f"{d.name}: window {n_done}/{len(wis)}", flush=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    bm_rows, bm_ratio = summarize({k: v for k, v in acc.items() if k[0] == "bm"}, 3, args.bootstrap, args.seed)
    cx_rows, cx_ratio = summarize({k: v for k, v in acc.items() if k[0] == "context"}, nb, args.bootstrap, args.seed)
    cx_rows["edge_low"] = cx_rows.bin.map(dict(enumerate(np.r_[np.nanmin(ctx_all), edges])))
    bm_rows.to_csv(args.out_dir / "bm_tertiles.tsv", sep="\t", index=False)
    cx_rows.to_csv(args.out_dir / "context_deciles.tsv", sep="\t", index=False)
    summary = pd.concat([bm_ratio, cx_ratio], ignore_index=True)
    summary.to_csv(args.out_dir / "summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(bm_rows.pivot_table(index=["model", "level", "element"], columns="bin", values="share").round(4).to_string())
    print(summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
