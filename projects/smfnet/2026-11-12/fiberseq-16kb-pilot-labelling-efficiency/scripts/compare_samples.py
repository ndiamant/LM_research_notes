"""Observed vs generated reads in a samples.h5 file (scripts/sample.py on a data split).

Writes to --out:
- example_<region>.png for a few regions: observed and generated reads (binned accessible
  fraction, one row per read, sorted by each read's mean), their mean profiles, and their
  bin-by-bin co-accessibility (correlation across reads) matrices. Generated reads are
  subsampled to the region's observed read count, so both matrices carry the same noise.
- summary.png and summary.json: co-accessibility by distance over all regions (raw and
  adjusted for each read's labelling level), observed vs generated; with --elements (candidates.bed),
  the correlation of element pairs >= 2 kb apart (CTCF-bound pairs and other active
  pairs), observed vs generated.

Usage: compare_samples.py --samples samples.h5 --out dir [--elements candidates.bed]
"""

import argparse
import json
import os
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter

from smf_net.data import parse_region_string
from smf_net.evaluation.sample_io import read_samples_h5

ap = argparse.ArgumentParser()
ap.add_argument("--samples", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--elements", default=None, help="candidates.bed, for element pairs")
ap.add_argument("--n-examples", type=int, default=6)
ap.add_argument("--regions", nargs="*", default=None,
                help="only these region IDs (examples and summary), e.g. to redraw one")
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--no-summary", action="store_true",
                help="draw the examples only, leaving summary.png/json in --out as they are")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)
rng = np.random.default_rng(args.seed)

SURFACE, INK, INK2, MUTED, GRID, AXIS = ("#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7")
OBSERVED, GENERATED = "#2a78d6", "#eb6834"  # categorical slots 1-2
ACCESS = LinearSegmentedColormap.from_list(  # sequential blue ramp, 100 -> 700
    "access", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
ACCESS.set_bad(SURFACE)
CORR = LinearSegmentedColormap.from_list("corr", ["#1c5cab", "#f0efec", "#e34948"])
CORR.set_bad(SURFACE)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})


def dense(smf, positions, length):
    out = np.full((smf.shape[0], length), -1, dtype=np.int8)
    out[:, positions] = smf
    return out


def binned(calls, bin_size):
    """[reads, bins] accessible fraction of called sites (calls: 0 accessible, 1 protected)."""
    n = calls.shape[1] // bin_size
    c = calls[:, : n * bin_size].reshape(len(calls), n, bin_size)
    called = (c != -1).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(called > 0, (c == 0).sum(-1) / called, np.nan)


def corr_matrix(x, min_n=10):
    m = ~np.isnan(x)
    v, mf = np.where(m, x, 0.0), m.astype(float)
    n = mf.T @ mf
    s1, s2, sxy = v.T @ mf, (v * v).T @ mf, v.T @ v
    with np.errstate(invalid="ignore", divide="ignore"):
        var = s2 - s1 * s1 / n
        r = (sxy - s1 * s1.T / n) / np.sqrt(var * var.T)
    r[n < min_n] = np.nan
    return r


def adjust_for_read_level(x):
    """Regress each bin on the read's mean accessibility over the window (per bin, over
    the reads that call it) and return the residuals: correlating them gives the partial
    correlation given the read's overall level, which removes molecule-wide labelling
    efficiency. Unlike subtracting the read mean, the slope is fitted per bin, so open
    bins can respond more strongly than closed ones."""
    m = np.nanmean(x, axis=1, keepdims=True)  # read level over the window
    ok = ~np.isnan(x) & ~np.isnan(m)
    n = ok.sum(0).clip(1)
    xv = np.where(ok, x, 0.0)
    mv = np.where(ok, m, 0.0)
    xbar, mbar = xv.sum(0) / n, mv.sum(0) / n
    cov = (np.where(ok, (x - xbar) * (m - mbar), 0.0)).sum(0) / n
    var = (np.where(ok, (m - mbar) ** 2, 0.0)).sum(0) / n
    slope = np.where(var > 0, cov / np.where(var > 0, var, 1), 0.0)
    return np.where(ok, x - xbar - slope * (m - mbar), np.nan)


def span_label(bp):
    return f"{bp / 1000:g} kb" if bp >= 1000 else f"{bp:g} bp"


regions = list(read_samples_h5(args.samples))
if args.regions:  # plot just these; the summary then covers only them
    regions = [r for r in regions if r.region_id in set(args.regions)]
length = len(regions[0].dna)
corr_bin = 100 if length >= 8000 else max(4, length // 128)
heat_bin = 50 if length >= 8000 else max(2, length // 256)
print(f"{len(regions)} regions, output window {length} bp; co-accessibility bins {corr_bin} bp")


def region_arrays(region):
    obs = dense(region.observed_smf.numpy(), region.observed_positions.numpy(), length)
    gen = dense(region.predicted_smf.numpy(), region.predicted_positions.numpy(), length)
    chosen = rng.choice(len(gen), min(len(obs), len(gen)), replace=False)
    return obs, gen, gen[chosen]


def genomic_start(region):
    chrom, start, _ = parse_region_string(region.region_id)
    return chrom, start - 1 + region.window_start


# ---- summary: co-accessibility by distance, and element pairs ----------------------
n_corr_bins = length // corr_bin
sep = np.abs(np.subtract.outer(np.arange(n_corr_bins), np.arange(n_corr_bins)))
curves = {k: [np.zeros(n_corr_bins), np.zeros(n_corr_bins)] for k in
          ("observed", "generated", "observed_residual", "generated_residual")}
elements = defaultdict(list)
if args.elements:
    with open(args.elements) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            elements[f[0]].append(((int(f[1]) + int(f[2])) // 2, "CTCF-bound" in f[7], f[4] == "active"))
pairs = {"ctcf": [], "active": []}
structure = []
for i, region in enumerate(regions):
    obs, _, gen = region_arrays(region)
    for name, calls in (("observed", obs), ("generated", gen)):
        x = binned(calls, corr_bin)
        for suffix, data in (("", x), ("_residual", adjust_for_read_level(x))):
            r = corr_matrix(data)
            ok = ~np.isnan(r) & (sep > 0)
            np.add.at(curves[name + suffix][0], sep[ok], r[ok])
            np.add.at(curves[name + suffix][1], sep[ok], 1)
    profile = np.nanmean(binned(obs, heat_bin), axis=0)
    structure.append(np.nanstd(profile))
    if args.elements:
        chrom, g0 = genomic_start(region)
        inside = [e for e in elements[chrom] if g0 + 100 <= e[0] < g0 + length - 100]
        for a in range(len(inside)):
            for b in range(a + 1, len(inside)):
                (p, c1, a1), (q, c2, a2) = inside[a], inside[b]
                if abs(q - p) < 2000:
                    continue
                kind = "ctcf" if c1 and c2 else "active" if a1 and a2 and not (c1 or c2) else None
                if kind is None:
                    continue
                row = []
                for calls in (obs, gen):
                    x = binned(calls[:, p - g0 - 100 : p - g0 + 100], 200)[:, 0]
                    y = binned(calls[:, q - g0 - 100 : q - g0 + 100], 200)[:, 0]
                    mm = np.nanmean(binned(calls, corr_bin), axis=1)
                    ok = ~np.isnan(x) & ~np.isnan(y)
                    u, v = x[ok] - mm[ok], y[ok] - mm[ok]
                    row.append(float(np.corrcoef(u, v)[0, 1]) if ok.sum() >= 20 and u.std() and v.std() else np.nan)
                pairs[kind].append((abs(q - p), *row, region.region_id))

summary = {"n_regions": len(regions), "output_length": length, "corr_bin": corr_bin,
           "separation_bp": (np.arange(n_corr_bins) * corr_bin).tolist()}
for k, (s, n) in curves.items():
    with np.errstate(invalid="ignore"):
        summary[k] = (s / n).tolist()
for kind, rows in pairs.items():
    a = np.array([r[1:3] for r in rows], dtype=float)
    ok = ~np.isnan(a).any(axis=1) if len(a) else np.zeros(0, bool)
    summary[f"pairs_{kind}"] = {
        "n": int(ok.sum()),
        "observed_mean": float(a[ok, 0].mean()) if ok.any() else None,
        "generated_mean": float(a[ok, 1].mean()) if ok.any() else None,
        "observed_vs_generated_r": float(np.corrcoef(a[ok, 0], a[ok, 1])[0, 1]) if ok.sum() > 2 else None,
    }
if not args.no_summary:
    with open(f"{args.out}/summary.json", "w") as handle:
        json.dump(summary, handle)

d = np.array(summary["separation_bp"])[1:]
fmt = FuncFormatter(lambda v, _: span_label(v))
n_panels = 2 if any(summary[f"pairs_{k}"]["n"] for k in pairs) else 1
fig, axes = plt.subplots(1, n_panels, figsize=(5.8 * n_panels, 4.2), facecolor=SURFACE, squeeze=False)
ax = axes[0, 0]
for name, color in (("observed", OBSERVED), ("generated", GENERATED)):
    ax.plot(d, summary[name][1:], color=color, linewidth=2, label=f"{name.capitalize()}")
    ax.plot(d, summary[name + "_residual"][1:], color=color, linewidth=1.3,
            linestyle=(0, (4, 2)), label=f"{name.capitalize()}, adjusted for read labelling level")
ax.axhline(0, color=AXIS, linewidth=0.8)
ax.set_xscale("log")
ax.xaxis.set_major_formatter(fmt)
ax.set_xlabel(f"Distance between {corr_bin} bp bins")
ax.set_ylabel("Correlation of accessibility across reads")
ax.set_title("Co-accessibility by distance", loc="left", fontsize=10)
ax.legend(frameon=False, fontsize=8)
if n_panels == 2:
    ax2 = axes[0, 1]
    lim = [-0.6, 0.9]
    ax2.plot(lim, lim, color=AXIS, linewidth=0.8)
    for kind, color, label in (("ctcf", OBSERVED, "CTCF-bound pairs"), ("active", MUTED, "Other active pairs")):
        a = np.array([r[1:3] for r in pairs[kind]], dtype=float)
        if len(a):
            ax2.scatter(a[:, 0], a[:, 1], s=10, color=color, alpha=0.6, linewidths=0, label=label)
    ax2.set_xlim(lim)
    ax2.set_ylim(lim)
    ax2.set_xlabel("Observed r (each read's mean removed)")
    ax2.set_ylabel("Generated r")
    s = summary["pairs_ctcf"]
    ax2.set_title(f"Element pairs >= 2 kb apart (CTCF: {s['n']} pairs)", loc="left", fontsize=10)
    ax2.legend(frameon=False, fontsize=8, loc="upper left")
for a in axes.flat:
    a.set_facecolor(SURFACE)
    a.grid(axis="y", color=GRID, linewidth=0.6)
    a.set_axisbelow(True)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)
fig.text(0.01, 0.01, f"{len(regions)} regions, {length} bp output windows; generated reads subsampled to "
         "each region's observed read count.", color=MUTED, fontsize=7.5)
fig.tight_layout(rect=(0, 0.04, 1, 1))
if not args.no_summary:
    fig.savefig(f"{args.out}/summary.png", dpi=150, facecolor=SURFACE)
plt.close(fig)

# ---- example regions --------------------------------------------------------------
if args.regions:
    chosen = [i for i, r in enumerate(regions) if r.region_id in set(args.regions)]
else:
    order = np.argsort(structure)[::-1]
    chosen = list(order[: args.n_examples // 2])
    rest = [i for i in range(len(regions)) if i not in chosen]
    chosen += list(rng.choice(rest, args.n_examples - len(chosen), replace=False))
    ctcf_regions = sorted({r[3] for r in pairs["ctcf"]})
    if ctcf_regions:
        ids = [r.region_id for r in regions]
        best = max(ctcf_regions, key=lambda rid: structure[ids.index(rid)])
        if ids.index(best) not in chosen:
            chosen.append(ids.index(best))

ZOOM, ZOOM_BIN, CORR_LIMIT = 3000, 20, 0.6


def zoom_centre(obs, chrom, g0):
    """Most accessible active element at least ZOOM / 2 from the window's ends, or the
    observed profile's peak; in window coordinates."""
    half = ZOOM // 2
    candidates = [p - g0 for p, ctcf, active in elements.get(chrom, [])
                  if active and g0 + half <= p < g0 + length - half]
    # 200 bp means of the observed profile, every 50 bp
    profile = np.array([np.nanmean(binned(obs[:, k : k + 200], 200)) for k in range(0, length - 200, 50)])
    centres = np.arange(len(profile)) * 50 + 100
    allowed = (centres >= half) & (centres < length - half)
    if candidates:
        # candidates.bed centres are thinned and sit off the peak, so re-centre on the
        # observed peak within 1 kb of the most accessible active cCRE
        score = [np.nanmean(binned(obs[:, c - 100 : c + 100], 200)) for c in candidates]
        best = candidates[int(np.nanargmax(score))]
        allowed &= np.abs(centres - best) <= 1000
        reason = "accessibility peak near the most accessible active cCRE"
    else:
        reason = "observed accessibility peak"
    return int(centres[np.nanargmax(np.where(allowed, profile, np.nan))]), reason


for i in chosen:
    region = regions[i]
    obs, gen_all, gen = region_arrays(region)
    chrom, g0 = genomic_start(region)
    zc, zoom_reason = zoom_centre(obs, chrom, g0)
    z0, z1 = zc - ZOOM // 2, zc + ZOOM // 2
    fig = plt.figure(figsize=(13.5, 15), facecolor=SURFACE)
    grid = fig.add_gridspec(4, 6, height_ratios=[0.7, 2.0, 2.6, 2.6], hspace=0.42, wspace=0.5)
    prof = fig.add_subplot(grid[0, :])
    x = (np.arange(length // heat_bin) + 0.5) * heat_bin
    for name, calls, color in (("Observed", obs, OBSERVED), (f"Generated ({len(gen_all)} reads)", gen_all, GENERATED)):
        prof.plot(x, 1 - np.nanmean(binned(calls, heat_bin), axis=0), color=color, linewidth=1.5, label=name)
    prof.axvspan(z0, z1, color=GRID, alpha=0.6, linewidth=0, zorder=0)
    prof.set_xlim(0, length)
    prof.set_ylabel("Protected")
    prof.legend(frameon=False, fontsize=8, ncol=2, loc="lower left")
    prof.set_title(f"{chrom}:{g0 + 1:,}-{g0 + length:,} ({span_label(length)}), {len(obs)} observed reads; "
                   f"shaded: the zoomed {span_label(ZOOM)} around the {zoom_reason}", loc="left", fontsize=10)
    if args.elements:
        for p, ctcf, active in elements[chrom]:
            if g0 <= p < g0 + length and (ctcf or active):
                prof.axvline(p - g0, color=OBSERVED if ctcf else AXIS, linewidth=0.8,
                             linestyle="-" if ctcf else ":", zorder=0)
    heat_axes = []
    for col, (name, calls) in enumerate((("Observed reads", obs), (f"Generated reads (as many as observed)", gen))):
        ax = fig.add_subplot(grid[1, 3 * col : 3 * col + 3])
        heat_axes.append(ax)
        b = 1 - binned(calls, heat_bin)  # protected fraction
        b = b[np.argsort(np.nanmean(b, axis=1))]
        heat = ax.imshow(b, aspect="auto", cmap=ACCESS, vmin=0, vmax=1, interpolation="nearest",
                         extent=(0, length, len(b), 0))
        ax.set_title(name, loc="left", fontsize=9.5)
        ax.set_ylabel("Reads (sorted by mean protection)" if col == 0 else "")
        ax.set_xlabel(f"Position ({heat_bin} bp bins)")
    heat_bar = fig.colorbar(heat, ax=heat_axes, shrink=0.8, pad=0.02)
    heat_bar.set_label("Protected fraction (blank: no call)", color=INK2)
    corr_axes = []
    maps = (("Observed", obs), ("Generated, as many reads", gen), (f"Generated, all {len(gen_all)} reads", gen_all))
    for row, (lo, hi, bin_size, label) in enumerate(((0, length, corr_bin, "whole window"),
                                                       (z0, z1, ZOOM_BIN, "zoom"))):
        for col, (name, calls) in enumerate(maps):
            ax = fig.add_subplot(grid[2 + row, 2 * col : 2 * col + 2])
            corr_axes.append(ax)
            im = ax.imshow(corr_matrix(adjust_for_read_level(binned(calls[:, lo:hi], bin_size))), cmap=CORR,
                           vmin=-CORR_LIMIT, vmax=CORR_LIMIT, interpolation="nearest",
                           extent=(lo, hi, hi, lo))
            ax.set_title(f"{name} ({label})", loc="left", fontsize=9)
            ax.set_xlabel(f"Position ({bin_size} bp bins)")
            if row == 0:
                ax.add_patch(plt.Rectangle((z0, z0), ZOOM, ZOOM, fill=False, edgecolor=INK2, linewidth=1))
    cbar = fig.colorbar(im, ax=corr_axes, shrink=0.6, pad=0.02, extend="both")
    cbar.set_label(f"Correlation given read labelling level (clipped at +/-{CORR_LIMIT})", color=INK2)
    for a in [prof, *heat_axes, *corr_axes]:
        a.set_facecolor(SURFACE)
        a.xaxis.set_major_formatter(fmt)
    for a in corr_axes:
        a.yaxis.set_major_formatter(fmt)
    path = f"{args.out}/example_{region.region_id.replace(':', '_')}.png"
    fig.savefig(path, dpi=120, facecolor=SURFACE)
    plt.close(fig)
    print("wrote", path)

    # The whole window, large: raw and labelling-adjusted, observed vs all generated reads.
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 12), facecolor=SURFACE)
    for r, (row_label, transform) in enumerate((("raw", lambda b: b),
                                                ("adjusted for read labelling level", adjust_for_read_level))):
        for c, (name, calls) in enumerate(((f"Observed, {len(obs)} reads", obs),
                                           (f"Generated, {len(gen_all)} reads", gen_all))):
            ax = axes[r, c]
            im = ax.imshow(corr_matrix(transform(binned(calls, corr_bin))), cmap=CORR,
                           vmin=-CORR_LIMIT, vmax=CORR_LIMIT, interpolation="nearest",
                           extent=(0, length, length, 0))
            ax.set_title(f"{name}, {row_label}", loc="left", fontsize=10)
            ax.set_xlabel(f"Position ({corr_bin} bp bins)")
            ax.xaxis.set_major_formatter(fmt)
            ax.yaxis.set_major_formatter(fmt)
            ax.set_facecolor(SURFACE)
    cbar = fig.colorbar(im, ax=axes, shrink=0.6, pad=0.02, extend="both")
    cbar.set_label(f"Correlation across reads (clipped at +/-{CORR_LIMIT})", color=INK2)
    fig.suptitle(f"{chrom}:{g0 + 1:,}-{g0 + length:,} ({span_label(length)}), whole-window "
                 "co-accessibility", x=0.06, ha="left", fontsize=11)
    path = f"{args.out}/full_{region.region_id.replace(':', '_')}.png"
    fig.savefig(path, dpi=110, facecolor=SURFACE)
    plt.close(fig)
    print("wrote", path)
print(json.dumps({k: v for k, v in summary.items() if k.startswith("pairs")}, indent=1))
