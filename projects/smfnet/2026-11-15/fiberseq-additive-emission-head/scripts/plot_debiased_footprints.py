"""TF footprints of the additive-emission model with and without labelling bias.

From marginalize.py runs with sampling.record_debiased: per TF, read- and background-mean
P(m6A) around the inserted motif, as sampled ("with bias": the +/-3 bp context table and
the read's b_m) and the network's own ("bias removed": a mean-context site on a median
molecule), for the motif (solid) and its matched control (dashed). Reverse-orientation
insertions are flipped onto the forward orientation and averaged with them. Outside the
motif the curves are an 11 bp running mean (each position averages only the ~half of the
random backgrounds with an A/T there); inside it, single bases (C/G: no calls, gaps). The
last row zooms the first run to +/-40 bp, unsmoothed. Usage: plot_debiased_footprints.py OUT.png RUN_DIR[:LABEL] ...
(one row of panels per run; the first run's motifs set the panel order).
"""
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HALF = 150
SURFACE, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
ORANGE, BLUE = "#eb6834", "#2a78d6"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED})


def profiles(path):
    """{tf: {(kind, condition): mean profile over +/-HALF around the motif centre}}."""
    out = {}
    with h5py.File(path) as h:
        names = [n.decode() if isinstance(n, bytes) else n for n in h["motif_names"][:]]
        orient = [o.decode() if isinstance(o, bytes) else o for o in h["motif_orientations"][:]]
        seqs = [s.decode() if isinstance(s, bytes) else s for s in h["motif_sequences"][:]]
        starts = h["insertion_starts"][:]
        offset = int(h.attrs["output_offset"])
        for m, (name, o) in enumerate(zip(names, orient, strict=True)):
            tf = name.split("(")[0]
            centre = int(starts[m]) - offset + len(seqs[m]) // 2
            for kind in ("m6a", "m6a_debiased"):
                for cond in ("motif", "control"):
                    p = np.nanmean(h[f"profiles/{kind}/{cond}"][m], axis=0)  # over backgrounds
                    w = p[centre - HALF : centre + HALF + 1]
                    if o == "-":
                        w = w[::-1]
                    out.setdefault(tf, {}).setdefault((kind, cond), []).append(w)
            out[tf]["length"] = len(seqs[m])
    return {tf: {k: (np.nanmean(v, 0) if k != "length" else v) for k, v in d.items()} for tf, d in out.items()}


out_png, runs = sys.argv[1], sys.argv[2:]
loaded = []
for r in runs:
    d, _, label = r.partition(":")
    loaded.append((label or d, profiles(f"{d}/marginalizations.h5")))
tfs = list(loaded[0][1])
x = np.arange(-HALF, HALF + 1)
SMOOTH = 11


def smoothed(y, inside):
    """Running nanmean outside the motif; the motif's own positions as they are."""
    k = SMOOTH // 2
    out = np.array([np.nanmean(y[max(i - k, 0) : i + k + 1]) if not inside[i] else np.nan
                    for i in range(len(y))])
    return out


rows = loaded + [(f"{loaded[0][0]}, zoom", loaded[0][1])]
fig, axes = plt.subplots(len(rows), len(tfs), figsize=(2.1 * len(tfs), 2.2 * len(rows) + 0.4),
                         facecolor=SURFACE, squeeze=False)
for row, (label, prof) in enumerate(rows):
    zoom = row == len(rows) - 1
    for col, tf in enumerate(tfs):
        ax, p = axes[row, col], prof[tf]
        half_len = p["length"] / 2
        inside = np.abs(x) <= half_len
        ax.axvspan(-half_len, half_len, color=GRID, alpha=0.6, lw=0)
        for kind, color, name in (("m6a", ORANGE, "with bias"), ("m6a_debiased", BLUE, "bias removed")):
            for cond, ls, lw, alpha in (("motif", "-", 1.4, 1.0), ("control", "--", 0.9, 0.8)):
                y = p[(kind, cond)]
                lab = f"{cond}, {name}"
                if zoom:
                    ax.plot(x, y, color=color, ls=ls, lw=lw * 0.8, alpha=alpha, marker=".", ms=3, label=lab)
                else:
                    ax.plot(x, smoothed(y, inside), color=color, ls=ls, lw=lw, alpha=alpha, label=lab)
                    if cond == "motif":
                        ax.plot(x[inside], y[inside], color=color, ls="none", marker="o", ms=2.2)
        if zoom:
            ax.set_xlim(-40, 40)
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if row == 0:
            ax.set_title(tf, fontsize=9)
        if col == 0:
            ax.set_ylabel(f"{label}\nP(m6A)")
        if row >= len(rows) - 2:
            ax.set_xlabel("bp from motif centre")
axes[0, 0].legend(frameon=False, fontsize=6, loc="lower left")
fig.tight_layout()
fig.savefig(out_png, dpi=130, facecolor=SURFACE)
print("wrote", out_png)
