"""Schematic of the current model's transition topology (FixedBlockTransition), for the report.

Usage: python plot_transition_topology.py OUT.png
Hand-drawn schematic, not computed from a checkpoint; block ranges follow the hist190_tf60 run
(TF 15-60, other_long 40-100, histone 100-190 with a 40-state exit tail; 10 bp TF spacers, no chaining).
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
C = {"open": "#8a8a86", "TF": "#2a78d6", "other_long": "#eb6834", "histone": "#4a3aa7"}
fig, ax = plt.subplots(figsize=(12, 5.6), facecolor=SURF)
ax.set_facecolor(SURF); ax.set_xlim(0, 12.4); ax.set_ylim(-0.6, 5.0); ax.axis("off")


def box(x0, x1, y, text, group, sub=None, h=0.62):
    ax.add_patch(FancyBboxPatch((x0, y - h / 2), x1 - x0, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                facecolor=C[group], alpha=0.16, edgecolor=C[group], linewidth=1.6))
    ax.text((x0 + x1) / 2, y + (0.08 if sub else 0), text, ha="center", va="center", fontsize=10, color=INK)
    if sub:
        ax.text((x0 + x1) / 2, y - 0.17, sub, ha="center", va="center", fontsize=8, color=INK2)


def arrow(p, q, rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=12, color=INK2, lw=1.2,
                                 connectionstyle=f"arc3,rad={rad}"))


rows = {"TF": 3.6, "other_long": 2.1, "histone": 0.6}
box(0.1, 1.2, 2.1, "open", "open", "hub", h=3.6)
box(11.2, 12.3, 2.1, "open", "open", "same state", h=3.6)
# TF row: entry spacer, TF chain, exit spacer
y = rows["TF"]
box(1.8, 3.6, y, "entry spacer", "open", "10 open states")
box(4.0, 8.2, y, "TF chain", "TF", "15-60 bp; exit from ≥ 15, forced at 60")
box(8.6, 10.6, y, "exit spacer", "open", "10 open states → open only")
arrow((1.2, y), (1.8, y)); arrow((3.6, y), (4.0, y)); arrow((8.2, y), (8.6, y)); arrow((10.6, y), (11.2, y))
# other_long row
y = rows["other_long"]
box(1.8, 10.6, y, "other_long chain", "other_long", "40-100 bp; exit from ≥ 40, forced at 100")
arrow((1.2, y), (1.8, y)); arrow((10.6, y), (11.2, y))
# histone row: variable part then exit tail
y = rows["histone"]
box(1.8, 3.6, y, "entry taper", "histone", "first 40 states")
box(3.6, 7.6, y, "variable core", "histone", "length set by exits into the tail")
box(7.6, 10.6, y, "exit tail (taper)", "histone", "last 40 states, every nucleosome")
arrow((1.2, y), (1.8, y)); arrow((10.6, y), (11.2, y))
for x in (5.0, 6.0, 7.0):
    arrow((x, y + 0.31), (7.9, y + 0.31), rad=-0.5)
ax.text(6.1, y - 0.62, "histone 100-190 bp; zero-bias P(d) ∝ exp(−|d − 147| / 15)", ha="center", fontsize=8.5, color=INK2)
# prior
ax.text(0.1, 4.75, "Sequence prior (frozen AlphaGenome encoder + GELU head) adds one potential per state group and position:",
        fontsize=9.5, color=INK, ha="left")
ax.text(0.1, 4.42, "every state of a block shares its block's potential (tied); spacer states use the open potential. "
        "Transition probabilities are frozen; the prior only reweights branch points.", fontsize=8.5, color=INK2, ha="left")
fig.tight_layout()
fig.savefig(sys.argv[1], dpi=120, facecolor=SURF)
