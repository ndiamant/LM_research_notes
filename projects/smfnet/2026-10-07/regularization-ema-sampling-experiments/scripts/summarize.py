"""Summarize the overnight regularization runs from their CSV logs.

Writes summary.tsv (one row per run), loss_curves.png (train and validation
loss per epoch, one panel per run) and param_norms.png (parameter norm growth
by part of the network, relative to the first epoch).

    python analyses/overnight_regularization_20261006/summarize.py
"""

from pathlib import Path
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(os.path.expandvars("$SCRATCH/smf_models/overnight_20261006"))
VAL_COLOR, TRAIN_COLOR, INK, GRID = "#2a78d6", "#8a8984", "#52514e", "#e1e0d9"
PART_COLORS = {"DNA encoder": "#2a78d6", "read blocks": "#eb6834",
               "DNA conditioning": "#1baf7a", "input layer": "#eda100", "output layer": "#e87ba4"}


def part(group: str) -> str:
    if group.startswith("dna_encoder"):
        return "DNA encoder"
    if group.endswith(".cond"):
        return "DNA conditioning"
    if group.startswith("body"):
        return "read blocks"
    return {"input_conv": "input layer", "output_conv": "output layer"}.get(group, group)


def load(run: Path) -> pd.DataFrame | None:
    path = run / "csv" / "metrics.csv"
    return pd.read_csv(path) if path.exists() else None


def train_regions(run: Path) -> int | None:
    """Training regions per H5, from the run's saved Hydra config."""
    import yaml

    data = yaml.safe_load((run / ".hydra" / "config.yaml").read_text())["data"]["cfg"]
    return data.get("max_train_regions_per_h5") or data.get("max_regions_per_h5_per_split")


def complete_epochs(metrics: pd.DataFrame) -> pd.DataFrame:
    return metrics.groupby("epoch").agg({"val/loss": "last", "train/loss_epoch": "last"}).dropna()


def summarize(name: str, metrics: pd.DataFrame) -> dict:
    epochs = complete_epochs(metrics)
    best = epochs["val/loss"].idxmin()
    after = epochs.loc[best:]
    row = dict(run=name, train_regions=train_regions(ROOT / name),
               epochs=int(epochs.index.max()) + 1, best_epoch=int(best),
               best_val_loss=epochs.loc[best, "val/loss"],
               train_loss_at_best=epochs.loc[best, "train/loss_epoch"],
               final_val_loss=epochs["val/loss"].iloc[-1],
               final_train_loss=epochs["train/loss_epoch"].iloc[-1])
    row["gap_at_best"] = row["best_val_loss"] - row["train_loss_at_best"]
    row["val_rise_after_best"] = row["final_val_loss"] - row["best_val_loss"]
    # Single-epoch validation losses are noisy (about +-0.003), so also report
    # steadier summaries: the lowest 5-epoch rolling mean and the mean of the
    # last 20 epochs, with the train/validation gap over the same epochs.
    rolling = epochs["val/loss"].rolling(5, center=True).mean()
    row["best_rolling5_val"] = rolling.min()
    row["best_rolling5_epoch"] = int(rolling.idxmin()) if rolling.notna().any() else -1
    last = epochs.tail(20)
    row["val_last20"] = last["val/loss"].mean()
    row["gap_last20"] = (last["val/loss"] - last["train/loss_epoch"]).mean()
    # The same epochs for every run, because early stopping or faster epochs
    # change which epochs the last 20 are.
    window = epochs.loc[70:89]
    row["val_epochs70_89"] = window["val/loss"].mean() if len(window) >= 15 else np.nan
    row["gap_epochs70_89"] = ((window["val/loss"] - window["train/loss_epoch"]).mean()
                              if len(window) >= 15 else np.nan)
    # Rise per 10 epochs after the minimum, from a line fit (0 if too few points).
    row["val_rise_per_10_epochs"] = (10 * np.polyfit(after.index, after["val/loss"], 1)[0]
                                     if len(after) >= 5 else np.nan)
    return row


def style(axis):
    axis.grid(color=GRID, linewidth=0.6)
    axis.set_axisbelow(True)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)


def main():
    runs = {p.name: load(p) for p in sorted(ROOT.iterdir())
            if p.is_dir() and p.name not in {"logs", "wandb", "smoke"}}
    # Runs that are just starting have no complete epochs yet.
    runs = {k: v for k, v in runs.items() if v is not None and len(complete_epochs(v)) >= 5}
    summary = pd.DataFrame([summarize(k, v) for k, v in runs.items()]).sort_values(
        ["train_regions", "val_epochs70_89"])
    summary.to_csv(ROOT / "summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(summary.round(4).to_string(index=False))

    names = list(summary.run)
    cols = min(4, len(names))
    rows = int(np.ceil(len(names) / cols))
    figure, axes = plt.subplots(rows, cols, figsize=(3.6 * cols, 2.9 * rows), sharey=True,
                                squeeze=False, layout="constrained")
    for axis, name in zip(axes.flat, names):
        epochs = runs[name].groupby("epoch").agg({"val/loss": "last", "train/loss_epoch": "last"})
        axis.plot(epochs.index, epochs["train/loss_epoch"], color=TRAIN_COLOR, lw=1.5, ls="--")
        axis.plot(epochs.index, epochs["val/loss"], color=VAL_COLOR, lw=2)
        row = summary.set_index("run").loc[name]
        axis.scatter([row.best_epoch], [row.best_val_loss], s=40, color=VAL_COLOR,
                     edgecolors="white", linewidths=1.5, zorder=3)
        axis.set_title(f"{name} ({int(row.train_regions):,} train regions)\n"
                       f"last-20 mean {row.val_last20:.4f}, gap {row.gap_last20:.3f}", fontsize=9)
        axis.set_xlabel("Epoch", fontsize=8, color=INK)
        style(axis)
    for axis in axes.flat[len(names):]:
        axis.set_visible(False)
    low = min(summary.train_loss_at_best.min(), summary.best_val_loss.min())
    axes[0, 0].set_ylim(low - 0.01, summary.best_val_loss.max() + 0.06)
    axes[0, 0].set_ylabel("Loss", fontsize=9)
    figure.legend(handles=[plt.Line2D([], [], color=VAL_COLOR, lw=2, label="Validation"),
                           plt.Line2D([], [], color=TRAIN_COLOR, lw=1.5, ls="--", label="Training")],
                  loc="outside lower center", ncols=2, frameon=False, fontsize=9)
    figure.suptitle("Loss by epoch, ES proxies (same 15,000 validation regions)",
                    fontsize=11, x=0.01, ha="left")
    figure.savefig(ROOT / "loss_curves.png", dpi=110)
    plt.close(figure)

    figure, axes = plt.subplots(rows, cols, figsize=(3.6 * cols, 2.9 * rows), sharey=True,
                                squeeze=False, layout="constrained")
    for axis, name in zip(axes.flat, names):
        metrics = runs[name]
        columns = [c for c in metrics.columns if c.startswith("param_norm/")]
        norms = metrics.groupby("epoch")[columns].last().dropna(how="all")
        if norms.empty:
            axis.set_visible(False)
            continue
        squares = (norms**2).T.groupby(lambda c: part(re.sub("^param_norm/", "", c))).sum().T**0.5
        relative = squares / squares.iloc[0]
        for label in PART_COLORS:
            if label in relative:
                axis.plot(relative.index, relative[label], color=PART_COLORS[label], lw=2)
        best = summary.set_index("run").loc[name, "best_epoch"]
        axis.axvline(best, color=INK, lw=0.8, ls=":")
        axis.set_title(name, fontsize=9)
        axis.set_xlabel("Epoch", fontsize=8, color=INK)
        style(axis)
    for axis in axes.flat[len(names):]:
        axis.set_visible(False)
    axes[0, 0].set_ylabel("Parameter norm / first epoch", fontsize=9)
    figure.legend(handles=[plt.Line2D([], [], color=c, lw=2, label=k) for k, c in PART_COLORS.items()],
                  loc="outside lower center", ncols=5, frameon=False, fontsize=8)
    figure.suptitle("Parameter norm growth by part of the network (dotted = best epoch)",
                    fontsize=11, x=0.01, ha="left")
    figure.savefig(ROOT / "param_norms.png", dpi=110)
    plt.close(figure)
    print(f"Wrote {ROOT / 'summary.tsv'}, loss_curves.png and param_norms.png")


if __name__ == "__main__":
    main()
