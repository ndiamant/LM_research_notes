"""Does the lack of a per-read m6A-rate filter hurt the Fiber-seq HMMs or the ChIP benchmark?

Usage: analyze_read_rate.py TABLE_DIR   (the output of build_read_tables.py for chr8, chr20)

Rate classes use the whole-molecule m6A rate: low < 2%, in-range 2-8% (Bubb et al. 2026's
recommended range), high >= 8% (split at 12% where noted). Writes summary.txt and figures to
TABLE_DIR:

A. who the benchmark reads are: class shares, molecule vs window rate, class mix at bound vs
   unbound sites;
B. model state occupancy by molecule rate (are states used to absorb labeling rate?);
C. m6A autocorrelation along reads by rate (nucleosome periodicity; Bubb et al.'s concern is
   that heavy labeling methylates nucleosomal DNA);
D. m6A profiles at IDR-bound vs unbound CTCF / NRF1 motifs by rate (footprint depth);
E. per-class information about binding at matched read counts: sites with >= k reads of every
   compared class, k sampled per class, AUROC bound vs unbound of read-only and model scores;
F. simulated filtering: site scores from in-range reads only vs all reads vs random subsets of
   the same per-site size, scored with the benchmark's IDR metrics;
G. mean P(TF) on footprint reads at bound sites by rate class (footprints called TF?).
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/users/diamant/repos/smf_hmm")
from evaluation.tf_chip.analyze_covariation import auroc, partial_spearman, stratified_auroc  # noqa: E402

D = Path(sys.argv[1])
CHROMS = ("chr8", "chr20")
MODELS = ["phase_type", "ag_kl10", "ag_ol001", "ag_ol001_s05", "ag_kl10_s0"]
FOOTPRINT_TFS = ["ctcf", "nrf1", "usf2", "crem"]
TFS = FOOTPRINT_TFS + ["rfx5", "ebf1", "pax5", "spi1", "yy1"]
N_DRAWS = 20
rng = np.random.default_rng(0)
out = open(D / "summary.txt", "w")


def say(*args):
    print(*args)
    print(*args, file=out, flush=True)


def rate_class(rate, split_high=False):
    c = np.where(rate < 0.02, "low", np.where(rate < 0.08, "in_range", "high"))
    if split_high:
        c = np.where(rate >= 0.12, "high12", c)
    return np.where(np.isnan(rate), "unknown", c)


reads = pd.concat([pd.read_csv(D / f"reads_{c}.tsv.gz", sep="\t") for c in CHROMS], ignore_index=True)
rs = pd.concat([pd.read_csv(D / f"read_sites_{c}.tsv.gz", sep="\t") for c in CHROMS], ignore_index=True)
sites = pd.concat([pd.read_csv(f"/scratch/users/diamant/gm12878_validation/covariation_windows_v2/{c}-512-valid90-reads10-flank128/analysis_idr/unique_motif_summaries.tsv.gz",
                               sep="\t", usecols=["site_key", "tf", "idr_label", "chip_pooled_sum", "accessibility"])
                   for c in CHROMS], ignore_index=True).set_index("site_key")
reads["cls"] = rate_class(reads.mol_rate.to_numpy())
rs["rate_cls"] = rate_class(rs.mol_rate.to_numpy())
rs["rate_cls4"] = rate_class(rs.mol_rate.to_numpy(), split_high=True)
pd.set_option("display.width", 250)

# A. Who the reads are.
say("== A. Benchmark reads (val chr8 + chr20)")
say(f"reads {len(reads):,}; without molecule rate {reads.mol_rate.isna().mean():.2%}")
say("molecule-rate percentiles 1/5/25/50/75/95/99:",
    np.round(np.nanpercentile(reads.mol_rate, [1, 5, 25, 50, 75, 95, 99]), 3).tolist())
say("class shares:", reads.cls.value_counts(normalize=True).round(3).to_dict())
win_rate = reads.win_m6a / reads.win_called.replace(0, np.nan)
say(f"Spearman(molecule rate, 512 bp window rate) = {reads.mol_rate.rank().corr(win_rate.rank()):.2f}")
say("Mean window rate by molecule class:", win_rate.groupby(reads.cls).mean().round(3).to_dict())
mix = rs[rs.idr_label != "ambiguous"].groupby(["tf", "idr_label"]).rate_cls.value_counts(normalize=True).unstack().round(3)
say("Class mix of reads at bound vs unbound sites:\n" + mix.to_string())

# B. State occupancy by molecule rate.
bins = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 1.01])
reads["bin"] = pd.cut(reads.mol_rate, bins, right=False)
say("\n== B. Mean posterior occupancy over the 512 bp window, by molecule rate")
fig, axes = plt.subplots(1, len(MODELS), figsize=(4 * len(MODELS), 3.4), sharey=True)
centers = (bins[:-1] + np.minimum(bins[1:], 0.3)) / 2
for ax, m in zip(axes, MODELS):
    cols = [c for c in reads.columns if c.startswith(m + "_") and not c.startswith(m + "_s")
            and c[len(m) + 1:] in ("open", "TF", "other_long", "histone")]
    t = reads.groupby("bin", observed=False)[cols].mean()
    t.columns = [c[len(m) + 1:] for c in cols]
    t["n_reads"] = reads.groupby("bin", observed=False).size()
    say(f"-- {m}\n" + t.round(3).to_string())
    for g in t.columns[:-1]:
        ax.plot(centers, t[g], marker="o", label=g)
    ax.axvspan(0.02, 0.08, color="0.9", zorder=0)
    ax.set(title=m, xlabel="molecule m6A rate")
axes[0].set_ylabel("mean posterior occupancy")
axes[0].legend(fontsize=7, frameon=False)
fig.suptitle("State occupancy by molecule m6A rate (gray: Bubb et al. 2-8%)", fontsize=10)
fig.tight_layout()
fig.savefig(D / "occupancy_by_rate.png", dpi=110)
plt.close(fig)

# C. Autocorrelation.
acc = {c: np.load(D / f"accumulators_{c}.npz") for c in CHROMS}
num = sum(a["acf_num"] for a in acc.values())
den = sum(a["acf_den"] for a in acc.values())
acf = num / np.where(den > 0, den, np.nan)
lags = np.arange(acf.shape[1])
say("\n== C. m6A autocorrelation along reads (standardized by each read's window rate)")
fig, ax = plt.subplots(figsize=(6.5, 3.6))
cmap = plt.get_cmap("viridis")
rows = []
for i in range(acf.shape[0]):
    label = f"{bins[i]:.2f}-{min(bins[i + 1], 1):.2f}"
    ax.plot(lags[5:], acf[i, 5:], color=cmap(i / (acf.shape[0] - 1)), label=label)
    peak = 140 + int(np.nanargmax(acf[i, 140:241]))
    trough = 60 + int(np.nanargmin(acf[i, 60:131]))
    rows.append(dict(bin=label, peak_lag=peak, peak=acf[i, peak], trough_lag=trough, trough=acf[i, trough],
                     amplitude=acf[i, peak] - acf[i, trough], acf_lag10=acf[i, 10]))
say(pd.DataFrame(rows).round(3).to_string(index=False))
ax.axhline(0, color="0.5", linewidth=0.6)
ax.set(xlabel="lag (bp)", ylabel="m6A autocorrelation", title="m6A autocorrelation by molecule m6A rate")
ax.legend(fontsize=6, ncols=2, frameon=False)
fig.tight_layout()
fig.savefig(D / "autocorrelation_by_rate.png", dpi=110)
plt.close(fig)

# D. Motif profiles.
say("\n== D. m6A profiles at CTCF / NRF1 motifs (strand oriented), by molecule rate class")
groups = {"low": [0], "in_range": [1, 2, 3], "high 8-12%": [4, 5], "high >=12%": [6, 7, 8]}
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
rows = []
for ax, tf in zip(axes, ("ctcf", "nrf1")):
    for label in ("bound", "unbound"):
        k = f"prof_m6a_{tf}_{label}"
        if k not in acc["chr8"]:
            continue
        m6a = sum(a[k] for a in acc.values())
        called = sum(a[f"prof_called_{tf}_{label}"] for a in acc.values())
        x = np.arange(m6a.shape[1]) - m6a.shape[1] // 2
        for gi, (g, idx) in enumerate(groups.items()):
            r = m6a[idx].sum(0) / np.maximum(called[idx].sum(0), 1)
            flank = (np.abs(x) >= 20) & (np.abs(x) <= 60)
            core = np.abs(x) <= 8
            rel = r / r[flank].mean()
            rows.append(dict(tf=tf, sites=label, cls=g, calls=int(called[idx].sum()), core_rate=r[core].mean(),
                             flank_rate=r[flank].mean(), footprint_depth=1 - r[core].mean() / r[flank].mean()))
            if label == "bound":
                ax.plot(x, pd.Series(rel).rolling(5, center=True).mean(), color=plt.get_cmap("plasma")(gi / 4), label=g)
    ax.set(title=f"{tf.upper()} IDR-bound: m6A rate / 20-60 bp flank rate", xlabel="bp from motif center")
    ax.legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(D / "motif_profiles_by_rate.png", dpi=110)
plt.close(fig)
say(pd.DataFrame(rows).round(3).to_string(index=False))


# E. Matched read counts per class.
def site_scores(sub: pd.DataFrame) -> pd.DataFrame:
    g = sub.groupby("site_key")
    fp = (sub.cls == "footprint").groupby(sub.site_key).sum()
    op = (sub.cls == "open").groupby(sub.site_key).sum()
    s = pd.DataFrame({m: g[m].mean() for m in MODELS})
    s["ref_occupancy"] = fp / (fp + op).replace(0, np.nan)
    s["ref_protected_fraction"] = g.protected_in_motif.mean()
    return s


def matched_auroc(frame: pd.DataFrame, column: str, classes: list[str], k: int, tfs: list[str]) -> pd.DataFrame:
    rows = []
    for tf in tfs:
        t = frame[(frame.tf == tf) & (frame.idr_label != "ambiguous")]
        counts = t.groupby(["site_key", column]).size().unstack(fill_value=0)
        eligible = counts.index[(counts.reindex(columns=classes, fill_value=0) >= k).all(1)]
        t = t[t.site_key.isin(eligible)]
        labels = sites.loc[eligible, "idr_label"]
        n_bound, n_unbound = int((labels == "bound").sum()), int((labels == "unbound").sum())
        if n_bound < 10 or n_unbound < 10:
            continue
        for cls in classes:
            u = t[t[column] == cls]
            draws = []
            for _ in range(N_DRAWS):
                picked = u.groupby("site_key", group_keys=False).sample(n=k, random_state=int(rng.integers(1 << 31)))
                sc = site_scores(picked)
                bound = (sites.loc[sc.index, "idr_label"] == "bound").to_numpy()
                draws.append({c: auroc(sc[c].to_numpy(), bound) for c in sc.columns})
            d = pd.DataFrame(draws)
            row = dict(tf=tf, cls=cls, k=k, n_bound=n_bound, n_unbound=n_unbound)
            row.update({c: d[c].mean() for c in d.columns})
            rows.append(row)
    return pd.DataFrame(rows)


say("\n== E. AUROC bound vs unbound from k reads per class (same sites for every class; mean of "
    f"{N_DRAWS} draws)")
for classes, k, column, label in [
    (["in_range", "high"], 5, "rate_cls", "all reads"),
    (["low", "in_range", "high"], 1, "rate_cls", "all reads"),
    (["in_range", "high", "high12"], 3, "rate_cls4", "all reads (high = 8-12%)"),
]:
    say(f"-- classes {classes}, k = {k}, {label}")
    say(matched_auroc(rs, column, classes, k, TFS).round(3).to_string(index=False))
prot = rs[rs.protected_in_motif]
say("-- reads with a protected call in the motif (TF vs other protection), classes in_range/high, k = 3")
say(matched_auroc(prot, "rate_cls", ["in_range", "high"], 3, TFS).round(3).to_string(index=False))


# F. Simulated filter.
def idr_metrics(sc: pd.DataFrame, columns: list[str]) -> dict:
    s = sites.loc[sc.index]
    bound = (s.idr_label == "bound").to_numpy()
    access = s.accessibility.to_numpy()
    known = np.isfinite(access)
    res = {}
    for c in columns:
        v = sc[c].to_numpy()
        res[f"{c}|auroc"] = auroc(v, bound)
        res[f"{c}|auroc_access"] = stratified_auroc(v[known], bound[known], access[known])
        w = bound & known & np.isfinite(v)
        res[f"{c}|rho_bound"] = partial_spearman(s.chip_pooled_sum.to_numpy()[w], v[w], access[w])
    return res


say("\n== F. Simulated 2-8% filter: site scores from in-range reads vs all reads vs random subsets of "
    "the same per-site size (sites with >= 3 in-range reads)")
labeled = rs[rs.idr_label != "ambiguous"]
cols = MODELS + ["ref_occupancy", "ref_protected_fraction"]
rows = []
for tf in TFS:
    t = labeled[labeled.tf == tf]
    n_in = t[t.rate_cls == "in_range"].groupby("site_key").size()
    keep = n_in.index[n_in >= 3]
    labels = sites.loc[keep, "idr_label"]
    if (labels == "bound").sum() < 10 or (labels == "unbound").sum() < 10:
        continue
    t = t[t.site_key.isin(keep)]
    full = idr_metrics(site_scores(t), cols)
    filt = idr_metrics(site_scores(t[t.rate_cls == "in_range"]), cols)
    rand = []
    for _ in range(N_DRAWS):
        picked = pd.concat([g.sample(n=int(n_in[s]), random_state=int(rng.integers(1 << 31)))
                            for s, g in t.groupby("site_key")])
        rand.append(idr_metrics(site_scores(picked), cols))
    rand = pd.DataFrame(rand).mean()
    for c in cols:
        for metric in ("auroc", "auroc_access", "rho_bound"):
            k = f"{c}|{metric}"
            rows.append(dict(tf=tf, score=c, metric=metric, n_sites=len(keep), all_reads=full[k],
                             in_range_only=filt[k], random_same_size=rand[k]))
F = pd.DataFrame(rows)
F["filter_minus_random"] = F.in_range_only - F.random_same_size
F.to_csv(D / "simulated_filter.tsv", sep="\t", index=False)
for metric in ("auroc_access", "rho_bound", "auroc"):
    say(f"-- {metric}: in-range only minus random same-size subset (positive = filtering helps beyond read count)")
    say(F[F.metric == metric].pivot(index="tf", columns="score", values="filter_minus_random").reindex(TFS).round(3).to_string())
    say(f"-- {metric}: all reads / in-range only")
    p = F[F.metric == metric].pivot(index="tf", columns="score", values="all_reads").reindex(TFS).round(2).astype(str) + " / " + \
        F[F.metric == metric].pivot(index="tf", columns="score", values="in_range_only").reindex(TFS).round(2).astype(str)
    say(p.to_string())

# G. Footprint reads.
say("\n== G. Mean P(TF) on footprint reads at IDR-bound sites, by molecule rate class")
fpb = rs[(rs.cls == "footprint") & (rs.idr_label == "bound") & rs.tf.isin(FOOTPRINT_TFS)]
say(fpb.groupby(["tf", "rate_cls4"])[MODELS].mean().round(3).assign(
    n=fpb.groupby(["tf", "rate_cls4"]).size()).to_string())
say("Share of footprint reads among footprint + open reads at bound vs unbound sites, by class:")
fo = rs[rs.cls.isin(["footprint", "open"]) & (rs.idr_label != "ambiguous") & rs.tf.isin(FOOTPRINT_TFS)]
say((fo.cls == "footprint").groupby([fo.tf, fo.idr_label, fo.rate_cls4]).mean().unstack().round(3).to_string())
out.close()
