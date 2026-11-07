"""Figure for CTCF partial motif protection (reads from partial_motif_positions.py).

    python plot_partial.py OUT_DIR
"""
import sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import binom

OUT = sys.argv[1]
c = pd.read_csv(f"{OUT}/calls.tsv.gz", sep="\t"); c = c[c.tf == "ctcf"]
reads = c.groupby(["site_key", "read"]).agg(p=("chip_pct", "first"), n=("n_calls", "first"), k=("n_prot", "first")).reset_index()
s = lambda v: 1 / (1 + np.exp(-v))
fig, ax = plt.subplots(1, 4, figsize=(19, 4.2))
for a, n in zip(ax[:2], (3, 4)):
    bg = reads[(reads.p <= .25) & (reads.n == n)].k.value_counts(normalize=True).reindex(range(n + 1), fill_value=0).to_numpy()
    obs = reads[(reads.p > .75) & (reads.n == n)].k.value_counts().reindex(range(n + 1), fill_value=0).to_numpy()
    fits = {}
    for name, extra in (("unbound + bound Binomial(n, q)", False), ("+ partial class (uniform rate)", True)):
        def mix(t):
            w = np.exp(np.r_[0, t[:2 if extra else 1]]); w /= w.sum()
            m = w[0] * bg + w[1] * binom.pmf(range(n + 1), n, s(t[-1]))
            return m + (w[2] / (n + 1) if extra else 0)
        r = minimize(lambda t: -(obs * np.log(mix(t) + 1e-12)).sum(), np.zeros(3 if extra else 2), method="Nelder-Mead", options=dict(maxiter=5000))
        fits[name] = mix(r.x) * obs.sum()
    x = np.arange(n + 1)
    a.bar(x, obs, color="lightgray", label="observed (top-ChIP-quartile CTCF)")
    for (name, f), mk in zip(fits.items(), ("o-", "s--")):
        a.plot(x, f, mk, label=name)
    a.set_xlabel(f"protected calls out of {n} in motif"); a.set_ylabel("non-nucleosomal reads")
    a.set_title(f"CTCF motifs with {n} calls"); a.legend(fontsize=7)
x = c[(c.chip_pct > .75)]
a = ax[2]
for label, sel, col in (("all protected", x.n_prot == x.n_calls, "k"),
                         ("partial (1..n-1 protected)", (x.n_prot > 0) & (x.n_prot < x.n_calls), "tab:red")):
    g = x[sel].groupby("rel_pos").protected.agg(["mean", "size"]); g = g[g["size"] >= 30]
    a.plot(g.index, g["mean"], "o-", c=col, label=label)
a.set_xlabel("call position in CTCF motif (5'->3' on motif strand)"); a.set_ylabel("P(call protected)")
a.set_ylim(-.02, 1.02); a.legend(fontsize=8); a.set_title("Which motif calls stay protected")
# example site with most partial reads
pr = reads[(reads.p > .75) & (reads.k > 0) & (reads.k < reads.n)]
site = pr.site_key.value_counts().index[0]
y = c[c.site_key == site]
mat = y.pivot_table(index="read", columns="rel_pos", values="protected")
mat = mat.loc[mat.sum(1).sort_values().index]
a = ax[3]
a.imshow(mat.to_numpy(), aspect="auto", cmap=matplotlib.colors.ListedColormap(["mediumseagreen", "k"]), interpolation="nearest")
a.set_xticks(range(mat.shape[1]), mat.columns); a.set_xlabel("motif position"); a.set_ylabel("non-nucleosomal reads")
a.set_title(f"{site.split(':', 2)[2]} (black protected)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{OUT}/ctcf_partial_protection.png", dpi=90)
print(site)
