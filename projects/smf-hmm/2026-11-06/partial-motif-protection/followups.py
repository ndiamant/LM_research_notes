"""Follow-up tables for the partial motif protection note.

    PYTHONPATH=<smf_hmm repo> python followups.py OUT_DIR CHR8_DATASET_DIR FOOTPRINT_CHR8_DATASET_DIR

OUT_DIR must hold calls.tsv.gz from partial_motif_positions.py. Prints:
  1. CTCF read-class mixture fits (n = 3, 4 motif calls)
  2. per-site orientation consistency of one-of-two-protected CTCF reads
  3. OCT4 / NF-Y / CTCF read classes with the motif widened by 0, 3, 5 bp
  4. callable sites per motif: SMF (CGC-merged calls in data) vs m6A (A/T bp)
"""
import sys
import h5py, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import binom
from evaluation.tf_chip.analyze_covariation import read_table
from smf_hmm.data import merge_cgc_sparse

OUT, CHR8, FP = sys.argv[1:4]
pd.set_option("display.width", 250)
sig = lambda v: 1 / (1 + np.exp(-v))
c = pd.read_csv(f"{OUT}/calls.tsv.gz", sep="\t"); c = c[c.tf == "ctcf"]
reads = c.groupby(["site_key", "read"]).agg(p=("chip_pct", "first"), n=("n_calls", "first"), k=("n_prot", "first")).reset_index()

print("1. CTCF mixture fits, top-ChIP-quartile non-nucleosomal reads")
for n in (3, 4):
    bg = reads[(reads.p <= .25) & (reads.n == n)].k.value_counts(normalize=True).reindex(range(n + 1), fill_value=0).to_numpy()
    obs = reads[(reads.p > .75) & (reads.n == n)].k.value_counts().reindex(range(n + 1), fill_value=0).to_numpy()
    for extra in (False, True):
        def mix(t):
            w = np.exp(np.r_[0, t[:2 if extra else 1]]); w /= w.sum()
            return w[0] * bg + w[1] * binom.pmf(range(n + 1), n, sig(t[-1])) + (w[2] / (n + 1) if extra else 0), w
        r = minimize(lambda t: -(obs * np.log(mix(t)[0] + 1e-12)).sum(), np.zeros(3 if extra else 2),
                     method="Nelder-Mead", options=dict(maxiter=5000))
        m, w = mix(r.x)
        name = "bg + bound Binomial + partial(uniform)" if extra else "bg + bound Binomial"
        print(f"  n={n} obs={obs.tolist()} {name}: fit={np.round(m * obs.sum()).astype(int).tolist()} "
              f"weights={np.round(w, 2).tolist()} q={sig(r.x[-1]):.2f} NLL={r.fun:.1f}")

print("\n2. Orientation consistency, CTCF 2-call motifs, reads with exactly 1 of 2 protected")
one = c[(c.n_calls == 2) & (c.n_prot == 1) & (c.protected == 1)]
rows = []
for key, g in one.groupby("site_key"):
    if len(g) >= 8:
        vc = g.rel_pos.value_counts()
        rows.append(dict(chip_pct=g.chip_pct.iloc[0], n=len(g), dominant=vc.iloc[0] / len(g), pos=vc.index[0]))
t = pd.DataFrame(rows)
rng = np.random.default_rng(0)
for lab, sel in (("Q4", t.chip_pct > .75), ("Q1-3", t.chip_pct <= .75)):
    x = t[sel]
    null = np.median([max(b, n - b) / n for n in x.n for b in [rng.binomial(n, .5)]])
    print(f"  {lab}: sites={len(x)} median dominant share={x.dominant.median():.2f} share>=0.9={(x.dominant >= .9).mean():.2f} "
          f"null median={null:.2f}")
print("  dominant protected motif position (Q4):", t[t.chip_pct > .75].pos.value_counts().head(8).to_dict())

print("\n3./4. Widened motif windows and callable sites")
res = []
for D, tfs in ((CHR8, ["oct4", "ctcf"]), (FP, ["nfy"])):
    u = pd.read_csv(f"{D}/analysis_clip_ema/unique_motif_summaries.tsv.gz", sep="\t"); u = u[u.tf.isin(tfs)].copy()
    u["q"] = u.groupby("tf").chip_pooled_sum.rank(pct=True).gt(.75).map({True: "Q4", False: "Q1-3"})
    with h5py.File(f"{D}/dataset.h5") as f:
        win = read_table(f["window_metadata"]).set_index("window_id")
        for wid, rs in u.groupby("window_id"):
            w = f["windows"][wid]; dna = w["dna"][:].astype(np.int64)
            obs, pos = merge_cgc_sparse(dna, w["observations"][:], w["positions"][:])
            ws = int(win.loc[wid, "window_start"])
            for r in rs.itertuples():
                l0, r0 = r.start - ws, r.end - ws
                for ext in (0, 3, 5):
                    l, rr = l0 - ext, r0 + ext
                    if l < 1 or rr > len(dna) - 1:
                        continue
                    inm = (pos >= l) & (pos < rr)
                    sides = [obs[:, (pos >= lo) & (pos < hi)] for lo, hi in ((l - 25, l - 5), (rr + 5, rr + 25))]
                    nonnuc = np.all([(s >= 0).any(1) & ~(s == 1).any(1) for s in sides], axis=0)
                    m = obs[:, inm]; no = (m >= 0).sum(1); k = (m == 1).sum(1); keep = nonnuc & (no >= 2)
                    res.append(dict(tf=r.tf, q=r.q, ext=ext, n_calls=int(inm.sum()), m6a=int(np.isin(dna[l:rr], [0, 3]).sum()),
                                    reads=int(keep.sum()), mixed=int(((k > 0) & (k < no))[keep].sum()), prot=int((k == no)[keep].sum())))
t = pd.DataFrame(res)
g = t.groupby(["tf", "ext", "q"]).agg(sites=("n_calls", "size"), sites_2calls=("n_calls", lambda x: (x >= 2).sum()),
                                      reads=("reads", "sum"), mixed=("mixed", "sum"), all_prot=("prot", "sum"))
g["mixed_frac"] = g.mixed / g.reads; g["all_prot_frac"] = g.all_prot / g.reads
print(g.round(3).to_string())
h = t[t.ext.isin([0, 5])].groupby(["tf", "ext"]).agg(smf_median=("n_calls", "median"), smf_mean=("n_calls", "mean"),
    smf_ge2=("n_calls", lambda x: (x >= 2).mean()), m6a_median=("m6a", "median"), m6a_mean=("m6a", "mean"))
print(h.round(2).to_string())
