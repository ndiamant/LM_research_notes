"""Accessibility of held-out calls inside MAP TF calls, by the held-out position's context efficiency (2026-11-24).

Reads heldout_calls.npz (heldout_inside_tf.py). O/E = observed m6A / sum of e (expected rate if accessible).
Relative accessibility A = O/E of a class / O/E of held-out calls inside MAP open calls, same model and
position-efficiency quintile (corrects residual miscalibration of e by context). Real protection: A near the
histone level at every efficiency; residual context bias: A rising towards 1 at efficient positions.
Classes: open, histone, TF on a strong-footprint motif, TF off motif with the call's mean context efficiency
low (< 0.523, context deciles 0-2 of robustness_checks) / mid / high (>= 0.631, deciles 6-9).
95% intervals: 200 bootstrap resamples of windows.
"""
import numpy as np, pandas as pd
d = np.load("heldout_calls.npz", allow_pickle=True)
D = pd.DataFrame({k: d[k] for k in d.files if k != "models"}); models = list(d["models"])
G = ["open", "TF", "other_long", "histone"]
D["cls"] = np.select(
    [D.group == 0, D.group == 3, (D.group == 1) & D.call_strong, (D.group == 1) & (D.call_ctx < 0.523),
     (D.group == 1) & (D.call_ctx >= 0.631), D.group == 1],
    ["open", "histone", "TF strong motif", "TF off, low-eff call", "TF off, high-eff call", "TF off, mid-eff call"], "other_long")
edges = np.quantile(D.ctx, [0.2, 0.4, 0.6, 0.8])
D["qbin"] = np.searchsorted(edges, D.ctx, side="right")
order = ["open", "TF strong motif", "TF off, high-eff call", "TF off, mid-eff call", "TF off, low-eff call", "histone"]
rng = np.random.default_rng(0)
rows = []
for mi, m in enumerate(models):
    x = D[D.model == mi]
    wins = np.unique(x.window)
    agg = x.groupby(["cls", "qbin", "window"]).agg(obs=("truth", "sum"), exp=("e", "sum"), n=("truth", "size")).reset_index()
    cube = {c: {q: agg[(agg.cls == c) & (agg.qbin == q)].set_index("window")[["obs", "exp", "n"]].reindex(wins, fill_value=0) for q in range(5)} for c in order}
    boot = rng.integers(0, len(wins), (200, len(wins)))
    for c in order:
        for q in range(5):
            t, o = cube[c][q], cube["open"][q]
            oe = t.obs.sum() / t.exp.sum(); oe_open = o.obs.sum() / o.exp.sum()
            bt = (t.obs.to_numpy()[boot].sum(1) / t.exp.to_numpy()[boot].sum(1)) / (o.obs.to_numpy()[boot].sum(1) / o.exp.to_numpy()[boot].sum(1))
            rows.append(dict(model=m, cls=c, qbin=q, n=int(t.n.sum()), m6a_rate=t.obs.sum() / max(t.n.sum(), 1), OE=oe,
                             A=oe / oe_open, A_lo=np.nanpercentile(bt, 2.5), A_hi=np.nanpercentile(bt, 97.5)))
R = pd.DataFrame(rows); R.to_csv("heldout_inside_tf.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print("position-efficiency quintile edges (ctx):", np.round(edges, 3))
for m in models:
    r = R[R.model == m]
    t = r.pivot(index="cls", columns="qbin", values="A").loc[order]
    lo = r.pivot(index="cls", columns="qbin", values="A_lo").loc[order]; hi = r.pivot(index="cls", columns="qbin", values="A_hi").loc[order]
    n = r.pivot(index="cls", columns="qbin", values="n").loc[order]
    cell = t.round(2).astype(str) + " [" + lo.round(2).astype(str) + "," + hi.round(2).astype(str) + "]"
    print(f"\n{m}: relative accessibility A of held-out calls (open = 1) by held-out position efficiency quintile (0 = least efficient)")
    print(cell.to_string()); print("n held-out calls:"); print(n.to_string())
    oe = r[r.cls == "open"].set_index("qbin")
    print("open O/E (calibration of e):", oe.OE.round(3).to_dict(), " open m6A rate:", oe.m6a_rate.round(3).to_dict())
