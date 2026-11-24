"""Summarize runs/*.json from bench_train_step.py into a table."""
import glob, json, re
from pathlib import Path
import numpy as np

B = Path("/scratch/users/diamant/2026-11-23_jax_container_bench/runs")
WARM = 20


def compile_log(log: Path) -> dict:
    out = {}
    for line in log.read_text(errors="ignore").splitlines():
        if "train_step_jitted" not in line:
            continue
        m = re.search(r"Finished (tracing \+ transforming|tracing|jaxpr to MLIR module conversion|XLA compilation)[^\d]* in ([\d.e+-]+) ?(sec|s)", line)
        if m:
            key = {"tracing + transforming": "trace", "tracing": "trace", "jaxpr to MLIR module conversion": "lower", "XLA compilation": "xla"}[m.group(1)]
            out[key] = out.get(key, 0.0) + float(m.group(2))
    return out


rows = {}
for f in sorted(B.glob("*_rep*.json")):
    r = json.loads(f.read_text())
    env, rep = f.stem.split("_rep")
    st = r["steps"][WARM:]
    step = np.array([s["step_s"] for s in st])
    data = np.array([s["data_s"] for s in st])
    disp = np.array([s["dispatch_s"] for s in st])
    tot = step + data
    pure = np.array([p["step_s"] for p in r["pure"][10:]])
    val = r["val"]
    full = [v["step_s"] for v in val[1:] if v["batch_size"] == val[0]["batch_size"]]
    rows[(env, int(rep))] = dict(
        r=r,
        first=r["steps"][0]["step_s"],
        med=np.median(step), q1=np.percentile(step, 25), q3=np.percentile(step, 75),
        data=np.median(data), disp=np.median(disp), tot=np.median(tot), n=len(step),
        val_first=val[0]["step_s"], val_med=np.median(full), val_last=val[-1]["step_s"],
        val_last_bs=val[-1]["batch_size"], nval=len(val),
        peak=r["train_peak_bytes"] / 2**30, pure_peak=r["pure_peak_bytes"] / 2**30,
        pure_med=np.median(pure), pure_pipe=r["pure_pipelined_s"],
        pure_lower=r["pure_lower_s"], pure_compile=r["pure_compile_s"],
        load=r["data_load_s"], build=r["model_build_s"],
        clog=compile_log(f.with_suffix(".log")),
    )

v = rows[next(iter(rows))]["r"]
for k, x in sorted(rows.items()): print(k, f"data load {x['load']:.0f}s, model build {x['build']:.1f}s")
print("train batches/epoch:", v["num_train_batches_per_epoch"], " steady steps per run:", rows[next(iter(rows))]["n"])
for env in ("baseline", "container"):
    print(env, rows[(env, 1)]["r"]["versions"])
print()
hdr = "env        rep first_s  [trace lower xla]   step_ms med [IQR]        data_ms dispatch_ms total_ms | pure_ms pipe_ms pure_compile_s | val_first_s val_ms last_s(bs) | peakGiB pure_peak"
print(hdr)
for (env, rep), x in sorted(rows.items()):
    c = x["clog"]
    print(f"{env:10s} {rep}  {x['first']:6.1f}  [{c.get('trace',0):4.1f} {c.get('lower',0):4.1f} {c.get('xla',0):5.1f}]   "
          f"{x['med']*1e3:6.1f} [{x['q1']*1e3:6.1f},{x['q3']*1e3:6.1f}]   {x['data']*1e3:5.1f}   {x['disp']*1e3:6.1f}    {x['tot']*1e3:6.1f} | "
          f"{x['pure_med']*1e3:6.1f}  {x['pure_pipe']*1e3:6.1f}  {x['pure_compile']:5.1f} | "
          f"{x['val_first']:5.1f}  {x['val_med']*1e3:6.1f}  {x['val_last']:5.1f}({x['val_last_bs']}) | {x['peak']:5.2f}  {x['pure_peak']:5.2f}")

print()
def agg(env, key):
    a = np.array([x[key] for (e, _), x in rows.items() if e == env])
    return a
for key, scale, unit in [("first", 1, "s"), ("med", 1e3, "ms"), ("tot", 1e3, "ms"), ("pure_pipe", 1e3, "ms"), ("val_med", 1e3, "ms"), ("val_first", 1, "s"), ("peak", 1, "GiB")]:
    b, c = agg("baseline", key) * scale, agg("container", key) * scale
    print(f"{key:10s} baseline {b.mean():8.2f} ± {b.std(ddof=1):.2f} {unit}   container {c.mean():8.2f} ± {c.std(ddof=1):.2f} {unit}   baseline/container = {b.mean()/c.mean():.3f}")

# Speedup of the steady-state step with a bootstrap over steps and runs.
rng = np.random.default_rng(0)
def steps_of(env):
    return [np.array([s["step_s"] for s in x["r"]["steps"][WARM:]]) for (e, _), x in rows.items() if e == env]
bs, cs = steps_of("baseline"), steps_of("container")
ratios = []
for _ in range(2000):
    bb = np.median(np.concatenate([b[rng.integers(0, len(b), len(b))] for b in [bs[i] for i in rng.integers(0, len(bs), len(bs))]]))
    cc = np.median(np.concatenate([c[rng.integers(0, len(c), len(c))] for c in [cs[i] for i in rng.integers(0, len(cs), len(cs))]]))
    ratios.append(bb / cc)
print(f"step speedup (baseline/container median, pooled), bootstrap 95% CI: {np.median(np.concatenate(bs))/np.median(np.concatenate(cs)):.3f} [{np.percentile(ratios,2.5):.3f}, {np.percentile(ratios,97.5):.3f}]")

# Loss agreement over the first 20 steps (same init, same batches).
print()
ref = rows[("baseline", 1)]["r"]["steps"]
for (env, rep), x in sorted(rows.items()):
    s = x["r"]["steps"]
    same_batches = all(a["batch"] == b["batch"] for a, b in zip(s[:25], ref[:25]))
    l = np.array([a["loss"] for a in s[:20]]); l0 = np.array([a["loss"] for a in ref[:20]])
    print(f"{env} rep{rep} vs baseline rep1: same batches={same_batches} init={x['r']['init_params']} "
          f"max rel loss diff (20 steps)={np.max(np.abs(l-l0)/np.abs(l0)):.2e}  step1 rel diff={abs(l[0]-l0[0])/abs(l0[0]):.2e}  "
          f"max rel diff steps 1-280={np.max(np.abs(np.array([a['loss'] for a in s])-np.array([a['loss'] for a in ref]))/np.abs(np.array([a['loss'] for a in ref]))):.2e}")
