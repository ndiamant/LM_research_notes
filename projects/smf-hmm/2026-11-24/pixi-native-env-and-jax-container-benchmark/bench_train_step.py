"""Time smf_hmm training steps for one run's config in the current environment.

Rebuilds the run's data module (with a region subset), model, objective,
optimizer and EMA exactly as smf_hmm.cli.train does, then times real jitted
train steps on real batches (the same batches in every environment) and
validation steps. Writes one JSON with per-step timings, losses and versions.

Initial trainable params are shared across environments through --init-params:
the first run saves them, later runs load them, so losses start from the same
init even if the PRNG or Flax's Rngs differ between versions.

Never evaluates the test chromosomes (no test_batches call).
"""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from dataclasses import asdict
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import optax
from flax import nnx
from omegaconf import OmegaConf

import smf_hmm
from smf_hmm import train
from smf_hmm.cli import train as cli
from smf_hmm.data import ChromosomeDataModule, ChromosomeDataModuleFromH5Config, SMFItem


def section(cfg, name: str) -> dict:
    d = OmegaConf.to_container(cfg[name], resolve=True)
    d.pop("_target_", None)
    return d


def flat_params(model: nnx.Module):
    state = nnx.state(model, nnx.Param)
    leaves_with_path, treedef = jax.tree_util.tree_flatten_with_path(state)
    keys = [jax.tree_util.keystr(p) for p, _ in leaves_with_path]
    return state, treedef, keys, [leaf for _, leaf in leaves_with_path]


def save_or_load_init(model: nnx.Module, path: Path) -> str:
    state, treedef, keys, leaves = flat_params(model)
    if not path.exists():
        np.savez(path, keys=np.array(keys), *[np.asarray(x) for x in leaves])
        return "saved"
    data = np.load(path)
    saved_keys = list(data["keys"])
    saved = [data[f"arr_{i}"] for i in range(len(saved_keys))]
    if saved_keys == keys:
        how = "loaded by key"
    else:
        # Path formatting can differ between Flax versions; fall back to order.
        assert len(saved) == len(leaves), (len(saved), len(leaves))
        how = "loaded by order (keys differ)"
    for s, leaf in zip(saved, leaves):
        assert s.shape == leaf.shape and s.dtype == leaf.dtype, (s.shape, leaf.shape)
    new_leaves = [jnp.asarray(s) for s in saved]
    nnx.update(model, jax.tree_util.tree_unflatten(treedef, new_leaves))
    return how


def batch_hash(batch) -> str:
    h = hashlib.sha1()
    h.update(np.ascontiguousarray(batch.smf).tobytes())
    h.update(np.ascontiguousarray(batch.dna).tobytes())
    return h.hexdigest()[:16]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--init-params", required=True)
    parser.add_argument("--region-subset", type=int, default=40000)
    parser.add_argument("--steps", type=int, default=280)
    parser.add_argument("--val-batches", type=int, default=1000)
    parser.add_argument("--pure-steps", type=int, default=200)
    args = parser.parse_args()

    jax.config.update("jax_enable_compilation_cache", False)
    # As smf_hmm.cli.train.main.
    jax.config.update("jax_default_matmul_precision", "float32")

    cfg = OmegaConf.load(Path(args.run_dir) / ".hydra" / "config.yaml")
    data_cfg = cli.DataConfig(**section(cfg, "data"))
    data_cfg.region_subset = args.region_subset
    data_cfg.verbose = False
    transition_cfg = cli.SingleTFTaperedHistoneFixedBlockTransitionCLIConfig(
        **section(cfg, "transition")
    )
    prior_cfg = cli.PriorConfig(**section(cfg, "prior"))
    emission_cfg = cli.EmissionConfig(**section(cfg, "emission"))
    opt_cfg = cli.OptimizerConfig(**section(cfg, "optimizer"))
    obj_cfg = cli.ObjectiveConfig(**section(cfg, "objective"))
    trainer_cfg = cli.TrainerConfig(**section(cfg, "trainer"))

    t = time.perf_counter()
    data_module = ChromosomeDataModule.from_h5(
        ChromosomeDataModuleFromH5Config(**asdict(data_cfg))
    )
    data_load_s = time.perf_counter() - t

    t = time.perf_counter()
    model = cli._build_hmm_model(
        transition_module=transition_cfg.build(),
        prior=prior_cfg,
        emission=emission_cfg,
        assay_type=data_cfg.assay_type,
    )
    init_how = save_or_load_init(model, Path(args.init_params))
    objective_fn = train.build_smf_read_objective(
        entropy_regularization=train.EntropyRegularizationConfig(
            weight=obj_cfg.entropy_weight
        ),
        prior_occupancy_kl_regularization=train.PriorOccupancyKLRegularizationConfig(
            weight=obj_cfg.prior_occupancy_kl_weight,
            group_labels=obj_cfg.prior_occupancy_group_labels,
            target_rates=obj_cfg.prior_occupancy_target_rates,
        ),
    )
    trainer = train.Trainer(
        model=model,
        optimizer=nnx.Optimizer(
            model,
            optax.chain(
                optax.clip_by_global_norm(opt_cfg.clip_norm),
                optax.adam(opt_cfg.lr),
            ),
            wrt=nnx.Param,
        ),
        objective_fn=objective_fn,
        config=train.TrainerConfig(
            num_epochs=1,
            print_every=None,
            verbose=False,
            target_dtype=cli._dtype_from_name(trainer_cfg.target_dtype),
            checkpoint_dir=None,
            ema_decay=opt_cfg.ema_decay,
        ),
    )
    model_build_s = time.perf_counter() - t

    train_iter = data_module.train_batches(0)
    n_steps = min(args.steps, len(train_iter))
    steps = []
    staged = []
    for i in range(n_steps):
        t0 = time.perf_counter()
        batch = next(train_iter)
        t1 = time.perf_counter()
        # Trainer.train_step, split at the dispatch and at the host sync.
        loss, metrics = trainer._train_step(trainer.model, trainer.optimizer, batch)
        t2 = time.perf_counter()
        trainer._ema_update(trainer.eval_model, trainer.model)
        loss = float(loss)
        metrics = {k: float(v) for k, v in metrics.items()}
        t3 = time.perf_counter()
        steps.append(
            dict(
                data_s=t1 - t0,
                dispatch_s=t2 - t1,
                step_s=t3 - t1,
                loss=loss,
                nll=metrics["nll"],
                batch=batch_hash(batch) if i < 25 else None,
            )
        )
        if len(staged) < 8:
            staged.append(batch)
    # Drain queued EMA work before reading memory and timing validation.
    jax.block_until_ready(nnx.state(trainer.eval_model, nnx.Param))
    train_peak = jax.devices()[0].memory_stats().get("peak_bytes_in_use")

    val = []
    val_iter = data_module.validation_batches(0)
    for i, batch in enumerate(val_iter):
        if i >= args.val_batches:
            break
        t0 = time.perf_counter()
        loss, _ = trainer.eval_step(batch)
        val.append(
            dict(
                step_s=time.perf_counter() - t0,
                loss=loss,
                batch_size=int(np.shape(batch.smf)[0]),
            )
        )
    stats = jax.devices()[0].memory_stats()

    # Pure-JAX train step on device-resident batches: same objective and
    # optimizer update, without nnx.jit's per-call graph traversal or host
    # batch transfer, so its time is (close to) the device time of one step.
    import flax

    graphdef, state = nnx.split((trainer.model, trainer.optimizer))

    def pure_step(state, smf, dna, avg, efficiency):
        model, optimizer = nnx.merge(graphdef, state)
        batch = SMFItem(smf=smf, dna=dna, avg=avg, efficiency=efficiency)
        (loss, _), grads = nnx.value_and_grad(
            lambda m: objective_fn(m, batch), has_aux=True
        )(model)
        if flax.__version__.startswith("0.10"):
            optimizer.update(grads)
        else:
            optimizer.update(model, grads)
        return nnx.state((model, optimizer)), loss

    dev_batches = [
        tuple(
            jax.device_put(np.asarray(x))
            for x in (b.smf, b.dna, b.avg, train._efficiency_or_zeros(b))
        )
        for b in staged
    ]
    t = time.perf_counter()
    lowered = jax.jit(pure_step).lower(state, *dev_batches[0])
    pure_lower_s = time.perf_counter() - t
    t = time.perf_counter()
    compiled = lowered.compile()
    pure_compile_s = time.perf_counter() - t
    pure = []
    for i in range(args.pure_steps + 10):
        t0 = time.perf_counter()
        state, loss = compiled(state, *dev_batches[i % len(dev_batches)])
        t1 = time.perf_counter()
        loss.block_until_ready()
        t2 = time.perf_counter()
        pure.append(dict(dispatch_s=t1 - t0, step_s=t2 - t0))
    # Back-to-back steps without a per-step sync: the device-bound throughput.
    t = time.perf_counter()
    for i in range(args.pure_steps):
        state, loss = compiled(state, *dev_batches[i % len(dev_batches)])
    loss.block_until_ready()
    pure_pipelined_s = (time.perf_counter() - t) / args.pure_steps
    pure_peak = jax.devices()[0].memory_stats().get("peak_bytes_in_use")

    smi = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"],
        capture_output=True,
        text=True,
    ).stdout.strip()
    import jaxlib.version
    import orbax.checkpoint

    result = dict(
        args=vars(args),
        versions=dict(
            python=platform.python_version(),
            jax=jax.__version__,
            jaxlib=jaxlib.version.__version__,
            flax=flax.__version__,
            optax=optax.__version__,
            orbax=getattr(orbax.checkpoint, "__version__", None),
            numpy=np.__version__,
            glibc=" ".join(platform.libc_ver()),
            smf_hmm=list(smf_hmm.__path__),
            gpu=smi,
            device=str(jax.devices()[0].device_kind),
        ),
        hostname=platform.node(),
        init_params=init_how,
        num_train_batches_per_epoch=len(train_iter),
        data_load_s=data_load_s,
        model_build_s=model_build_s,
        steps=steps,
        val=val,
        train_peak_bytes=train_peak,
        final_peak_bytes=stats.get("peak_bytes_in_use"),
        bytes_limit=stats.get("bytes_limit"),
        pure_lower_s=pure_lower_s,
        pure_compile_s=pure_compile_s,
        pure=pure,
        pure_pipelined_s=pure_pipelined_s,
        pure_peak_bytes=pure_peak,
    )
    Path(args.out).write_text(json.dumps(result, indent=1))
    s = np.array([x["step_s"] for x in steps[min(20, len(steps) // 2) :]])
    p = np.array([x["step_s"] for x in pure[10:]])
    print(
        f"{jax.__version__} first={steps[0]['step_s']:.1f}s "
        f"median={np.median(s)*1e3:.1f}ms IQR=[{np.percentile(s,25)*1e3:.1f},"
        f"{np.percentile(s,75)*1e3:.1f}] n={len(s)} "
        f"peak={train_peak/2**30 if train_peak else float('nan'):.2f}GiB "
        f"init={init_how} pure: lower={pure_lower_s:.1f}s compile={pure_compile_s:.1f}s "
        f"median={np.median(p)*1e3:.1f}ms pipelined={pure_pipelined_s*1e3:.1f}ms"
    )


if __name__ == "__main__":
    main()
