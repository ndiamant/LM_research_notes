---
title: "Native pixi env for smf_hmm (JAX capped at 0.4.35) and a JAX 0.12 container benchmark: the container is 5% slower"
date: 2026-11-24
project: smf-hmm
agent: Claude Code (claude-opus-5-5), benchmark run by a Claude Code subagent
status: draft
sources:
  - smf_hmm branch nd/pixi-migration (48e9f27)
  - smf_hmm branch nd/jax-container-bench (80aa05c, b91e7af)
  - /scratch/users/diamant/2026-11-23_jax_container_bench (benchmark outputs; scripts copied into this note's folder)
  - https://api.anaconda.org/package/conda-forge/<package>/files (build metadata)
tags: [environment, pixi, jax, flax, orbax, checkpointing, apptainer, benchmark, sherlock]
---

# Summary

This note follows up the [pixi/Orbax planning note](../2026-11-19/pixi-migration-and-orbax-v1-checkpointing.md) and **corrects its main claim**. On Sherlock (CentOS 7, glibc 2.17) the newest conda-forge GPU JAX is **jax/jaxlib 0.4.35, the version the old `hmm` conda env already had**, not 0.11. Each jaxlib release requires a minimum cuDNN, and cuDNN 9.3 is the newest conda-forge build for glibc 2.17; only jaxlib ≤ 0.4.35 accepts it. Newer cuBLAS/cuSOLVER builds have the same glibc problem.

Two things followed:

1. **Native pixi env (`nd/pixi-migration`, 48e9f27):** jax 0.4.35, flax 0.10.4, **orbax 0.9.0** (up from 0.4.4), plus the BioPython/BigWig libraries, so the `smf_clean` env is no longer needed. Checkpoints were rewritten as one orbax composite directory (state plus JSON spec). All 388 core tests pass on an H100.
2. **JAX 0.12 in an Apptainer container (`nd/jax-container-bench`):** on the most recent AlphaGenome-prior training config, the container is **about 5% slower per step**, with a speedup of 0.949× (95% CI 0.944–0.951). It is also about 10 s slower to compile and uses 15% less GPU memory. **Decision: stay on the native stack.** The ~8% rule of thumb for added complexity is not even reached in the other direction.

# Key Points

- **Why the native GPU stack is capped (glibc 2.17):**

  | chain link | newest glibc-2.17 conda-forge build | consequence |
  |---|---|---|
  | cuDNN | 9.3.0.75 | jaxlib ≤ 0.4.35; 0.4.36 and later need cuDNN ≥ 9.7 |
  | cuBLAS / cuSOLVER | CUDA 12.6 | jaxlib CUDA 12.9 builds (0.6.2 and later) are excluded |
  | jaxlib 0.6.2 CUDA 12.6 build | labelled `broken` on conda-forge | |
  | tensorstore ≥ 0.1.84, needed by orbax 0.12 | built against a newer abseil than jaxlib 0.4.35–0.6 | orbax capped at 0.9.0 |
  | PyPI jaxlib / tensorstore wheels | manylinux_2_27 | cannot be installed on glibc 2.17 |

  The CPU-only env could take jax 0.10 and flax 0.12. It is pinned to the GPU versions anyway, so the code targets a single Flax API.
- **Container benchmark** (same H100 on sh04-03n07, driver 550.163, 3 alternating runs per env):

  | | native | container |
  |---|---|---|
  | versions | jax/jaxlib 0.4.35 (CUDA 12.6, cuDNN 9.3), flax 0.10.4, orbax 0.9.0, glibc 2.17 | jax/jaxlib 0.12.0.1 (CUDA 12.9, cuDNN 9.10), flax 0.12.12, orbax 0.12.12, glibc 2.34 (Rocky 9) |
  | first train step (compile + run) | 25.7 ± 3.0 s | 35.3 ± 2.5 s |
  | train step, median | **160.5 ± 0.6 ms** | **169.3 ± 0.2 ms** |
  | GPU-only step (batches already on GPU) | 155.8 ms | 164.6 ms |
  | host data loading per batch | about 20 ms | about 21 ms |
  | validation step | 95.8 ms | 95.1 ms |
  | peak GPU memory (training) | 15.7 GiB | 13.3 GiB |

- **Step time is dominated by the GPU in both envs.** The host-side overhead of `nnx.jit` is only about 4 ms per step, so the gap lies in the compiled computation itself.
- **CUDA 12.9 ran on driver 550 (CUDA 12.4) without `cuda-compat`**, through CUDA minor-version compatibility.
- **Untested caveat:** jaxlib 0.12.1 (cuDNN 9.28) could not be solved, because its cuDNN needs the `__cuda_arch` virtual package that Sherlock's pixi 0.53 does not provide. So this benchmark cannot rule out that the very latest build closes the gap. Nothing suggests a large win.
- **Flax 0.11+ behaviour change, relevant whenever we upgrade:** the optimizer no longer holds a reference to the model. The trainer's dtype cast went through the optimizer, so it **would silently stop casting the model's parameters**. The port in `b91e7af` casts model and optimizer state separately.

# Details

## Native pixi env and checkpoint rewrite (`nd/pixi-migration`, 48e9f27)

- **Config:** `pyproject.toml` has a `[tool.pixi]` workspace with conda-forge and bioconda, `linux-64` only, `libc = "2.17"`, and `jax == 0.4.35`, `flax == 0.10.4`, `orbax-checkpoint 0.9.*`, `numpy 2.4.*`.
  - `gpu` feature: jaxlib `*cuda126*`. `cpu` feature: jaxlib `*cpu*`.
  - Environments `default` (GPU) and `cpu`.
  - pybigtools and pysam come from bioconda.
  - macOS ARM was dropped: jaxlib 0.4.35/0.6 and tensorstore ≥ 0.1.84 conflict over abseil on osx-arm64.
- **Location:** the env is detached to `$GROUP_HOME/pixi-envs` through `.pixi/config.toml`. `.pixi` is gitignored, so every checkout or worktree must run `pixi config set --local detached-environments ...` again. Installed size is 5.0 GB.
- **pytest:** `[tool.pytest.ini_options] pythonpath = ["."]` was added so `pixi run test` resolves `baselines/` and `evaluation/`.
- **Checkpoints:** `ocp.Checkpointer(ocp.CompositeCheckpointHandler())` saves `Composite(state=StandardSave(nnx.state(model)), spec=JsonSave({format_version: 2, model_spec, meta}))` with `force=True`.
  - Loading restores `spec` only, rebuilds the model, then restores `state` against `jax.tree.map(ocp.utils.to_shape_dtype_struct, nnx.state(model))`, so the arrays land on the current device.
  - A checkpoint saved on the H100 restored in a `JAX_PLATFORMS=cpu` process with the same log-probability. That made the `HostCheckpointer` workaround in `evaluation/emission_context_bias/` unnecessary.
  - Removed: the `best.model_spec.json` sidecar, the manual delete, the `Checkpointer` protocol, the pickle test fake and `make_hmm_checkpoint_payload`. `checkpointing.py` went from 509 to 464 lines.
  - Most of what remains (~290 lines) is the per-module spec registry. A config-driven replacement (rebuild from the run's hydra-zen config) was proposed but is not done.
- **Orbax 0.9 behaviour change:** `StandardCheckpointer` saves **asynchronously**, so load-after-save failed in two tests with the old code. The rewrite uses the synchronous `ocp.Checkpointer`.
- **Test fixes:**
  - two subprocess tests had a hard-coded Mac interpreter path; they now use `sys.executable`;
  - the child processes of CLI tests run with `JAX_PLATFORMS=cpu`, because Sherlock GPUs are in `Exclusive_Process` mode and pytest already holds the device;
  - one float32 equality test (`assert_allclose` at the default `rtol=1e-7`, below float32 precision) now uses `rtol=1e-6`.
- **Not migrated yet:** the `tf_chip` evaluation code, docs and sbatch scripts still point at the conda envs. No full training run has been done in the pixi env yet.

## Container benchmark (`nd/jax-container-bench`)

- **Image:** `/home/groups/btrippe/diamant/containers/rockylinux9.sif` (72 MB, Rocky Linux 9, glibc 2.34). The pixi env was solved with `libc = "2.34"` and a `*cuda129*` jaxlib (`80aa05c`), and stored detached in `$GROUP_HOME/pixi-envs`. It was 6.7 GB and is being deleted, since the decision is not to use it.
- **Flax 0.12 port** (`b91e7af`, 4 files, +35/−33):
  - `optimizer.update(model, grads)`;
  - separate model and optimizer dtype casts;
  - `nnx.List` for `blocks`, `dilated_residual_convs`, `norms` and `adapters` in `prior_architectures.py`;
  - optional submodules assigned once, as a module or `None`;
  - the AlphaGenome weights dict wrapped in `nnx.data(...)`, with reads changed from `.value` to `[...]`.

  Eight core test files (262 tests) pass in the container. The branch is **kept** as the starting point for an upgrade once Sherlock's OS moves past glibc 2.17.
- **Benchmark protocol** (`bench_train_step.py`):
  - Rebuilds the model, objective, optimizer and EMA of `/scratch/users/diamant/2026-11-17_fiber_hist190_tf60_ol61_kl300_alphagenome_bf16` from its `.hydra/config.yaml`. That config is the AlphaGenome prior (3 blocks, bf16), the fixed-block tapered-histone transition with TF 15–60 and other_long 61–100, and batch 128 × 32 reads × 512 bp.
  - Uses `data.region_subset=40000`, giving 221 train batches.
  - Times 280 train steps, taking the median over steps 21–221, and 200 GPU-only steps.
  - Initial params are shared across envs (saved once to `init.npz`) and batches are hash-verified identical. wandb is off and the test chromosomes are never touched.
  - Runs alternate between envs on one GPU, 3 per env.
- **Loss agreement:** the step-1 loss differs by about 3e-7 relative. The maximum over the first 20 steps is 7e-4 across envs, against 1e-5 between native runs. That is plausibly different bf16 numerics in the newer XLA/cuDNN. Over hundreds of steps, runs diverge by similar amounts within and across envs (12–19%), as expected for chaotic training trajectories. **Inferred, not verified:** the 20-step gap is numerical noise rather than a semantic change.
- **Test suite:** it runs faster in the container (310 s against 591 s), but that suite is dominated by compile time, not training throughput.

Commands (scripts copied to [`pixi-native-env-and-jax-container-benchmark/`](pixi-native-env-and-jax-container-benchmark/); outputs in `/scratch/users/diamant/2026-11-23_jax_container_bench/runs/`, with the full table in [`summary.txt`](pixi-native-env-and-jax-container-benchmark/summary.txt)):

```bash
module load devel pixi/0.53.0; export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
B=/scratch/users/diamant/2026-11-23_jax_container_bench
# Container env, once, from the jax-container-bench worktree (needs the deleted env rebuilt):
pixi config set --local detached-environments /home/groups/btrippe/diamant/pixi-envs
$B/ctr.sh pixi install        # ctr.sh = apptainer exec --nv with Sherlock binds (see copy in note folder)
$B/run_all.sh                 # 3 alternating runs per env, about 25 min on one H100
<native env python> $B/analyze.py   # writes runs/summary.txt
```

Native env setup and tests (from `.claude/worktrees/pixi-migration`, inside a GPU job):

```bash
module load devel pixi/0.53.0; export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
pixi config set --local detached-environments /home/groups/btrippe/diamant/pixi-envs
pixi lock && pixi install -e default
MPLCONFIGDIR=$PWD/.mplconfig pixi run pytest -q -p no:cacheprovider tests/test_emission_modules.py ...  # 388 core tests
```

How the glibc caps were found: build metadata from `curl -s https://api.anaconda.org/package/conda-forge/<pkg>/files` (fields `attrs.depends`, `labels`) for jaxlib, cudnn, libcublas, libcusolver, tensorstore and orbax-checkpoint, plus the `pixi lock` solver errors.

# Related Notes

- [Pixi migration and Orbax v1 checkpointing (planning)](../2026-11-19/pixi-migration-and-orbax-v1-checkpointing.md): the plan this note follows up. Its claim that GPU jaxlib up to 0.11 runs on glibc 2.17 is **wrong**, because it checked only jaxlib's own glibc requirement and not its CUDA library dependencies. The Orbax v1 sketch there needs orbax 0.12, which is unavailable natively; the orbax 0.9 composite design above replaces it.
- [Fiber-seq model components report](../2026-11-18/fiberseq-model-components.md): describes the AlphaGenome-prior model family whose training config was benchmarked.

# Open Questions

- Would jaxlib 0.12.1 with cuDNN 9.28 (needs a pixi newer than Sherlock's 0.53) change the result? Probably not by much, but it is untested.
- Should checkpoints switch to rebuilding from the run's hydra-zen config (Lightning `save_hyperparameters` style), removing the ~290-line spec registry? This was proposed and the user agreed in principle; not implemented.
- `nd/pixi-migration` still needs the `tf_chip` scripts and docs migrated, plus one real training run, before merging.
- Is the 7e-4 20-step loss gap across envs purely numerical? It would matter only if the upgrade is revisited.

# Sources

- smf_hmm `nd/pixi-migration` @ 48e9f27: `pyproject.toml`, `pixi.lock`, `smf_hmm/checkpointing.py`, `AGENTS.md`, tests.
- smf_hmm `nd/jax-container-bench` @ 80aa05c (container pixi solve) and b91e7af (Flax 0.12 port).
- Benchmark scripts (copied): [`bench_train_step.py`](pixi-native-env-and-jax-container-benchmark/bench_train_step.py), [`run_all.sh`](pixi-native-env-and-jax-container-benchmark/run_all.sh), [`ctr.sh`](pixi-native-env-and-jax-container-benchmark/ctr.sh), [`analyze.py`](pixi-native-env-and-jax-container-benchmark/analyze.py), [`summary.txt`](pixi-native-env-and-jax-container-benchmark/summary.txt).
- Reference training run: `/scratch/users/diamant/2026-11-17_fiber_hist190_tf60_ol61_kl300_alphagenome_bf16` (`.hydra/config.yaml`, `train.sbatch`).
- conda-forge build metadata: https://api.anaconda.org/package/conda-forge/jaxlib/files (and cudnn, libcublas, libcusolver, tensorstore, orbax-checkpoint).
- Flax 0.11/0.12 release notes: https://github.com/google/flax/releases; JAX changelog: https://github.com/jax-ml/jax/blob/main/CHANGELOG.md
