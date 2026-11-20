---
title: "Pixi migration and Orbax v1 checkpointing for smf_hmm: feasible, not yet done"
date: 2026-11-19
project: smf-hmm
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - ~/repos/SMFNet_public/pyproject.toml (pixi setup used as the template)
  - https://api.anaconda.org/package/<channel>/<package>/files (conda-forge/bioconda build metadata)
  - https://github.com/google/orbax (checkpoint/ at tag v0.12.8)
  - https://github.com/google/flax/blob/main/docs_nnx/guides/checkpointing.ipynb
tags: [environment, pixi, jax, flax, orbax, checkpointing, sherlock, planning]
---

# Summary

We looked into moving smf_hmm from two conda envs (`hmm` for JAX, `smf_clean` for BioPython and BigWig) to one pixi env, modeled on SMFNet_public, and into replacing `smf_hmm/checkpointing.py` with Orbax's new v1 API. Both are feasible. **Neither has been implemented.** As of 2026-11-19 the `hmm` env is unchanged (jax/jaxlib 0.4.35 CUDA 12.6, flax 0.10.4, orbax-checkpoint 0.4.4), and `pyproject.toml` has no pixi section.

- **The JAX version is limited by the GPU driver, not by Sherlock's glibc (2.17).** conda-forge has glibc-2.17 CUDA 12 jaxlib builds up to 0.11.2. But jaxlib 0.7 and later needs `cuda-version >=12.9`, and the H100 node checked (sh04-03n07) has driver 550.163, which supports CUDA 12.4. The current env already runs CUDA 12.6 on that driver. That probably works through CUDA minor-version compatibility, with XLA turning off parallel compilation. This is **untested** for 12.9.
- **One env is enough.** `hmm` already contains biopython 1.86 and pyBigWig. Merging only needs `pybigtools` (now only in `smf_clean`) and `pysam` (in neither env; only used by `20260420_create_h5_from_mBED.py`).
- **Orbax v1 (`orbax.checkpoint.v1`, still under `experimental`) would simplify checkpointing.** `custom_metadata` can hold the model spec inside the checkpoint, saves are atomic, and loads are checked against the model's expected structure.

# Key Points

- **Versions with glibc-2.17 builds (checked 2026-10-21):**

  | package | in `hmm` now | newest with a glibc-2.17 build |
  |---|---|---|
  | jax / jaxlib CUDA 12 | 0.4.35 | 0.11.2 (CUDA 13 builds need glibc 2.28) |
  | flax | 0.10.4 | 0.12.8 (needs jax ≥ 0.8.1) |
  | orbax-checkpoint | 0.4.4 | 0.12.8 |
  | optax | 0.2.7 | 0.2.8 |
  | biopython / pybigwig / pybigtools / pysam | – | 1.88 / 0.3.26 / 0.3.0 / 0.24.1 |

- **Fallback if CUDA 12.9 fails on driver 550:** jaxlib 0.6.2 (CUDA 12.6, the same CUDA version as now). That caps flax at about 0.11.
- **Known Flax 0.12 break:** `optimizer.update(grads)` in `smf_hmm/train.py` (line 545 on 2026-11-18, `1d86e57`) must become `optimizer.update(model, grads)`. Every `nnx.Optimizer(...)` call already passes `wrt=nnx.Param`. Other NNX changes may surface only in tests.
- **Storage:** the env must be a detached env in `$GROUP_HOME`, because `$HOME` was 89% full (13.4 of 15 GB). `$GROUP_HOME` was 94% full (973 GB of 1 TB). A CUDA JAX env is roughly 5–10 GB. Removing `hmm` and `smf_clean` afterwards gets that back.
- **Recommended order (inferred):** run `pixi lock` (solves versions without installing), then `pixi install` and the tests in a GPU job, then the Flax and Orbax code changes, then repoint the docs and sbatch scripts. A git worktree suits this because it leaves the conda workflow on the main checkout untouched.

# Details

## Pixi layout to copy from SMFNet_public

SMFNet_public's `pyproject.toml` has:

- `[tool.pixi.workspace]` with conda-forge and bioconda, platforms `linux-64` and `osx-arm64`;
- `[tool.pixi.system-requirements] libc = "2.17"`;
- `pypi-options.no-build = true`;
- features `gpu` (`system-requirements = { cuda = "12" }`), `cpu`, `dev`, `training` and `notebook`, each in its own solve group.

Its `.pixi/config.toml` sets:

```toml
detached-environments = "/home/groups/btrippe/diamant/pixi-envs"
```

`~/.bashrc` has a `pixi_smf()` helper that runs `module load devel pixi/0.53.0`, sets `PIXI_CACHE_DIR=$SCRATCH/pixi-cache`, and runs `eval "$(pixi shell-hook --frozen ...)"`.

Proposed for smf_hmm:

- **Base dependencies:** the JAX stack, numpy/pandas/scipy/h5py/pyfaidx/matplotlib/seaborn/tqdm/hydra-zen, and the bio libraries.
- **Features:**
  - `gpu`: jaxlib CUDA build;
  - `cpu`: CPU jaxlib;
  - `dev`: pytest, black;
  - `training`: wandb;
  - `notebook`: ipykernel, jupyterlab.
- **Environments:** `default` (GPU) and `cpu`.

The bio libraries go in the base deps rather than a separate feature because they are small and already sit in `hmm`.

Setup steps not to forget:

- `.pixi` is gitignored in smf_hmm. Run `pixi config set --local detached-environments /home/groups/btrippe/diamant/pixi-envs` in every checkout or worktree before the first install.
- On a node without a GPU driver, set `CONDA_OVERRIDE_CUDA=12` for `pixi install`.
- About 15 docs and sbatch files hard-code `.../miniforge/envs/hmm/bin/python` or `.../envs/smf_clean/bin/python`. Find them with `grep -rn 'envs/hmm\|envs/smf_clean\|conda activate' .`. The sbatch files already accept a `PYTHON` override. AGENTS.md and two tests still use Mac paths (`/Users/ndiamant/...`).

How the version and glibc data were gathered:

```bash
curl -s https://api.anaconda.org/package/conda-forge/jaxlib/files -o /tmp/jaxlib_files.json
# Then, per linux-64 py312 build: the __glibc and cuda-version entries of attrs.depends.
# Same for conda-forge/{flax,orbax-checkpoint,jax,optax,biopython} and bioconda/{pysam,pybigwig,pybigtools}.
nvidia-smi --query-gpu=name,driver_version --format=csv   # H100 80GB HBM3, 550.163.01
```

## Current checkpointing and its problems

`save_best_hmm_checkpoint` saves `nnx.to_pure_dict(nnx.state(model))` plus format version and meta with `StandardCheckpointer` to `checkpoints/best/`. It writes the model spec to a separate `checkpoints/best.model_spec.json`; this sidecar came from commit `9aa5767` (2026-02-25, "fix checkpointing bug due to orbax version"). An example on disk is `/scratch/users/diamant/2026-09-08_fixed_phase_type_kernel_9/checkpoints/best`, 1.2 MB. It stores the tree at the top level, with no `_CHECKPOINT_METADATA` and no `commit_success.txt`.

Problems:

1. **A crash mid-save loses the best checkpoint.** Both the old `best/` and the sidecar are deleted before the new ones are written.
2. **The spec lives outside the checkpoint**, so `best/` alone is incomplete.
3. **Restore doesn't tell Orbax what to expect.** `restore(path)` is called without an abstract tree, which Orbax's docs call the brittle path.
4. **A `Checkpointer` protocol and `_make_standard_checkpointer` exist only so tests can swap in a fake.**

## Orbax v1 API (read at v0.12.8 source; HEAD matched the tag on 2026-10-21)

- `ocp.v1.save(path, tree, overwrite=True, custom_metadata=<json>)` saves atomically: the tmp-then-commit mode is set by `AtomicityOptions`.
- `ocp.v1.metadata(path).custom_metadata` reads the JSON without loading arrays.
- `ocp.v1.load(path, abstract_state)` checks the structure and dtypes against `abstract_state`.
- `save_pytree` and `load_pytree` are deprecated aliases of `save` and `load`. The API is still changing, so pin the version.
- `ocp.v1.training.Checkpointer(directory, save_decision_policy=..., preservation_policy=...)` manages one checkpoint per step. It is only useful for resumable training; the optimizer state is not saved today. `BestN(get_metric_fn, n, reverse)` keeps the **highest** metric values unless `reverse=True`, so validation loss needs `reverse=True`.
- Old-format checkpoints: v1 has an `OrbaxV0Layout` that recognises a pytree stored at the top level, and `AtomicityOptions.allow_legacy_atomic_rename` for checkpoints with no `commit_success.txt`. **Untested on our checkpoints.**

Proposed rewrite (sketch, not applied):

```python
from orbax.checkpoint import v1 as ocp

def save_best_hmm_checkpoint(checkpoint_dir, model, *, epoch=None, step=None, validation_loss=None) -> None:
    ocp.save(
        _best_checkpoint_path(checkpoint_dir),
        nnx.state(model),
        overwrite=True,
        custom_metadata={
            "format_version": CHECKPOINT_FORMAT_VERSION,
            "model_spec": _to_plain_python(model_spec_from_hmm(model)),
            "meta": {...},
        },
    )

def load_best_hmm_checkpoint(checkpoint_dir, *, rng_seed=0):
    path = _best_checkpoint_path(checkpoint_dir)
    custom = ocp.metadata(path).custom_metadata
    model = hmm_from_model_spec(custom["model_spec"], rng_seed=rng_seed)
    abstract_state = jax.tree.map(lambda x: jax.ShapeDtypeStruct(x.shape, x.dtype), nnx.state(model))
    nnx.update(model, ocp.load(path, abstract_state))
    return model, custom
```

What the rewrite changes:

- **Removed:** the protocol, the fake-checkpointer hook, the sidecar helpers, the delete logic and `make_hmm_checkpoint_payload`. The model-spec registry stays.
- **Callers:** callers doing `model, _ = ...` are unaffected. `tests/test_checkpointing.py` and `tmp_plot_fixed_blocks_posteriors.py` use the second return value and would get the metadata dict instead of the old payload.
- **Old checkpoints:** a one-off conversion script is preferred over a permanent second load path, per the smf_hmm AGENTS.md style guide.

# Related Notes

- [Project index](../README.md): the experiments whose checkpoints a format change would affect.

# Open Questions

- Does jaxlib ≥ 0.7 (CUDA 12.9) run correctly and fast enough on driver 550.163? This needs a GPU-node test of `pytest` plus a short training run. Does any Sherlock GPU partition have a newer driver?
- Which existing `checkpoints/best` directories must stay loadable after the switch? That decides whether to write a conversion script.
- Do we want resumable training? Only that would justify `training.Checkpointer` and saving the optimizer state.
- Which other Flax 0.10 → 0.12 NNX changes break smf_hmm? This is unknown until the tests run in the new env.

# Sources

- `~/repos/SMFNet_public/pyproject.toml`, `README.md` (Installation), `.pixi/config.toml`
- `smf_hmm/checkpointing.py`, `smf_hmm/train.py`, `sherlock_env.yaml` in `~/repos/smf_hmm` (commit `1d86e57`)
- conda package metadata: `https://api.anaconda.org/package/conda-forge/jaxlib/files`, and the same endpoint for the other packages above
- Orbax source, tag v0.12.8: `checkpoint/orbax/checkpoint/experimental/v1/_src/{saving/saving.py, loading/loading.py, metadata/loading.py, layout/orbax_v0_layout.py, training/checkpointer.py}`, `checkpoint/orbax/checkpoint/_src/checkpoint_managers/preservation_policy.py`, `checkpoint/CHANGELOG.md` — https://github.com/google/orbax
- Flax NNX checkpointing guide: https://github.com/google/flax/blob/main/docs_nnx/guides/checkpointing.ipynb
