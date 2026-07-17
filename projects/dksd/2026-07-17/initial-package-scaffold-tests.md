---
title: Initial reusable package scaffold and tests
date: 2026-07-17
project: dksd
agent: Codex
status: draft
sources:
  - /home/users/diamant/repos/DKSD/AGENTS.md
  - /home/users/diamant/repos/DKSD/pyproject.toml
  - /home/users/diamant/repos/DKSD/src/dksd/kernels.py
  - /home/users/diamant/repos/DKSD/src/dksd/losses.py
  - /home/users/diamant/repos/DKSD/tests/test_kernels.py
  - /home/users/diamant/repos/DKSD/tests/test_losses.py
tags:
  - dksd
  - phase-0
  - kernels
  - losses
  - testing
---

# Summary

We created the first reusable Python package scaffold for the DKSD experiments and verified the basic kernel and loss utilities with a small pytest suite. This is not yet the full Phase 0 mathematical validation; it is the supporting code layer needed before implementing the fixed-target Monte Carlo identity checks.

# Key Points

- Added a `src/dksd` package with the required public loss API: DSM loss, DKSD U-statistic, DKSD V-statistic, and hybrid DSM-DKSD loss.
- Added initial kernel classes with a common `gram(x, y=None, *, t=None)` API: ambient RBF, multiscale RBF, anchored energy, normalized anchored energy, direct angular, radial, and angular-times-radial product kernels.
- Added tests for kernel Gram shape, symmetry, and PSD behavior on random batches.
- Added tests for the exact V-statistic decomposition when the kernel diagonal is one.
- Verified the test suite with `/home/groups/btrippe/diamant/miniforge/envs/dm/bin/python -m pytest`: 10 tests passed.

# Details

The scaffold follows the project guidance in `AGENTS.md`: keep the initial implementation small, readable, and focused on the algebra needed for Phase 0 before starting neural training experiments.

The current loss implementation computes the DKSD U-statistic using the matrix identity

```text
(tr(Delta.T K Delta) - sum_i K_ii ||delta_i||^2) / (B(B - 1))
```

and avoids constructing an explicit `B x B x d` residual tensor. The V-statistic is computed as `tr(Delta.T K Delta) / B^2`. The tested decomposition is

```text
L_V = L_DSM / B + ((B - 1) / B) L_U
```

for unit-diagonal RBF Gram matrices.

The kernel tests currently check numerical PSD with `torch.linalg.eigvalsh` on fixed random batches. This establishes a basic implementation sanity check for the kernels, but it is not a substitute for the planned Phase 0 tests comparing Stein, score-error, and denoising representations.

# Related Notes

- None yet.

# Open Questions

- Which exact fixed target should be used first for Phase 0: a single Gaussian or the Phase 1-style eight-component Gaussian mixture?
- Should the Phase 0 Stein-vs-score-error comparison use analytic RBF derivatives or PyTorch autograd for the first implementation?
- What Monte Carlo batch sizes and confidence interval rule should be standardized for the `1e-3` identity checks?

# Sources

- `/home/users/diamant/repos/DKSD/AGENTS.md`: Defines the DKSD mission, Phase 0 acceptance criteria, required APIs, style guidance, and Python environment.
- `/home/users/diamant/repos/DKSD/pyproject.toml`: Defines the package metadata, `src` layout, pytest config, and development tooling.
- `/home/users/diamant/repos/DKSD/src/dksd/kernels.py`: Implements the initial kernel Gram APIs.
- `/home/users/diamant/repos/DKSD/src/dksd/losses.py`: Implements DSM, U-statistic, V-statistic, and hybrid losses.
- `/home/users/diamant/repos/DKSD/tests/test_kernels.py`: Tests Gram shape, symmetry, and PSD checks.
- `/home/users/diamant/repos/DKSD/tests/test_losses.py`: Tests the V-statistic decomposition and hybrid loss bookkeeping.

