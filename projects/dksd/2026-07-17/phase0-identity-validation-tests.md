---
title: Phase 0 identity validation tests
date: 2026-07-17
project: dksd
agent: Codex
status: draft
sources:
  - /home/users/diamant/repos/DKSD/AGENTS.md
  - /home/users/diamant/repos/DKSD/src/dksd/phase0.py
  - /home/users/diamant/repos/DKSD/src/dksd/exact_scores.py
  - /home/users/diamant/repos/DKSD/tests/test_phase0_identities.py
  - /home/users/diamant/repos/DKSD/tests/test_exact_scores.py
  - /home/users/diamant/repos/DKSD
tags:
  - dksd
  - phase-0
  - identities
  - stein
  - denoising
---

# Summary

The DKSD repository now has deterministic pytest coverage for the full Phase 0 identity checks on a tractable Gaussian-mixture target. The implementation was committed in the DKSD repo as `38d5dae` with all 19 tests passing.

# Key Points

- Added validation-only utilities for an analytic imperfect score of the form `s_theta(x) = s_t(x) + A x + b`.
- Implemented the paired RBF Stein kernel formula used only for Phase 0 validation.
- Added tests comparing the Stein representation, score-error representation, denoising representation, and U-statistic estimator.
- Extended exact GMM score utilities with posterior component weights and posterior clean means, enabling a direct denoising identity check.
- Marked the Phase 0 identity-test checklist item complete in the DKSD `AGENTS.md` progress section.

# Details

The Phase 0 target is an eight-component Gaussian mixture with fixed diffusion noise. The exact marginal score is available in closed form from GMM component responsibilities. The validation tests use a fixed linear perturbation `A x + b` to create a nonzero score error without involving a neural model.

The implemented tests cover:

- Denoising identity: posterior clean mean implies the same score as the closed-form GMM marginal score.
- Stein identity: the RBF Stein kernel expectation matches the derivative-free score-error DKSD expectation within a Monte Carlo confidence rule.
- Denoising DKSD identity: the observable corruption-score residual representation matches the exact score-error representation within a Monte Carlo confidence rule.
- U-statistic unbiasedness: minibatch DKSD U-statistics agree with a high-precision population Monte Carlo estimate within combined standard error.

This closes the main gap identified in the scaffold note: the project now has actual Phase 0 algebraic validation tests, not just kernel/loss unit tests. The remaining open implementation work moves to the next project-plan item: a baseline score model.

# Related Notes

- [Initial reusable package scaffold and tests](initial-package-scaffold-tests.md): Records the package and kernel/loss scaffold that these Phase 0 tests build on.

# Open Questions

- The current Monte Carlo identity tests use fixed sample sizes and a five-standard-error acceptance rule. This is deterministic under the fixed seeds, but a later validation report may want to record numerical estimates and confidence intervals explicitly.
- The Stein-vs-score-error test currently uses analytic RBF derivatives. This is preferable for clarity, but an autograd cross-check could be added if derivative formulas are extended to other kernels.

# Sources

- `/home/users/diamant/repos/DKSD/AGENTS.md`: Defines Phase 0 tests and records the updated implementation status.
- `/home/users/diamant/repos/DKSD/src/dksd/phase0.py`: Contains validation-only score perturbation, RBF Stein kernel, score-error integrand, corruption score, and Monte Carlo helper code.
- `/home/users/diamant/repos/DKSD/src/dksd/exact_scores.py`: Provides exact GMM score and posterior clean mean calculations used in the denoising identity check.
- `/home/users/diamant/repos/DKSD/tests/test_phase0_identities.py`: Contains the deterministic Phase 0 identity tests.
- `/home/users/diamant/repos/DKSD/tests/test_exact_scores.py`: Checks the exact GMM score against finite differences.
- DKSD commit `38d5dae`: Local project commit containing the Phase 0 implementation and tests.
