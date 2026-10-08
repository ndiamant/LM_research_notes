---
title: What ends observed SMF runs, and whether to bridge missing gaps in run-length metrics
date: 2026-10-08
project: smfnet
agent: Claude Code
status: draft
sources:
  - SMFNet repo, src/smf_net/evaluation/run_lengths.py and metrics.py (branch nd/more_assays, uncommitted)
  - /scratch/users/diamant/evaluate_SMFNet_test_set/old_full_model/{sample/samples.h5,metrics_v3,metrics_v4}
  - /scratch/users/diamant/evaluate_SMFNet_fiber/test_samples_ont_polarity/samples.h5
  - /scratch/groups/btrippe/ndiamant/gch_SMF/20260920_SMF_MM_{C2C12,ES,MEL,NP}_NO_combined_baitregions_2048bp_reassigned_compressed.h5
  - User discussion, 2026-10-08
tags: [smfnet, evaluation, run-lengths, wasserstein, missingness, ont, fiberseq, gch]
---

# Summary

SMFNet's run-length Wasserstein metric ends a run at any missing call. We checked
what actually ends observed runs on three data types, and tried **bridging**:
letting a run continue across a missing gap when the calls on both sides agree
and lie at most a capped distance apart (30 bp).

- **ONT GpC (mESC):** 59% of interior runs end at an internal missing gap,
  usually a single scattered missing site. Bridging changed the results a lot
  (see below).
- **Fiber-seq m6A (GM12878):** 97% of runs end at a real 0↔1 change. Bridging
  would touch about 3% of runs.
- **GCH short reads (C2C12, ES, MEL, NP):** about 60% of runs end at a **read
  end**. Internal gaps are rare and mostly long. Bridging would touch under 1%
  of runs.

**Decision: bridging was not adopted.** It helps only ONT. It costs extra
explanation, and with a default of 0 the code was not worth keeping. The
implementation was removed before commit and is kept in
[bridge_missing.py](run-length-missingness-gap-bridging/bridge_missing.py),
along with instructions for restoring it.

**The bigger unresolved issue is read-end truncation on GCH.** It is not
missingness inside reads. See Open Questions.

# Key Points

- **Missing calls are what break runs, and runs ending at a missing call are
  kept.** `run_length_bp` drops runs that touch the window edge, but keeps runs
  ending at a missing call or a read end. Generated reads get observed
  per-read missingness donated to them, so the comparison is fair. But on ONT,
  and on GCH because of read ends, much of the distribution reflects
  missingness rather than chromatin.
- **Dropping missing-bounded runs is a bad fix.** It throws away most runs
  (63% on ONT), raising the noise floor by about √(1/0.37) ≈ 1.6×, and long
  runs are the most likely to be cut. That is why we tried bridging instead.
- **The cap has to be small.** The chance that a bridged gap hides a call of
  the other value rises steadily with gap span. On ONT it is under 5% up to
  20 bp and about 12% at 30–40 bp (zeros).
- **On ONT, bridging changed the gap closed in opposite directions for the
  two run types:**
  - Zeros got worse, from 95.2% to 93.4%. Bridging exposed a real deficit: the
    model produces too few long zero runs (≥200 bp: 1.3% of runs vs 1.6%
    observed).
  - Ones improved, from 90.4% to 95.4%. A likely explanation, **not tested**:
    observed missing calls are not placed at random relative to the calls
    around them, so randomly donated missingness breaks up model runs
    differently from real ones.
- **The pseudo-replicate floor didn't move relative to scale.** Floor ÷ shuffle
  null was 0.185 → 0.186 for zeros and 0.159 → 0.157 for ones. Only the model
  moved, so the changes are not a floor artefact.

# Details

## Setup

Run-length definitions are in `smf_net.evaluation.run_lengths`.
`donate_missingness` masks each generated read with a random observed read's
missing pattern. `run_length_bp` measures the bp span of maximal runs and
drops runs touching the window edges.

The analysis script is
[analyze_run_boundaries.py](run-length-missingness-gap-bridging/analyze_run_boundaries.py).
It classifies each side of every interior observed run as one of:

- a call of the other value;
- a read end (missing from the read's last call to the window edge);
- an internal missing gap.

It also simulates the bridging error rate: hide the interior calls of fully
called stretches whose two end calls agree, and check whether any hidden call
had the other value.

## What ends observed runs

| | ONT GpC (mESC) | Fiber-seq m6A | GCH C2C12 | GCH ES | GCH MEL | GCH NP |
|---|---|---|---|---|---|---|
| Median site spacing | 12.5 bp | 1 bp | 14.5 bp | 15 bp | 15 bp | 15 bp |
| Median read span | ~970 bp | n/a | 171 bp | 174 bp | 164 bp | 167 bp |
| Complete runs | 37.5% | 96.7% | 33% | 32% | 33% | 34% |
| Ends at an internal gap only | **58.7%** | 3.2% | 7.7% | 8.5% | 6.9% | 7.1% |
| Ends at a read end | 3.8% | 0.2% | **59%** | **60%** | **61%** | **59%** |
| Single-site share of internal gaps | 76% | 92% | 37% | 31% | 41% | 46% |
| Median internal-gap span | 23 bp | 3 bp | 93 bp | 106 bp | 88 bp | 78 bp |
| Internal gaps ≤30 bp | 64% | ~99% | 11% | 8% | 13% | 16% |

ONT and Fiber-seq used the first 300 and 150 regions of the samples files (60
reads per region for Fiber-seq). GCH used 1,500 random regions per file, 80
reads per region.

GCH internal gaps are long (29–43% are over 120 bp). Many are probably the
unsequenced stretch between paired-end mates rather than low-confidence calls
(**inferred**, not checked against read pairing).

## How often bridging would be wrong, by gap span

Share of simulated gaps whose hidden calls include the other value:

| Gap span | ONT zeros | ONT ones | Fiber-seq zeros | Fiber-seq ones | GCH zeros (4 cells) | GCH ones (4 cells) |
|---|---|---|---|---|---|---|
| ≤10 bp | 0.8% | 0.3% | 38% (≤5 bp) | 4% (≤5 bp) | 6–11% | 5–6% |
| 10–20 bp | 4.6% | 2.4% | 67% | 19% | 6–10% | 6–7% |
| 20–30 bp | 7.5% | 4.6% | 81% | 28% | 7–12% | 7–9% |
| 30–40 bp | 12% | 6% | 92% | 36% | 9–13% | 8–10% |
| 60–80 bp | 24% | 18% | — | — | 16–26% | 23–27% |

On Fiber-seq, zero runs are short and broken up by scattered calls of the
other value, so bridging zeros is unsafe even across a single site. On GCH the
error is higher than on ONT at every span.

## ONT re-run with a 30 bp cap

`metrics_v3` has no bridging and `metrics_v4` bridges with a 30 bp cap. Both
are on `old_full_model`, 13,744 test regions. Bridging uses no random draws,
so every other column is identical between the two.

| | Model | Floor | Null | Gap closed |
|---|---|---|---|---|
| Zeros, v3 | 3.64 | 3.01 | 16.3 | 95.2% |
| Zeros, v4 | 5.13 | 3.98 | 21.4 | 93.4% |
| Ones, v3 | 7.59 | 5.04 | 31.6 | 90.4% |
| Ones, v4 | 7.95 | 6.36 | 40.6 | 95.4% |

- Runs per region (median) fell from 773 to 655 for zeros and from 831 to 680
  for ones.
- Mean observed run length rose from 23.9 to 30.0 bp for zeros and from 54.4
  to 66.5 bp for ones.
- The share of runs ≥100 bp rose from 4.7% to 8.1% for zeros and from 19.9% to
  28.5% for ones (observed).
- The difference between the pooled model and observed histograms (total
  variation) went from 0.009 to 0.018 for zeros and from 0.059 to 0.043 for
  ones.

A Fiber-seq metrics comparison with bridging on and off was started and then
stopped once the decision was made, so there are no Fiber-seq metric numbers.

## Reproduce

```bash
cd /home/users/diamant/repos/LM_research_notes/projects/smfnet/2026-10-08/run-length-missingness-gap-bridging
PY="$GROUP_HOME/uv-envs/SMFNet-cu124/bin/python"   # SMFNet must be importable
CUDA_VISIBLE_DEVICES="" $PY analyze_run_boundaries.py samples \
    /scratch/users/diamant/evaluate_SMFNet_test_set/old_full_model/sample/samples.h5 --regions 300
CUDA_VISIBLE_DEVICES="" $PY analyze_run_boundaries.py samples \
    /scratch/users/diamant/evaluate_SMFNet_fiber/test_samples_ont_polarity/samples.h5 \
    --regions 150 --reads-per-region 60
for CELL in C2C12 ES MEL NP; do
  CUDA_VISIBLE_DEVICES="" $PY analyze_run_boundaries.py gch \
    /scratch/groups/btrippe/ndiamant/gch_SMF/20260920_SMF_MM_${CELL}_NO_combined_baitregions_2048bp_reassigned_compressed.h5 \
    --regions 1500 --reads-per-region 80
done
```

The script draws its own random stretches, so the error-rate cells can differ
slightly from the tables above. The ONT metrics runs used `scripts/metrics.py`
from the SMFNet repo:

```bash
cd ~/repos/SMFNet
CUDA_VISIBLE_DEVICES="" MPLBACKEND=Agg $GROUP_HOME/uv-envs/SMFNet-cu124/bin/python scripts/metrics.py \
  samples_path=/scratch/users/diamant/evaluate_SMFNet_test_set/old_full_model/sample/samples.h5 \
  metrics.num_workers=8 plot.enabled=True \
  hydra.run.dir=/scratch/users/diamant/evaluate_SMFNet_test_set/old_full_model/metrics_v4
# metrics_v4 also had the bridging code with metrics.run_length_max_bridge_bp=30 (the default then)
```

# Related Notes

- [Train-set positional baselines and noise-matched correlation-of-correlations](../2026-09-29/positional-covariation-baselines-and-noise-matched-corr-of-corr.md):
  The same evaluation-noise thinking (matched read counts, pseudo-replicate
  floors) behind the run-length floor comparison here.
- [Fiber-seq GM12878 data preparation](../2026-10-05/fiberseq-gm12878-data-preparation.md):
  Source of the Fiber-seq reads whose run boundaries are summarized here.

# Open Questions

- **GCH read-end truncation.** About 60% of GCH runs are cut off by a ~170 bp
  read ending inside the window, so nucleosome-length zero runs almost never
  appear uncut. Dropping those runs would remove nearly every long run. A
  Kaplan–Meier-style estimate that treats them as censored is the principled
  fix, but it hasn't been tried.
- **Is the ones improvement really caused by where missing calls fall?** To
  test: compare the observed calls next to real missing sites with those next
  to random sites. If the placement isn't random, donated missingness is a
  biased stand-in, with or without bridging.
- **Are the long GCH internal gaps the stretch between paired-end mates?**
  Check against read pairing in the source BAMs.

# Sources

- SMFNet `src/smf_net/evaluation/run_lengths.py` (`donate_missingness`, `run_length_bp`).
- Analysis script: [analyze_run_boundaries.py](run-length-missingness-gap-bridging/analyze_run_boundaries.py).
- Removed implementation: [bridge_missing.py](run-length-missingness-gap-bridging/bridge_missing.py).
- Metrics outputs: `/scratch/users/diamant/evaluate_SMFNet_test_set/old_full_model/metrics_v3/` (no bridging) and `metrics_v4/` (30 bp cap).
- Samples: `/scratch/users/diamant/evaluate_SMFNet_test_set/old_full_model/sample/samples.h5` (ONT mESC, source `/scratch/users/diamant/smf_data/ONT/ONT_mESC_acRegions.h5`) and `/scratch/users/diamant/evaluate_SMFNet_fiber/test_samples_ont_polarity/samples.h5`.
- GCH data: `/scratch/groups/btrippe/ndiamant/gch_SMF/20260920_SMF_MM_{C2C12,ES,MEL,NP}_NO_combined_baitregions_2048bp_reassigned_compressed.h5`.
