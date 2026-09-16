---
title: Multi-TF motif caches and state-call visualization plan
date: 2026-09-16
project: smf-hmm
agent: Codex
status: complete
sources:
  - /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/configs/mesc_validation.yaml
  - /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/benchmark_cache.py
  - /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/select_examples.py
  - /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/plot_examples.py
  - /scratch/users/diamant/mesc_validation/benchmark_cache/mesc_tf_chip_validation_v1
tags:
  - smf-hmm
  - visualization
  - motif-annotation
  - chip-nexus
  - benchmark-cache
  - slurm
  - mesc
---

# Summary

The chr8 TF/ChIP benchmark caches for ESRRB, KLF4, and OCT4 were built to
support a qualitative state-call visualization that shows the anchor motif
and other nearby high-confidence ChIP-supported motifs, including motifs for
other TFs. Three independent Slurm jobs were submitted to the
`hns,stat,akundaje` partition list on 2026-09-16. All completed successfully,
their stderr logs were empty, and all three published cache shards passed
`validate_cache_shard`.

The intended visualization design keeps secondary motif annotations separate
from the frozen example/read-selection manifest. A derived annotation table
will join the selected plot windows against named TF cache shards, and the
plotter will render a dedicated motif track below the hard state calls. This
prevents additional annotations from changing the selected examples or reads.

# Key Points

- ESRRB job `43840467` completed on `akundaje` in 1 minute 59 seconds.
- KLF4 job `43840469` completed on `akundaje` in 2 minutes 27 seconds.
- OCT4 job `43840470` completed on `stat` in 1 minute 29 seconds.
- Each job requested one CPU, 32 GB RAM, and four hours, and was eligible for
  `hns`, `stat`, or `akundaje`.
- Validated chr8 caches contain 46,709 ESRRB sites, 106,944 KLF4 sites, and
  41,436 OCT4 sites, with six ChIP–nexus replicates per TF.
- The existing CTCF chr8 cache remains the anchor source and contains seven
  replicates. Under the current strong-site filters, it contains 497 strong
  CTCF sites; two of the five selected CTCF panels contain one additional
  strong CTCF site inside the displayed ±250 bp window.
- Candidate secondary motifs should initially use TF- and chromosome-relative
  confidence filters rather than comparing raw ChIP sums between TFs.

# Details

## Goal and proposed representation

The goal is to make hard-state plots interpretable in their local regulatory
context. Each page should retain its designated anchor motif while also
showing nearby canonical motifs that have strong matching ChIP signal. Sites
may come from CTCF, ESRRB, KLF4, or OCT4.

The minimal implementation is:

1. Add an `annotate_examples.py` step that accepts the frozen example manifest
   and repeated named cache shards such as `CTCF=.../ctcf/chr8` and
   `ESRRB=.../esrrb/chr8`.
2. Write one annotation row per motif-window overlap, including
   `example_rank`, TF, site ID, interval, strand, motif score and percentile,
   pooled ChIP sum and percentile, replicate support, and `is_anchor`.
3. Add a `--motif-annotations` input to `plot_examples.py`.
4. Render TF-colored motif boxes in a small track below the bottom state-call
   row. Use a heavier outline for the anchor, stack colliding boxes into lanes,
   and put TF names in the legend rather than repeatedly labeling boxes.

The annotations should not be added directly to the long-form read manifest:
that table repeats site information once per read, whereas nearby motifs are
page-level metadata and may change as confidence rules are refined.

For the first pass, reuse the qualitative-example confidence criteria:

- pooled ChIP signal at or above the 99th percentile for that TF/chromosome;
- motif score at or above the 90th percentile; and
- at least half the replicates individually at or above their 95th signal
  percentile.

These percentile and replicate-support criteria are preferable to a shared
raw-count threshold because TF datasets differ in sequencing depth and number
of replicates. The thresholds remain provisional and should be relaxed only
if the resulting tracks are too sparse.

## Preflight checks

Before submission, all configured MEME files and positive/negative BigWig
pairs were confirmed present. The target shard directories did not exist, the
shared genome background cache was present, and the three requested Slurm
partitions were up.

## Exact Slurm submission

The following was run from
`/home/users/diamant/repos/smf_hmm`:

```bash
set -euo pipefail
cache_root=/scratch/users/diamant/mesc_validation/benchmark_cache
benchmark_root="$cache_root/mesc_tf_chip_validation_v1"
for tf in esrrb klf4 oct4; do
  shard="$benchmark_root/$tf/chr8"
  if [[ -e "$shard" ]]; then
    echo "Existing target prevents submission: $shard" >&2
    exit 1
  fi
done
[[ -s "$benchmark_root/genome_background.json" ]]
mkdir -p /scratch/users/diamant/mesc_validation/logs/cache
for tf in esrrb klf4 oct4; do
  sbatch \
    --job-name="cache_${tf}_chr8" \
    --partition=hns,stat,akundaje \
    --requeue \
    --time=04:00:00 \
    --cpus-per-task=1 \
    --mem=32G \
    --chdir=/home/users/diamant/repos/smf_hmm \
    --output="/scratch/users/diamant/mesc_validation/logs/cache/${tf}-chr8-%j.out" \
    --error="/scratch/users/diamant/mesc_validation/logs/cache/${tf}-chr8-%j.err" \
    --wrap="MPLCONFIGDIR=/home/users/diamant/repos/smf_hmm/.mplconfig /home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python -m evaluation.tf_chip.benchmark_cache --config /home/users/diamant/repos/smf_hmm/evaluation/tf_chip/configs/mesc_validation.yaml --cache-root /scratch/users/diamant/mesc_validation/benchmark_cache --tf $tf --chrom chr8"
done
```

Slurm returned:

```text
Submitted batch job 43840467
Submitted batch job 43840469
Submitted batch job 43840470
```

The jobs can be monitored with:

```bash
squeue -j 43840467,43840469,43840470 \
  -o '%.18i %.24j %.18P %.10T %.10M %.4D %R'
```

Logs are under:

```text
/scratch/users/diamant/mesc_validation/logs/cache/
```

## Completion and cache validation

Final accounting was checked with:

```bash
for job in 43840467 43840469 43840470; do
  sacct -j "$job" \
    --format=JobIDRaw,JobName%24,Partition,State,Elapsed,ExitCode \
    -n -P | head -1
done
```

This reported:

| Job | TF | Partition | State | Elapsed | Exit code |
| ---: | --- | --- | --- | ---: | ---: |
| 43840467 | ESRRB | akundaje | completed | 00:01:59 | 0:0 |
| 43840469 | KLF4 | akundaje | completed | 00:02:27 | 0:0 |
| 43840470 | OCT4 | stat | completed | 00:01:29 | 0:0 |

The cache contract was then checked directly:

```bash
/home/groups/btrippe/diamant/miniforge/envs/smf_clean/bin/python - <<'PY'
from pathlib import Path
from evaluation.tf_chip.benchmark_cache import validate_cache_shard

root = Path(
    "/scratch/users/diamant/mesc_validation/benchmark_cache/"
    "mesc_tf_chip_validation_v1"
)
for tf in ("esrrb", "klf4", "oct4"):
    sites, metadata = validate_cache_shard(root / tf / "chr8")
    print(
        f"{tf}\t{len(sites)} sites\t"
        f"{len(metadata['replicates'])} replicates\t"
        f"{metadata['table']['sha256']}"
    )
PY
```

Validated tables:

| TF | Sites | Replicates | Table SHA-256 |
| --- | ---: | ---: | --- |
| ESRRB | 46,709 | 6 | `60873d83c67242a6b0ef0e162d37e15e9c94d141f191794d891cf08f9cbdf027` |
| KLF4 | 106,944 | 6 | `9432f88659f034ddfc3ef363e215b383d468d253d845b95ee45c72a2c9a5661f` |
| OCT4 | 41,436 | 6 | `50af395087b3b3bb1b64f7bf6af228fc75dcd7132030d948cd7109ca11d96e89` |

# Related Notes

- [CTCF motif and local-control ChIP–nexus signal collection](../2026-09-02/ctcf-motif-chip-nexus-local-controls.md): Establishes the motif scanning, strand-paired ChIP signal, and mm10 blacklist conventions reused by the shared caches.
- [Fixed-HMM posterior-sampling baseline comparison](../2026-09-04/fixed-hmm-posterior-sampling-comparison.md): Documents the state callers whose hard calls are being compared in the qualitative visualization.
- [Simple interpolation baseline against chr1 CTCF ChIP–nexus](../2026-09-04/simple-interpolation-ctcf-chip-benchmark.md): Records the first state-caller/ChIP benchmark and its motif-overlap interpretation.

# Open Questions

- Are the initial 99th-percentile ChIP and 90th-percentile motif thresholds too
  sparse for ESRRB, KLF4, or OCT4 within 501 bp panels?
- Should motif boxes encode strand with a chevron, or is strand sufficiently
  represented in annotation metadata and hover-free PDF legends?
- Should overlapping motif boxes be stacked globally by interval or in fixed
  TF-specific lanes?
- Should the same annotation-table design be extended to chr16 immediately,
  or first be validated visually on the existing five chr8 CTCF examples?

# Sources

- `/home/users/diamant/repos/smf_hmm/evaluation/tf_chip/configs/mesc_validation.yaml`
- `/home/users/diamant/repos/smf_hmm/evaluation/tf_chip/benchmark_cache.py`
- `/home/users/diamant/repos/smf_hmm/evaluation/tf_chip/CACHE_SPEC.md`
- `/home/users/diamant/repos/smf_hmm/evaluation/tf_chip/select_examples.py`
- `/home/users/diamant/repos/smf_hmm/evaluation/tf_chip/plot_examples.py`
- `/scratch/users/diamant/mesc_validation/benchmark_cache/mesc_tf_chip_validation_v1/{ctcf,esrrb,klf4,oct4}/chr8/`
- `/scratch/users/diamant/mesc_validation/logs/cache/`
