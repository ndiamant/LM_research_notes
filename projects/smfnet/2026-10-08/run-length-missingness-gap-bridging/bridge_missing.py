"""Removed SMFNet implementation of capped gap bridging for run lengths.

This lived in smf_net/evaluation/run_lengths.py on branch nd/more_assays
(uncommitted, 2026-10-08) and was removed before commit. To restore it:

1. add this function to run_lengths.py;
2. give run_length_metrics a ``max_bridge_bp: float = 0`` argument and, before
   computing any runs, pass every read set (predicted, observed, the three
   matched samples and each RunLengthShuffles field) through
   ``bridge_missing(reads, positions, max_bridge_bp)``;
3. add ``run_length_max_bridge_bp: float`` to MetricsConfig and forward it from
   compute_region_metrics.

Bridging consumes no random draws, so turning it on changes only the
run-length columns of a metrics run.
"""

import torch


def bridge_missing(
    reads: torch.Tensor,
    positions: torch.Tensor,
    max_bridge_bp: float,
) -> torch.Tensor:
    """Fill short missing gaps whose flanking calls agree with that value.

    A gap is filled when the calls on either side match and lie at most
    ``max_bridge_bp`` apart. Gaps with disagreeing flanks and gaps at a read's
    ends stay missing, and a non-positive ``max_bridge_bp`` disables bridging.
    """

    if max_bridge_bp <= 0 or reads.numel() == 0:
        return reads
    num_positions = reads.shape[1]
    valid = reads != -1
    index = torch.arange(num_positions).expand_as(reads)
    previous = torch.where(valid, index, -1).cummax(dim=1).values
    following = (
        torch.where(valid, index, num_positions).flip(1).cummin(dim=1).values.flip(1)
    )
    inside = (previous >= 0) & (following < num_positions)
    previous = previous.clamp(0, num_positions - 1)
    following = following.clamp(0, num_positions - 1)
    previous_call = reads.gather(1, previous)
    fill = (
        ~valid
        & inside
        & (previous_call == reads.gather(1, following))
        & (positions[following] - positions[previous] <= max_bridge_bp)
    )
    return torch.where(fill, previous_call, reads)
