#!/usr/bin/env bash
# One overnight regularization run: ES gch_accessibility on a fixed 15,000-region
# subset per split, stopped after MAX_TIME of wall time. The learning-rate
# schedule keeps max_epochs=400, matching the production runs, so stopping
# early does not change how fast the learning rate decays.
#
# Usage: run_one.sh NAME [extra hydra overrides...]
# Run from the SMFNet repo root on a GPU node with the SMFNet env.
set -euo pipefail

NAME="$1"; shift
ROOT="${SCRATCH}/smf_models/overnight_20261006"
MAX_TIME="${MAX_TIME:-00:00:45:00}"
REGIONS="${REGIONS:-15000}"
TRAIN_REGIONS="${TRAIN_REGIONS:-5000}"
PYTHON="${PYTHON:-$GROUP_HOME/uv-envs/SMFNet-cu124/bin/python}"

mkdir -p "$ROOT/logs"
echo "[$(date '+%F %T')] start $NAME: $*" >> "$ROOT/logs/runs.txt"
WANDB_MODE=offline "$PYTHON" scripts/train.py \
  hydra.run.dir="$ROOT/$NAME" \
  wandb_logger.save_dir="$ROOT/wandb" wandb_logger.offline=true wandb_logger.name="overnight_$NAME" \
  assay_type=gch_accessibility \
  data.cfg.h5_path=/scratch/groups/btrippe/ndiamant/gch_SMF/20260920_SMF_MM_ES_NO_combined_baitregions_2048bp_reassigned_compressed.h5 \
  data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/mm10.fa \
  data.cfg.h5_region_length=2048 data.cfg.dna_sub_one_from_end_idx=False \
  data.cfg.rc_augmentation=false \
  data.cfg.max_regions_per_h5_per_split="$REGIONS" \
  data.cfg.max_train_regions_per_h5="$TRAIN_REGIONS" \
  "+trainer.max_time='$MAX_TIME'" log_layer_norms=true \
  "$@" > "$ROOT/logs/$NAME.log" 2>&1 && status=0 || status=$?
echo "[$(date '+%F %T')] end $NAME (exit $status)" >> "$ROOT/logs/runs.txt"
