#!/bin/bash
#SBATCH --job-name=smfnet_baseline_val_resample
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_ag_smoke/baseline_val_resample_%j.out

# Resample the DNA-only baseline on the val active regions with the current
# sampling code (same settings as the AlphaGenome and ATAC runs), then compute
# metrics with plots, to check the old samples (schema v1) aren't the reason
# for the AlphaGenome gap.
set -euo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
cd $HOME/repos/SMFNet_public
B=/scratch/users/diamant/smf_models/SMFNet_full_run_2026-08-11
git log -1 --oneline
pixi run --frozen python scripts/sample.py hydra.run.dir=$B/samples_val_ac_epoch377_current_code \
    "checkpoint='$B/epoch=377-step=171990-converted.ckpt'" \
    regions.source=split regions.split=val \
    data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/mm10.fa \
    "data.cfg.h5_path=[/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions.h5]"
pixi run --frozen python scripts/metrics.py hydra.run.dir=$B/metrics_val_ac_epoch377_current_code_plots \
    samples_path=$B/samples_val_ac_epoch377_current_code/samples.h5 plot.enabled=true
