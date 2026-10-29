#!/bin/bash
#SBATCH --job-name=smfnet_baseline_val_sample
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_ag_smoke/sample_baseline_val_%j.out

# Baseline (DNA-only) val active-region samples, matching the ATAC and
# AlphaGenome runs' samples_val_ac_epoch377 settings.
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
cd $HOME/repos/SMFNet_public
pixi run --frozen python scripts/sample.py hydra.run.dir=/scratch/users/diamant/smf_models/SMFNet_full_run_2026-08-11/samples_val_ac_epoch377 \
    "checkpoint='/scratch/users/diamant/smf_models/SMFNet_full_run_2026-08-11/epoch=377-step=171990-converted.ckpt'" \
    regions.source=split regions.split=val \
    data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/mm10.fa \
    "data.cfg.h5_path=[/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions.h5]"
