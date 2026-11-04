#!/bin/bash
#SBATCH --job-name=smfnet_cheri_lr0.01
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/users/diamant/smf_models/ONT_mESC_cheri_muon_lr0.01_%j.out

# Same as ONT_mESC_cheri_muon, with Muon lr 0.01 instead of 0.025: the first
# run's training loss stalled while Muon's lr warmed up to 0.025.
set -euo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
export PYTHONFAULTHANDLER=1
cd $HOME/repos/SMFNet_cheri
git log -1 --oneline
R=/scratch/users/diamant/smf_models/ONT_mESC_cheri_muon_lr0.01
PYTHONPATH=$PWD/src pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml \
    python scripts/train.py hydra.run.dir=$R \
    wandb_logger.save_dir=/scratch/users/diamant/wandb wandb_logger.name=ONT_mESC_cheri_muon_lr0.01 \
    assay_type=cg_gc_accessibility \
    "data.cfg.h5_path=[/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions.h5,/scratch/groups/btrippe/ndiamant/ONT_mESC_bgRegions.h5]" \
    data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/mm10.fa \
    model.cfg.cnn_config.block_type=cheri model.cfg.cnn_config.cheri_checkpoint=true \
    model.cfg.optimizer=muon model.cfg.muon_lr=0.01 \
    log_mask_rate_bins=true checkpoint_every_n_epochs=25
