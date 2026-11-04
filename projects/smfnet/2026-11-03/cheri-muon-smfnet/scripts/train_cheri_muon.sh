#!/bin/bash
# Cheri blocks (DNA encoder and SMF body) + Muon, otherwise the baseline
# config. Runs from the nd_cheri worktree with the main checkout's pixi env.
set -euo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
export PYTHONFAULTHANDLER=1
cd $HOME/repos/SMFNet_cheri
git log -1 --oneline
R=/scratch/users/diamant/smf_models/ONT_mESC_cheri_muon
PYTHONPATH=$PWD/src pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml \
    python scripts/train.py hydra.run.dir=$R \
    wandb_logger.save_dir=/scratch/users/diamant/wandb wandb_logger.name=ONT_mESC_cheri_muon \
    assay_type=cg_gc_accessibility \
    "data.cfg.h5_path=[/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions.h5,/scratch/groups/btrippe/ndiamant/ONT_mESC_bgRegions.h5]" \
    data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/mm10.fa \
    model.cfg.cnn_config.block_type=cheri model.cfg.cnn_config.cheri_checkpoint=true \
    model.cfg.optimizer=muon log_mask_rate_bins=true checkpoint_every_n_epochs=25
