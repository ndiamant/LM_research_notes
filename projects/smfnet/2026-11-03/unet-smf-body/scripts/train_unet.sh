#!/bin/bash
# Frozen AlphaGenome DNA encoder + U-Net SMF body (channels 48..144), otherwise
# the same as ONT_mESC_alphagenome_frozen. Runs from the nd_unet worktree with
# the main checkout's pixi env.
set -euo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
export PYTHONFAULTHANDLER=1
cd $HOME/repos/SMFNet_unet
git log -1 --oneline
R=/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_frozen_unet
PYTHONPATH=$PWD/src pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml \
    python scripts/train.py hydra.run.dir=$R \
    wandb_logger.save_dir=/scratch/users/diamant/wandb wandb_logger.name=ONT_mESC_alphagenome_frozen_unet \
    assay_type=cg_gc_accessibility \
    "data.cfg.h5_path=[/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions.h5,/scratch/groups/btrippe/ndiamant/ONT_mESC_bgRegions.h5]" \
    data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/mm10.fa \
    model.cfg.cnn_config.dna_encoder_type=alphagenome model.cfg.cnn_config.smf_body_type=unet \
    alphagenome_weights=/scratch/users/diamant/alphagenome_pytorch/model_all_folds.safetensors \
    log_mask_rate_bins=true checkpoint_every_n_epochs=25
