#!/bin/bash
# Full 400-epoch run: patched attention U-Net (4 bp patches, channels
# 128-256, two 32-channel 1 bp head blocks) + Muon (lr 0.01) + EMA 0.999,
# 256 regions x 32 reads, to compare with ONT_mESC_alphagenome_unet_attn_muon_ema.
# Run with: bash train_p4_wide_h2_muon_ema.sh [ckpt_path=<run dir>/last.ckpt]
set -euo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
CODE=$HOME/repos/SMFNet_patch
cd $CODE
git log -1 --oneline
export PYTHONPATH=$CODE/src
COMMON=$(cat $SCRATCH/smfnet_unet/common_args.txt)
R=/scratch/users/diamant/smf_models/ONT_mESC_unet_attn_p4_wide_h2_muon_ema
eval pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml python scripts/train.py \
    hydra.run.dir=$R wandb_logger.save_dir=/scratch/users/diamant/wandb \
    wandb_logger.name=ONT_mESC_unet_attn_p4_wide_h2_muon_ema checkpoint_every_n_epochs=25 \
    model.cfg.optimizer=muon model.cfg.muon_lr=0.01 ema_decay=0.999 $COMMON \
    model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true \
    model.cfg.cnn_config.output_conv_matmul=true model.cfg.cnn_config.unet_patch_size=4 \
    "model.cfg.cnn_config.unet_channels=[128,160,192,224,256]" model.cfg.cnn_config.unet_head_channels=32 \
    model.cfg.cnn_config.unet_head_blocks=2 "$@"
