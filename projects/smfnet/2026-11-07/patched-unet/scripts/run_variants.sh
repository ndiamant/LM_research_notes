#!/bin/bash
# Patched U-Net variants (nd_unet_patch), AdamW + EMA 0.999, first $EPOCHS
# epochs of the 400-epoch schedule, one after another on this GPU.
# Run with: bash run_variants.sh [variant ...]
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
CODE=$HOME/repos/SMFNet_patch
cd $CODE
git log -1 --oneline
export PYTHONPATH=$CODE/src
EPOCHS=${EPOCHS:-30}
COMMON=$(cat $SCRATCH/smfnet_unet/common_args.txt)
FAST="model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true model.cfg.cnn_config.output_conv_matmul=true"
declare -A ARGS=(
    [p4_narrow]="model.cfg.cnn_config.unet_patch_size=4 model.cfg.cnn_config.unet_channels=[80,96,112,128,144]"
    [p4_mid]="model.cfg.cnn_config.unet_patch_size=4 model.cfg.cnn_config.unet_channels=[96,128,160,192,224]"
    [p4_wide]="model.cfg.cnn_config.unet_patch_size=4 model.cfg.cnn_config.unet_channels=[128,160,192,224,256]"
    [p2]="model.cfg.cnn_config.unet_patch_size=2 model.cfg.cnn_config.unet_channels=[64,80,96,112,128,144]"
    [p4_wide_h2]="model.cfg.cnn_config.unet_patch_size=4 model.cfg.cnn_config.unet_channels=[128,160,192,224,256] model.cfg.cnn_config.unet_head_channels=32 model.cfg.cnn_config.unet_head_blocks=2"
)
for v in "$@"; do
    R=/scratch/users/diamant/smf_models/ONT_mESC_unet_attn_adamw_ema_$v
    echo "=== $v -> $R  $(date)"
    eval pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml python scripts/train.py \
        hydra.run.dir=$R wandb_logger.save_dir=/scratch/users/diamant/wandb \
        wandb_logger.name=ONT_mESC_unet_attn_adamw_ema_$v ema_decay=0.999 stop_after_epochs=$EPOCHS \
        $COMMON $FAST ${ARGS[$v]} > $R.log 2>&1
    echo "=== $v exit $?  $(date)"
done
