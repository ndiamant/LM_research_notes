#!/bin/bash
# p4_wide_h2 at 128 regions x 64 reads (same reads per step as 256 x 32),
# max_epochs=200 with 5 warmup epochs, so 15 epochs match the 256 x 32 run's
# first 30 epochs in steps, reads and learning rate. Then fixed-mask val loss
# (32 reads per region) of both. Run with: bash run_r64.sh
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
CODE=$HOME/repos/SMFNet_patch
cd $CODE && git log -1 --oneline
export PYTHONPATH=$CODE/src
M=/scratch/users/diamant/smf_models
COMMON=$(cat $SCRATCH/smfnet_unet/common_args.txt)
V=p4_wide_h2_r64
R=$M/ONT_mESC_unet_attn_adamw_ema_$V
eval $PIXI python scripts/train.py hydra.run.dir=$R wandb_logger.save_dir=/scratch/users/diamant/wandb \
    wandb_logger.name=ONT_mESC_unet_attn_adamw_ema_$V ema_decay=0.999 stop_after_epochs=15 \
    trainer.max_epochs=200 model.cfg.warmup_epochs=5 data.cfg.batch_size=128 data.cfg.reads_per_region=64 \
    $COMMON model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true \
    model.cfg.cnn_config.output_conv_matmul=true model.cfg.cnn_config.unet_patch_size=4 \
    "model.cfg.cnn_config.unet_channels=[128,160,192,224,256]" model.cfg.cnn_config.unet_head_channels=32 \
    model.cfg.cnn_config.unet_head_blocks=2 > $R.log 2>&1
echo "train exit $?"
S=$SCRATCH/smfnet_unet
for c in "p4_wide_h2_e29:$M/ONT_mESC_unet_attn_adamw_ema_p4_wide_h2/epoch=29-step=13650.ckpt" "p4_wide_h2_r64_e14:$(ls $R/epoch=14-*.ckpt)"; do
  $PIXI python $S/rescore_val.py --name ${c%%:*} --checkpoint "${c#*:}" --out $S/rescored/${c%%:*}.json || echo "FAILED ${c%%:*}"
done
echo all done
