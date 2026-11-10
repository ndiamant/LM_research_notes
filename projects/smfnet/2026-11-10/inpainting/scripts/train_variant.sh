#!/bin/bash
# One inpainting variant of the k=9 attention U-Net (Muon lr 0.01 + EMA
# 0.999, 256 regions x 32 reads), trained for the first 50 epochs of the
# 400-epoch schedule, then: fixed-mask val loss and span inpainting of the
# epoch-49 checkpoint. Resumes from the run's last.ckpt when requeued.
# Code: the detached checkout ~/repos/SMFNet_inpaint_runs (commit ec73658).
# Run with: bash train_variant.sh <variant>
set -uo pipefail
V=$1
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
CODE=$HOME/repos/SMFNet_inpaint_runs
cd $CODE && git log -1 --oneline; hostname; nvidia-smi --query-gpu=name --format=csv,noheader
export PYTHONPATH=$CODE/src
M=/scratch/users/diamant/smf_models
S=$SCRATCH/smfnet_inpaint
N=ONT_mESC_unet_attn_muon_ema_k9_$V
R=$M/$N
FUSED="model.cfg.cnn_config.unet_fused_norm=true"
case $V in
  control) ARGS="" ;;
  linear) ARGS="model.cfg.cnn_config.unet_upsample=linear" ;;
  dilated) ARGS="model.cfg.cnn_config.unet_down_dilations=[4,4,2,1,1,1,1]" ;;
  local_attn) ARGS="model.cfg.cnn_config.unet_local_attention_levels=[2]" ;;
  bottleneck16) ARGS="model.cfg.cnn_config.unet_channels=[48,64,80,96,112]" ;;
  read_masked) ARGS="model.cfg.cnn_config.unet_norm=read_masked"; FUSED="model.cfg.cnn_config.unet_fused_norm=false" ;;
  window) ARGS="model.cfg.cnn_config.unet_norm=window"; FUSED="model.cfg.cnn_config.unet_fused_norm=false" ;;
  *) echo "unknown variant $V"; exit 1 ;;
esac
RESUME=""
if [ -f $R/last.ckpt ] && [ ! -f $R/periodic/epoch=49.ckpt ]; then
  RESUME="ckpt_path=$R/last.ckpt"
  mv $R/csv $R/csv_before_resume_$(date +%s)  # the CSV logger would overwrite it
  echo "resuming from $R/last.ckpt"
fi
if [ ! -f $R/periodic/epoch=49.ckpt ]; then
  echo "=== train $(date)"
  COMMON=$(cat $SCRATCH/smfnet_unet/common_args.txt)
  eval $PIXI python scripts/train.py hydra.run.dir=$R wandb_logger.save_dir=/scratch/users/diamant/wandb \
      wandb_logger.name=$N checkpoint_every_n_epochs=25 stop_after_epochs=50 \
      model.cfg.optimizer=muon model.cfg.muon_lr=0.01 ema_decay=0.999 \
      model.cfg.cnn_config.unet_kernel_size=9 $COMMON $FUSED \
      model.cfg.cnn_config.unet_channels_last=true model.cfg.cnn_config.output_conv_matmul=true \
      $ARGS $RESUME >> $R.log 2>&1
  echo "=== train exit $? $(date)"
fi
CKPT=$R/periodic/epoch=49.ckpt
[ -f $CKPT ] || { echo "no epoch-49 checkpoint"; exit 1; }
[ -f $S/rescored_e49_$V.json ] || $PIXI python $SCRATCH/smfnet_unet/rescore_val.py --name ${V}_e49 \
    --checkpoint $CKPT --out $S/rescored_e49_$V.json || echo "FAILED rescore"
[ -f $S/pilot_e49_$V.json ] || $PIXI python $S/span_infill.py --checkpoint $CKPT \
    --out $S/pilot_e49_$V.json 2>&1 | grep -v "Warning\|it/s\]\|^batch\|weights_only" || echo "FAILED inpaint"
echo "=== all done $(date)"
