#!/bin/bash
# Full run: k=9 attention U-Net + dilated second down convs at the three
# finest levels (dilations 4, 4, 2), Muon (lr 0.01) + EMA 0.999, 256 regions
# x 32 reads, 400 epochs. Compare with ONT_mESC_alphagenome_unet_attn_muon_ema_k9
# (fixed-mask val loss 0.29259). Then: fixed-mask val loss (best, epoch 399),
# span inpainting of epochs 99/199/299/399 (the drift over training),
# memorization diagnostics, the 32-step eval and the 64-step subset.
# Resumes from last.ckpt if rerun. Run with: bash train_eval_k9_dilated.sh
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
export CODE=$HOME/repos/SMFNet_inpaint_runs  # detached at ec73658
cd $CODE && git log -1 --oneline; hostname; nvidia-smi --query-gpu=name --format=csv,noheader
export PYTHONPATH=$CODE/src
M=/scratch/users/diamant/smf_models
H=/scratch/groups/btrippe/ndiamant
S=$SCRATCH/smfnet_unet
G=$SCRATCH/smfnet_memorization
I=$SCRATCH/smfnet_inpaint
N=ONT_mESC_alphagenome_unet_attn_muon_ema_k9_dilated
R=$M/$N
RESUME=""
if [ -f $R/last.ckpt ] && [ ! -f $R/periodic/epoch=399.ckpt ]; then
  RESUME="ckpt_path=$R/last.ckpt"; mv $R/csv $R/csv_before_resume_$(date +%s); echo "resuming"
fi
if [ ! -f $R/periodic/epoch=399.ckpt ]; then
  echo "=== train $(date)"
  COMMON=$(cat $S/common_args.txt)
  eval $PIXI python scripts/train.py hydra.run.dir=$R wandb_logger.save_dir=/scratch/users/diamant/wandb \
      wandb_logger.name=$N checkpoint_every_n_epochs=25 \
      model.cfg.optimizer=muon model.cfg.muon_lr=0.01 ema_decay=0.999 \
      model.cfg.cnn_config.unet_kernel_size=9 "model.cfg.cnn_config.unet_down_dilations=[4,4,2,1,1,1,1]" \
      $COMMON model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true \
      model.cfg.cnn_config.output_conv_matmul=true $RESUME >> $R.log 2>&1
  echo "=== train exit $? $(date)"
fi
BEST=$(ls $R/epoch=*-step=*.ckpt | head -1); BEST_FILE=$(basename "$BEST")
TAG=${BEST_FILE%%-step*}; TAG=${TAG/=/}
echo "best $BEST_FILE ($TAG)"
$PIXI python $S/rescore_val.py --name k9_dilated --checkpoint "$BEST" --out $S/rescored/k9_dilated.json || echo "FAILED rescore best"
$PIXI python $S/rescore_val.py --name k9_dilated_e399 --checkpoint "$R/periodic/epoch=399.ckpt" --out $S/rescored/k9_dilated_e399.json || echo "FAILED rescore e399"
for e in 99 199 299 399; do
  $PIXI python $I/span_infill.py --checkpoint "$R/periodic/epoch=$e.ckpt" --out $I/full_k9_dilated_e$e.json 2>&1 \
      | grep "^span" || echo "FAILED inpaint e$e"
done
for kind in ac bg; do for split in train val; do
  $PIXI python $G/memorization_gap.py --checkpoint "$BEST" --h5 $H/ONT_mESC_${kind}Regions.h5 --split $split \
      --out $G/results/k9_dilated_${kind}_${split}.json || echo "FAILED memorization $kind $split"
done; done
bash -x $S/eval_run_code.sh $N "$BEST_FILE" $TAG || echo "FAILED eval_run"
out=$R/samples_val_ac_subset1000_steps64
$PIXI python scripts/sample.py hydra.run.dir=$out "checkpoint='$BEST'" \
    regions.source=split regions.split=val sampling.sample_steps=64 \
    data.cfg.fasta_path=$H/mm10.fa "data.cfg.h5_path=[$H/ONT_mESC_acRegions.h5]" \
    data.cfg.max_regions_per_h5_per_split=1000 \
&& $PIXI python scripts/metrics.py hydra.run.dir=${out/samples_/metrics_} \
    samples_path=$out/samples.h5 metrics.include_positional_baselines=false metrics.include_flank_baselines=false \
|| echo "FAILED 64-step"
echo "=== all done $(date)"
