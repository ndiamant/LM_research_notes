#!/bin/bash
# 1 bp attention U-Net + Muon (lr 0.01) + EMA 0.999 at 128 regions x 64
# reads, max_epochs=200 with 5 warmup epochs: the same steps, reads per step
# and LR schedule as ONT_mESC_alphagenome_unet_attn_muon_ema (256 x 32, 400
# epochs), with each region visited half as often. Does it memorize less DNA?
# Then: fixed-mask val loss (32 reads), memorization diagnostics, the 32-step
# eval and the 64-step subset. Run with: bash train_eval_attn_muon_ema_r64.sh
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
export CODE=$HOME/repos/SMFNet_patch
cd $CODE && git log -1 --oneline
export PYTHONPATH=$CODE/src
M=/scratch/users/diamant/smf_models
H=/scratch/groups/btrippe/ndiamant
S=$SCRATCH/smfnet_unet
G=$SCRATCH/smfnet_memorization
N=ONT_mESC_alphagenome_unet_attn_muon_ema_r64
R=$M/$N
COMMON=$(cat $S/common_args.txt)
echo "=== train $(date)"
eval $PIXI python scripts/train.py hydra.run.dir=$R wandb_logger.save_dir=/scratch/users/diamant/wandb \
    wandb_logger.name=$N checkpoint_every_n_epochs=25 \
    model.cfg.optimizer=muon model.cfg.muon_lr=0.01 ema_decay=0.999 \
    trainer.max_epochs=200 model.cfg.warmup_epochs=5 data.cfg.batch_size=128 data.cfg.reads_per_region=64 \
    $COMMON model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true \
    model.cfg.cnn_config.output_conv_matmul=true > $R.log 2>&1
echo "=== train exit $? $(date)"
BEST=$(ls $R/epoch=*-step=*.ckpt | head -1); BEST_FILE=$(basename "$BEST")
TAG=${BEST_FILE%%-step*}; TAG=${TAG/=/}
echo "best $BEST_FILE ($TAG)"
$PIXI python $S/rescore_val.py --name attn_muon_ema_r64 --checkpoint "$BEST" --out $S/rescored/attn_muon_ema_r64.json || echo "FAILED rescore best"
$PIXI python $S/rescore_val.py --name attn_muon_ema_r64_e199 --checkpoint "$R/periodic/epoch=199.ckpt" --out $S/rescored/attn_muon_ema_r64_e199.json || echo "FAILED rescore e199"
for kind in ac bg; do for split in train val; do
  $PIXI python $G/memorization_gap.py --checkpoint "$BEST" --h5 $H/ONT_mESC_${kind}Regions.h5 --split $split \
      --out $G/results/attn_muon_ema_r64_${kind}_${split}.json || echo "FAILED memorization $kind $split"
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
