#!/bin/bash
# GM12878 Fiber-seq through the figwig loader, on the H5's settings: the k=9 attention
# U-Net with dilated fine down convs (4, 4, 2), Muon (lr 0.01) + EMA 0.999, 256 regions x
# 32 reads, 400 epochs, as smfnet_dilated/train_eval_k9_dilated.sh. Regions are the H5's
# kept windows (kept_windows.bed), 2114 bp with 2048 bp crops, reads overlapping by 500
# bp, 32 molecules spanning each region, ChromBPNet fold 0; the old architecture's
# GM12878_fiberseq trained on the same windows. Then a 1,000-region val sample + metrics
# (no positional baselines yet). Resumes from last.ckpt if rerun.
# Run with: sbatch train_k9_dilated_fiber_2kb.sbatch
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml"
export CODE=$HOME/repos/SMFNet_unet_figwig_runs  # detached at 6ce74dd
cd $CODE && git log -1 --oneline; hostname; nvidia-smi --query-gpu=name --format=csv,noheader
export PYTHONPATH=$CODE/src
M=/scratch/users/diamant/smf_models
F=$SCRATCH/fiberseq/GM12878
FASTA=/scratch/groups/btrippe/brian/projects/SMF/fiberseq/reference/hg38.fa
N=GM12878_fiberseq_figwig_alphagenome_unet_attn_muon_ema_k9_dilated
R=$M/$N
DATA="data=figwig data.cfg.fasta_path=$FASTA \
  data.cfg.figwig_dir=$F/figwig data.cfg.regions_bed=$F/kept_windows.bed \
  data.cfg.region_length=2114 data.cfg.model_input_length=2048 \
  data.cfg.min_overlap=500 data.cfg.min_spanning_molecules=32 \
  \"data.cfg.train_chrm=[chr2,chr4,chr5,chr7,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr21,chr22,chrX]\" \
  \"data.cfg.val_chrm=[chr8,chr20]\" \"data.cfg.test_chrm=[chr1,chr3,chr6]\""
RESUME=""
if [ -f $R/last.ckpt ] && [ ! -f $R/periodic/epoch399.ckpt ]; then
  RESUME="ckpt_path=$R/last.ckpt"; mv $R/csv $R/csv_before_resume_$(date +%s); echo "resuming"
fi
if [ ! -f $R/periodic/epoch399.ckpt ]; then
  echo "=== train $(date)"
  eval $PIXI python scripts/train.py hydra.run.dir=$R wandb_logger.save_dir=/scratch/users/diamant/wandb \
      wandb_logger.name=$N checkpoint_every_n_epochs=25 log_mask_rate_bins=true assay_type=fiberseq_m6a \
      $DATA data.cfg.num_workers=7 \
      model.cfg.cnn_config.dna_encoder_type=alphagenome model.cfg.cnn_config.smf_body_type=unet \
      model.cfg.cnn_config.unet_norm=read model.cfg.cnn_config.unet_attention_layers=2 \
      model.cfg.cnn_config.unet_cross_attention=true \
      alphagenome_weights=/scratch/users/diamant/alphagenome_pytorch/model_all_folds.safetensors \
      model.cfg.optimizer=muon model.cfg.muon_lr=0.01 ema_decay=0.999 \
      model.cfg.cnn_config.unet_kernel_size=9 "model.cfg.cnn_config.unet_down_dilations=[4,4,2,1,1,1,1]" \
      model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true \
      model.cfg.cnn_config.output_conv_matmul=true $RESUME >> $R.log 2>&1
  echo "=== train exit $? $(date)"
fi
BEST=$(ls $R/epoch*-step*.ckpt | head -1)
echo "best $(basename "$BEST")"
out=$R/samples_val_subset1000
eval $PIXI python scripts/sample.py hydra.run.dir=$out "checkpoint='$BEST'" \
    $DATA data.cfg.max_regions_per_split=1000 regions.source=split regions.split=val \
&& $PIXI python scripts/metrics.py hydra.run.dir=${out/samples_/metrics_} samples_path=$out/samples.h5 \
    metrics.include_positional_baselines=false metrics.include_flank_baselines=false \
|| echo "FAILED val sample"
echo "=== all done $(date)"
