#!/bin/bash
# Smoke test: 16,384 bp SMF window + 2 x 1,024 bp DNA flanks, min_read_fraction 0.75, on chr21.
# Usage: bash smoke_16kb.sh <regions per batch> <tag> [reads per region, 32]
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH
cd $HOME/repos/SMFNet_unet_figwig && export PYTHONPATH=$PWD/src
BS=$1; TAG=$2; RPR=${3:-32}; F=$SCRATCH/fiberseq/GM12878
WANDB_MODE=offline PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True pixi run --frozen python scripts/train.py \
  hydra.run.dir=$SCRATCH/tm_inspect/smoke_16kb_$TAG assay_type=fiberseq_m6a \
  data=figwig data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/reference/hg38.fa \
  data.cfg.figwig_dir=$F/figwig data.cfg.regions_bed=$F/kept_windows.bed \
  "data.cfg.train_chrm=[chr21]" "data.cfg.val_chrm=[chr21]" "data.cfg.test_chrm=[chr21]" \
  data.cfg.region_length=18498 data.cfg.model_input_length=18432 data.cfg.dna_flank=1024 \
  data.cfg.min_read_fraction=0.75 data.cfg.batch_size=$BS data.cfg.reads_per_region=$RPR data.cfg.num_workers=7 \
  model.cfg.input_length=18432 model.cfg.cnn_config.dna_flank=1024 \
  model.cfg.cnn_config.dna_encoder_type=alphagenome model.cfg.cnn_config.smf_body_type=unet \
  model.cfg.cnn_config.unet_norm=read model.cfg.cnn_config.unet_attention_layers=2 \
  model.cfg.cnn_config.unet_cross_attention=true \
  alphagenome_weights=/scratch/users/diamant/alphagenome_pytorch/model_all_folds.safetensors \
  model.cfg.optimizer=muon model.cfg.muon_lr=0.01 ema_decay=0.999 \
  model.cfg.cnn_config.unet_kernel_size=9 "model.cfg.cnn_config.unet_down_dilations=[4,4,2,1,1,1,1]" \
  model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true \
  model.cfg.cnn_config.output_conv_matmul=true \
  trainer.max_epochs=1 +trainer.limit_train_batches=40 +trainer.limit_val_batches=4
