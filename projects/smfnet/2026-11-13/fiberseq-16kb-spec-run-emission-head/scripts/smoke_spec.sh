#!/bin/bash
# Smoke test of the next-run spec on chr21: 16 kb + 1,024 bp flanks, strong jitter (32,768 bp
# regions, 3 retries), efficiency + min_molecule_rate, CRE mask, 10% background, mixed short
# batches (50% here to exercise them), the emission head, and a warm start from the 16 kb pilot.
# Usage: bash smoke_spec.sh <tag> [regions per batch, 32]
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH
# figwig's numba cache on node-local disk (each worker loads it), and tracebacks on segfaults:
# the first smoke run lost all workers of one val loader to a segfault (not reproduced).
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache PYTHONFAULTHANDLER=1
cd $HOME/repos/SMFNet_unet_figwig && export PYTHONPATH=$PWD/src
TAG=$1; BS=${2:-32}; F=$SCRATCH/fiberseq/GM12878
PILOT=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_16kb_flank1024_alphagenome_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt
WANDB_MODE=offline PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True pixi run --frozen python scripts/train.py \
  hydra.run.dir=$SCRATCH/tm_inspect/smoke_spec_$TAG assay_type=fiberseq_m6a \
  data=figwig data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/reference/hg38.fa \
  data.cfg.figwig_dir=$F/figwig data.cfg.regions_bed=$F/kept_windows.bed \
  "data.cfg.train_chrm=[chr21]" "data.cfg.val_chrm=[chr21]" "data.cfg.test_chrm=[chr21]" \
  data.cfg.region_length=32768 data.cfg.model_input_length=18432 data.cfg.dna_flank=1024 \
  data.cfg.min_read_fraction=0.75 data.cfg.batch_size=$BS data.cfg.reads_per_region=32 data.cfg.num_workers=7 \
  data.cfg.jitter_retries=3 data.cfg.efficiency_dir=$F/efficiency data.cfg.min_molecule_rate=0.02 \
  data.cfg.cre_bed=$F/inputs/ENCFF614SMH.bed.gz data.cfg.background_fraction=0.1 \
  "data.cfg.short_smf_lengths=[1024,2048,4096,8192]" data.cfg.short_batch_fraction=0.5 \
  model.cfg.input_length=18432 model.cfg.cnn_config.dna_flank=1024 model.cfg.variable_input_length=true \
  model.cfg.emission_head=true model.cfg.cnn_config.input_channels=4 model.cfg.cnn_config.read_cond_features=1 \
  emission_prior=$F/enzyme/enzyme_prior.npz init_weights=$PILOT \
  model.cfg.cnn_config.dna_encoder_type=alphagenome model.cfg.cnn_config.smf_body_type=unet \
  model.cfg.cnn_config.unet_norm=read model.cfg.cnn_config.unet_attention_layers=2 \
  model.cfg.cnn_config.unet_cross_attention=true \
  alphagenome_weights=/scratch/users/diamant/alphagenome_pytorch/model_all_folds.safetensors \
  model.cfg.optimizer=muon model.cfg.muon_lr=0.01 ema_decay=0.999 \
  model.cfg.cnn_config.unet_kernel_size=9 "model.cfg.cnn_config.unet_down_dilations=[4,4,2,1,1,1,1]" \
  model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true \
  model.cfg.cnn_config.output_conv_matmul=true \
  trainer.max_epochs=1 +trainer.limit_train_batches=40 +trainer.limit_val_batches=4
