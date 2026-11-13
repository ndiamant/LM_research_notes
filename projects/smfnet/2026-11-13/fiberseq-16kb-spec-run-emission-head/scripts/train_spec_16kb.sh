#!/bin/bash
# GM12878 Fiber-seq, the next-run spec (FIGWIG_TODO.md, 2026-11-12) at 16,384 bp + 2 x 1,024 bp
# flanks, 25 epochs from scratch, otherwise the 16 kb pilot's architecture and optimizer:
# - mixed lengths: 10% of batches at a 1, 2, 4 or 8 kb SMF window (equal cost), a val loader each;
# - strong jitter: 32,768 bp regions (+/-7 kb), offsets redrawn (3 tries) if < 32 reads cover 75%;
# - 10% background regions (weighted sampler);
# - per-molecule labelling efficiency (b_m), molecules below a 2% background m6A rate dropped;
# - the enzyme emission head (frozen +/-3 bp context table, b_m FiLM, e input channel);
# - val/loss_cre and val/loss_bg from the IDR peaks.
# Resumes from last.ckpt if rerun. Run in an H100 allocation: bash train_spec_16kb.sh
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
# figwig's numba cache on node-local disk: each data worker loads it.
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_$$
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml"
export CODE=$HOME/repos/SMFNet_unet_figwig_runs  # detached at 75549a4
cd $CODE && git log -1 --oneline; hostname; nvidia-smi --query-gpu=name --format=csv,noheader
export PYTHONPATH=$CODE/src
M=/scratch/users/diamant/smf_models
F=$SCRATCH/fiberseq/GM12878
FASTA=/scratch/groups/btrippe/ndiamant/reference/hg38.fa
N=GM12878_fiberseq_figwig_16kb_spec_mixed_jitter_emission_unet_attn_muon_ema_k9_dilated
R=$M/$N
EPOCHS=25
DATA="data=figwig data.cfg.fasta_path=$FASTA \
  data.cfg.figwig_dir=$F/figwig data.cfg.regions_bed=$F/kept_windows.bed \
  data.cfg.region_length=32768 data.cfg.model_input_length=18432 data.cfg.dna_flank=1024 \
  data.cfg.min_read_fraction=0.75 data.cfg.min_overlap=500 data.cfg.min_spanning_molecules=32 \
  data.cfg.jitter_retries=3 data.cfg.background_fraction=0.1 \
  \"data.cfg.short_smf_lengths=[1024,2048,4096,8192]\" data.cfg.short_batch_fraction=0.1 \
  data.cfg.efficiency_dir=$F/efficiency data.cfg.min_molecule_rate=0.02 \
  data.cfg.cre_bed=$F/inputs/ENCFF614SMH.bed.gz \
  \"data.cfg.train_chrm=[chr2,chr4,chr5,chr7,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr21,chr22,chrX]\" \
  \"data.cfg.val_chrm=[chr8,chr20]\" \"data.cfg.test_chrm=[chr1,chr3,chr6]\""
RESUME=""
if [ -f $R/last.ckpt ] && [ ! -f $R/periodic/epoch$((EPOCHS - 1)).ckpt ]; then
  RESUME="ckpt_path=$R/last.ckpt"; mv $R/csv $R/csv_before_resume_$(date +%s); echo "resuming"
fi
PRIOR="emission_prior=$F/enzyme/enzyme_prior.npz"
echo "=== train $(date)"
eval $PIXI python scripts/train.py hydra.run.dir=$R wandb_logger.save_dir=/scratch/users/diamant/wandb \
    wandb_logger.name=$N checkpoint_every_n_epochs=5 log_mask_rate_bins=true assay_type=fiberseq_m6a \
    $DATA data.cfg.batch_size=32 data.cfg.reads_per_region=32 data.cfg.num_workers=7 \
    model.cfg.input_length=18432 model.cfg.cnn_config.dna_flank=1024 model.cfg.variable_input_length=true \
    model.cfg.emission_head=true model.cfg.cnn_config.input_channels=4 \
    model.cfg.cnn_config.read_cond_features=1 $PRIOR \
    model.cfg.cnn_config.dna_encoder_type=alphagenome model.cfg.cnn_config.smf_body_type=unet \
    model.cfg.cnn_config.unet_norm=read model.cfg.cnn_config.unet_attention_layers=2 \
    model.cfg.cnn_config.unet_cross_attention=true \
    alphagenome_weights=/scratch/users/diamant/alphagenome_pytorch/model_all_folds.safetensors \
    model.cfg.optimizer=muon model.cfg.muon_lr=0.01 ema_decay=0.999 \
    model.cfg.cnn_config.unet_kernel_size=9 "model.cfg.cnn_config.unet_down_dilations=[4,4,2,1,1,1,1]" \
    model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true \
    model.cfg.cnn_config.output_conv_matmul=true trainer.max_epochs=$EPOCHS $RESUME >> $R.log 2>&1
echo "=== train exit $? $(date)"
