#!/bin/bash
# GM12878 Fiber-seq, 16,384 bp SMF window with 2 x 1,024 bp DNA flanks (18,432 bp input),
# pilot of 25 epochs. The k=9 attention U-Net with dilated fine down convs, Muon (lr 0.01)
# + EMA 0.999 as the 2 kb run, at 32 regions x 32 reads (77 GB peak on an H100 in the smoke
# test). Regions are the H5's kept-window centres (kept_windows.bed) widened to 18,498 bp
# (66 bp of jitter, as at 2 kb); items draw reads covering >= 75% of the SMF window, and a
# region needs 32 of them. ChromBPNet fold 0. Resumes from last.ckpt if rerun.
# Run in an H100 allocation with: bash train_k9_dilated_fiber_16kb.sh
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml"
export CODE=$HOME/repos/SMFNet_unet_figwig_runs  # detached at 961d8ba
cd $CODE && git log -1 --oneline; hostname; nvidia-smi --query-gpu=name --format=csv,noheader
export PYTHONPATH=$CODE/src
M=/scratch/users/diamant/smf_models
F=$SCRATCH/fiberseq/GM12878
FASTA=/scratch/groups/btrippe/ndiamant/reference/hg38.fa
N=GM12878_fiberseq_figwig_16kb_flank1024_alphagenome_unet_attn_muon_ema_k9_dilated
R=$M/$N
EPOCHS=25
DATA="data=figwig data.cfg.fasta_path=$FASTA \
  data.cfg.figwig_dir=$F/figwig data.cfg.regions_bed=$F/kept_windows.bed \
  data.cfg.region_length=18498 data.cfg.model_input_length=18432 data.cfg.dna_flank=1024 \
  data.cfg.min_read_fraction=0.75 data.cfg.min_overlap=500 data.cfg.min_spanning_molecules=32 \
  \"data.cfg.train_chrm=[chr2,chr4,chr5,chr7,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr21,chr22,chrX]\" \
  \"data.cfg.val_chrm=[chr8,chr20]\" \"data.cfg.test_chrm=[chr1,chr3,chr6]\""
RESUME=""
if [ -f $R/last.ckpt ] && [ ! -f $R/periodic/epoch$((EPOCHS - 1)).ckpt ]; then
  RESUME="ckpt_path=$R/last.ckpt"; mv $R/csv $R/csv_before_resume_$(date +%s); echo "resuming"
fi
echo "=== train $(date)"
eval $PIXI python scripts/train.py hydra.run.dir=$R wandb_logger.save_dir=/scratch/users/diamant/wandb \
    wandb_logger.name=$N checkpoint_every_n_epochs=5 log_mask_rate_bins=true assay_type=fiberseq_m6a \
    $DATA data.cfg.batch_size=32 data.cfg.reads_per_region=32 data.cfg.num_workers=7 \
    model.cfg.input_length=18432 model.cfg.cnn_config.dna_flank=1024 \
    model.cfg.cnn_config.dna_encoder_type=alphagenome model.cfg.cnn_config.smf_body_type=unet \
    model.cfg.cnn_config.unet_norm=read model.cfg.cnn_config.unet_attention_layers=2 \
    model.cfg.cnn_config.unet_cross_attention=true \
    alphagenome_weights=/scratch/users/diamant/alphagenome_pytorch/model_all_folds.safetensors \
    model.cfg.optimizer=muon model.cfg.muon_lr=0.01 ema_decay=0.999 \
    model.cfg.cnn_config.unet_kernel_size=9 "model.cfg.cnn_config.unet_down_dilations=[4,4,2,1,1,1,1]" \
    model.cfg.cnn_config.unet_fused_norm=true model.cfg.cnn_config.unet_channels_last=true \
    model.cfg.cnn_config.output_conv_matmul=true trainer.max_epochs=$EPOCHS $RESUME >> $R.log 2>&1
echo "=== train exit $? $(date)"
