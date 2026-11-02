#!/bin/bash
#SBATCH --job-name=smfnet_ag_ftall
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_ftall_run2_%j.out

# Fully fine-tuned AlphaGenome DNA encoder (pretrained weights at 0.1x lr),
# otherwise the same as ONT_mESC_alphagenome_frozen. Code: nd_alpha_genome.
set -euo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
# Print a Python traceback if a process (e.g. a DataLoader worker) segfaults;
# run 1 (job 50228093) died from a worker segfault in epoch 17.
export PYTHONFAULTHANDLER=1
cd $HOME/repos/SMFNet_public
git log -1 --oneline
R=/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_ftall_run2
pixi run --frozen python scripts/train.py hydra.run.dir=$R \
    wandb_logger.save_dir=/scratch/users/diamant/wandb wandb_logger.name=ONT_mESC_alphagenome_ftall_run2 \
    assay_type=cg_gc_accessibility \
    "data.cfg.h5_path=[/scratch/groups/btrippe/ndiamant/ONT_mESC_acRegions.h5,/scratch/groups/btrippe/ndiamant/ONT_mESC_bgRegions.h5]" \
    data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/mm10.fa \
    model.cfg.cnn_config.dna_encoder_type=alphagenome \
    model.cfg.cnn_config.alphagenome_frozen=false \
    model.cfg.alphagenome_lr_mult=0.1 \
    log_mask_rate_bins=true \
    alphagenome_weights=/scratch/users/diamant/alphagenome_pytorch/model_all_folds.safetensors
