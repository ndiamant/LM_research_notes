#!/bin/bash
#SBATCH --job-name=smfnet_attn_muon_ema_k9
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --time=1-12:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_kernel/train_eval_attn_muon_ema_k9_%j.out

# Kernel-size-9 attention U-Net + Muon + EMA, 400 epochs, then evaluation;
# see train_eval_attn_muon_ema_k9.sh (~14 h training + ~2.5 h evaluation).
bash $SCRATCH/smfnet_kernel/train_eval_attn_muon_ema_k9.sh
