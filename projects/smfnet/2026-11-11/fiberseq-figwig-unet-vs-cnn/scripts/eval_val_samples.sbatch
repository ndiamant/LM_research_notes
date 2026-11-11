#!/bin/bash
#SBATCH --job-name=smfnet_fiber_eval
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/eval_val_samples_%j.out

# 1,000 val regions x 1,024 reads from the new U-Net and the old CNN, then metrics;
# see eval_val_samples.sh.
bash $SCRATCH/smfnet_figwig/eval_val_samples.sh
