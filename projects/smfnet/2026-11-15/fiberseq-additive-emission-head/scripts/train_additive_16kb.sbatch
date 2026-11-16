#!/bin/bash
#SBATCH --job-name=smfnet_additive_16kb
#SBATCH --partition=btrippe,hns
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --time=1-12:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/train_additive_16kb_%j.out

# The 16 kb pilot with the additive emission head, 25 epochs (~11.5 h on an H100, as the
# pilot); see train_additive_16kb.sh. Resumes from last.ckpt if resubmitted.
bash $SCRATCH/smfnet_figwig/train_additive_16kb.sh
