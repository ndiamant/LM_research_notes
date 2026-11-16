#!/bin/bash
#SBATCH --job-name=smfnet_additive_2kb
#SBATCH --partition=btrippe,hns
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --time=1-12:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/train_additive_2kb_%j.out

# The 2 kb figwig run with the additive emission head, 400 epochs (~14 h on an H100, as the
# original), then a val sample; see train_additive_2kb.sh.
bash $SCRATCH/smfnet_figwig/train_additive_2kb.sh
