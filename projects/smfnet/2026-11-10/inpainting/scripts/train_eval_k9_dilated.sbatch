#!/bin/bash
#SBATCH --job-name=smfnet_k9_dilated
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --time=1-12:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_dilated/train_eval_k9_dilated_%j.out

# k=9 + dilated fine down convs, Muon + EMA, 400 epochs, then evaluation;
# see train_eval_k9_dilated.sh (~14 h training + ~3 h evaluation).
bash $SCRATCH/smfnet_dilated/train_eval_k9_dilated.sh
