#!/bin/bash
#SBATCH --job-name=smfnet_memorization
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_memorization/run_all_%j.out

# Train vs val memorization diagnostics; see run_all.sh.
bash $SCRATCH/smfnet_memorization/run_all.sh
