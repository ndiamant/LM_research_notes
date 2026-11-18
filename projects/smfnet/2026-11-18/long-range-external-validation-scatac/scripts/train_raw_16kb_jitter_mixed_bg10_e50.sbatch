#!/bin/bash
#SBATCH --job-name=smfnet_raw16_jm_bg10_e50
#SBATCH --partition=btrippe,hns
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/train_raw_16kb_jitter_mixed_bg10_e50_%j.out

# The 16 kb jitter + short batches + 10% background run as raw discrete diffusion (no
# emission head, no b_m), 50 epochs (~24 h on an H100 at the 25-epoch runs' ~11.8 h); see
# train_raw_16kb_jitter_mixed_bg10_e50.sh. Resumes from last.ckpt if resubmitted.
bash $SCRATCH/smfnet_figwig/train_raw_16kb_jitter_mixed_bg10_e50.sh
