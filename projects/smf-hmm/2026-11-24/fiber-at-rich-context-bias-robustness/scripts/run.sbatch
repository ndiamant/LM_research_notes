#!/bin/bash
#SBATCH --job-name=heldout_tf
#SBATCH --partition=btrippe,akundaje,owners
#SBATCH --requeue
#SBATCH -G 1
#SBATCH -C GPU_SKU:L40S|GPU_SKU:H100_SXM5|GPU_SKU:A100_SXM4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
set -euo pipefail
cd /scratch/users/diamant/fiber_eval_work/describe_outputs/heldout_inside_tf
RATE=0.2 /home/groups/btrippe/diamant/miniforge/envs/hmm/bin/python heldout_inside_tf.py
