#!/bin/bash
#SBATCH --job-name=smfnet_inpaint_variant
#SBATCH --partition=owners
#SBATCH -G 1
#SBATCH -C GPU_SKU:H100_SXM5|GPU_SKU:H200_SXM5
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=/scratch/users/diamant/smfnet_inpaint/train_%x_%j.out

# One 50-epoch inpainting variant (owners: preemptible, so requeued; the
# script resumes from the run's last.ckpt). Submit with:
#   sbatch --job-name=inpaint_<variant> train_variant.sbatch <variant>
bash $SCRATCH/smfnet_inpaint/train_variant.sh "$1"
