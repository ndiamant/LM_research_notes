#!/bin/bash
#SBATCH --job-name=smfnet_read_spread
#SBATCH --partition=btrippe,hns
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/additive_eval/read_spread_%j.out
# Read-mean spread and co-accessibility of generated reads: additive emission model (matched
# b_m, b_m = 0) vs the original 2 kb model, 1,000 val regions x 128 reads; see read_spread.py.
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_$SLURM_JOB_ID
export CODE=$HOME/repos/SMFNet_unet_figwig_runs PYTHONPATH=$HOME/repos/SMFNet_unet_figwig_runs/src  # detached at 2514917
cd $CODE && git log -1 --oneline; hostname
M=/scratch/users/diamant/smf_models
pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml python $SCRATCH/smfnet_figwig/additive_eval/read_spread.py \
  --additive $M/GM12878_fiberseq_figwig_additive_emission_alphagenome_unet_attn_muon_ema_k9_dilated/epoch398-step187530.ckpt \
  --original $M/GM12878_fiberseq_figwig_alphagenome_unet_attn_muon_ema_k9_dilated/epoch399-step188000.ckpt \
  --out $SCRATCH/smfnet_figwig/additive_eval/read_spread.npz
