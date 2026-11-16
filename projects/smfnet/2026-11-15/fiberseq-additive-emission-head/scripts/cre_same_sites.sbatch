#!/bin/bash
#SBATCH --job-name=smfnet_cre_same_sites
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/additive_eval/cre_same_sites_%j.out
# Central-1 kb CRE vs other-site loss: 16 kb additive vs 2 kb additive (+ pilot, original 2 kb).
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_$SLURM_JOB_ID
export CODE=$HOME/repos/SMFNet_unet_figwig_runs PYTHONPATH=$HOME/repos/SMFNet_unet_figwig_runs/src  # detached at 759ba00
cd $CODE && git log -1 --oneline; hostname
M=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_
pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml python $SCRATCH/smfnet_figwig/additive_eval/cre_same_sites.py \
  --add16 ${M}16kb_flank1024_additive_emission_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt \
  --pilot16 ${M}16kb_flank1024_alphagenome_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt \
  --add2 ${M}additive_emission_alphagenome_unet_attn_muon_ema_k9_dilated/epoch398-step187530.ckpt \
  --orig2 ${M}alphagenome_unet_attn_muon_ema_k9_dilated/epoch399-step188000.ckpt \
  --out $SCRATCH/smfnet_figwig/additive_eval/cre_same_sites.npz
