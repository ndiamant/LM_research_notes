#!/bin/bash
#SBATCH --job-name=smfnet_marg_additive
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --array=0-2
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/marg/marginalize_additive_%A_%a.out
# TF marginalizations of the 2 kb additive-emission model (Muon + EMA) at fixed read
# efficiency b_m = 0 (median), -1.0 (10th percentile), +0.6 (90th), one per array task:
# the HMM repo's GM12878 panel (gm12878_hmm_tfs.tsv; make_motifs.py), 200 random
# backgrounds, both orientations, 1,024 reads, recording P(m6A) with and without the
# context and efficiency offsets (sampling.record_debiased).
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_${SLURM_JOB_ID}
export CODE=$HOME/repos/SMFNet_unet_figwig_runs PYTHONPATH=$HOME/repos/SMFNet_unet_figwig_runs/src  # detached at 759ba00
cd $CODE && git log -1 --oneline; hostname
B=(0.0 -1.0 0.6); TAG=(b0 bp10 bp90)
b=${B[$SLURM_ARRAY_TASK_ID]}; tag=${TAG[$SLURM_ARRAY_TASK_ID]}
M=/scratch/users/diamant/smf_models
pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml python scripts/marginalize.py \
  hydra.run.dir=$SCRATCH/smfnet_figwig/marg/additive_$tag \
  checkpoint=$M/GM12878_fiberseq_figwig_additive_emission_alphagenome_unet_attn_muon_ema_k9_dilated/epoch398-step187530.ckpt \
  motifs.path=$SCRATCH/smfnet_figwig/marg/gm12878_hmm_tfs.tsv \
  sampling.efficiency=$b sampling.record_debiased=true sampling.progress=false plot.enabled=true
