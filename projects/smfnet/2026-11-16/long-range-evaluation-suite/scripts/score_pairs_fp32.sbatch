#!/bin/bash
#SBATCH --job-name=smfnet_score_pairs_fp32
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:30:00
#SBATCH --array=0-1
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/long_range/score_pairs_fp32_%A_%a.out
# Paired CRE-pair likelihood scores (score_cre_pairs.py) on the val pair table, one 16 kb
# checkpoint per task: 0 the pilot, 1 the additive-emission run.
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_${SLURM_JOB_ID}
CODE=$HOME/repos/SMFNet_unet_figwig_eval  # detached at a9c6dec (cre_pairs code as at 92423ca); float32 check of the bf16 scores
export PYTHONPATH=$CODE/src
cd $CODE && git log -1 --oneline; hostname
M=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_16kb_flank1024_
F=$SCRATCH/fiberseq/GM12878; L=$SCRATCH/smfnet_figwig/long_range
NAMES=(pilot16 additive16)
CKPTS=(${M}alphagenome_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt
       ${M}additive_emission_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt)
i=$SLURM_ARRAY_TASK_ID
pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml python scripts/score_cre_pairs.py \
  hydra.run.dir=$L/scores_val_fp32_${NAMES[$i]} checkpoint=${CKPTS[$i]} pairs_path=$L/cre_pairs_val/cre_pairs.tsv \
  figwig_dir=$F/figwig efficiency_dir=$F/efficiency fasta_path=/scratch/groups/btrippe/ndiamant/reference/hg38.fa dtype=float32
