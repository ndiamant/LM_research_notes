#!/bin/bash
#SBATCH --job-name=smfnet_scan_ckpts
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --array=0-9
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/long_range/scan/scan_%A_%a.out
# CRE pair scores (noise masking) over training: the 16 kb pilot and the 16 kb additive run
# at the periodic checkpoints after epochs 5, 10, 15, 20, 25 (files epoch4..epoch24).
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_${SLURM_JOB_ID}
CODE=$HOME/repos/SMFNet_unet_figwig_eval  # detached at 8756ee6
export PYTHONPATH=$CODE/src
cd $CODE && git log -1 --oneline; hostname
M=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_16kb_flank1024_
F=$SCRATCH/fiberseq/GM12878; L=$SCRATCH/smfnet_figwig/long_range
NAMES=(pilot16 additive16)
DIRS=(${M}alphagenome_unet_attn_muon_ema_k9_dilated ${M}additive_emission_unet_attn_muon_ema_k9_dilated)
EPOCHS=(4 9 14 19 24)
m=$((SLURM_ARRAY_TASK_ID / 5)); e=${EPOCHS[$((SLURM_ARRAY_TASK_ID % 5))]}
pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml python scripts/score_cre_pairs.py \
  hydra.run.dir=$L/scan/${NAMES[$m]}_epoch$e checkpoint=${DIRS[$m]}/periodic/epoch$e.ckpt \
  pairs_path=$L/cre_pairs_val/cre_pairs.tsv figwig_dir=$F/figwig efficiency_dir=$F/efficiency \
  fasta_path=/scratch/groups/btrippe/ndiamant/reference/hg38.fa
