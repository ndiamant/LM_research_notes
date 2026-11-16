#!/bin/bash
#SBATCH --job-name=smfnet_long_range_eval
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --array=0-1
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/long_range/long_range_eval_%A_%a.out
# long_range_eval.py (read spread on 300 val regions, position invariance on 500 active
# elements, loss along the window on 300 regions) for one 16 kb checkpoint per task:
# 0 the pilot, 1 the additive-emission run.
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_${SLURM_JOB_ID}
CODE=$HOME/repos/SMFNet_unet_figwig_eval  # detached at a9c6dec
export PYTHONPATH=$CODE/src
cd $CODE && git log -1 --oneline; hostname
M=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_16kb_flank1024_
F=$SCRATCH/fiberseq/GM12878
NAMES=(pilot16 additive16)
CKPTS=(${M}alphagenome_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt
       ${M}additive_emission_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt)
i=$SLURM_ARRAY_TASK_ID
pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml python scripts/long_range_eval.py \
  hydra.run.dir=$SCRATCH/smfnet_figwig/long_range/eval_val_${NAMES[$i]} checkpoint=${CKPTS[$i]} \
  data.fasta_path=/scratch/groups/btrippe/ndiamant/reference/hg38.fa data.figwig_dir=$F/figwig \
  data.regions_bed=$F/kept_windows.bed data.efficiency_dir=$F/efficiency \
  data.region_length=18498 data.model_input_length=18432 data.dna_flank=1024 \
  data.min_read_fraction=0.75 data.min_overlap=500 data.min_spanning_molecules=32 \
  "data.val_chrm=[chr8,chr20]" spread.regions_per_batch=8
