#!/bin/bash
#SBATCH --job-name=smfnet_ctcf_additive_b0
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/ctcf_pairs/ctcf_additive_b0_%j.out
# Two inserted 3x CTCF elements (convergent, divergent, tandem, each alone; 2/5/8/12 kb), as
# the pilot's run50, with the 16 kb additive-emission model at fixed read b_m = 0 (no
# labelling variation between reads). Then the summary and figure.
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_$SLURM_JOB_ID
export PYTHONPATH=$HOME/repos/SMFNet_unet_figwig_runs/src  # detached at 759ba00
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml"
D=$SCRATCH/smfnet_figwig/ctcf_pairs; OUT=$D/run50_additive_b0; hostname
$PIXI python $D/ctcf_insertion.py \
  --checkpoint /scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_16kb_flank1024_additive_emission_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt \
  --out $OUT --efficiency 0 \
&& $PIXI python $D/analyse_ctcf_insertion.py $OUT > $OUT/summary.txt 2>&1 && cat $OUT/summary.txt
