#!/bin/bash
#SBATCH --job-name=inpaint_ref_e49
#SBATCH --partition=owners
#SBATCH -G 1
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=/scratch/users/diamant/smfnet_inpaint/score_ref_e49_%j.out

# Fixed-mask val loss and span inpainting of the original k=9 run's epoch-49
# checkpoint on an H100, so it uses the same (GPU-drawn) masks and spans as
# the 50-epoch variants: the run-to-run noise between two runs of one config.
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
export PYTHONPATH=$HOME/repos/SMFNet_inpaint_runs/src
S=$SCRATCH/smfnet_inpaint
CKPT=/scratch/users/diamant/smf_models/ONT_mESC_alphagenome_unet_attn_muon_ema_k9/periodic/epoch=49.ckpt
nvidia-smi --query-gpu=name --format=csv,noheader
$PIXI python $SCRATCH/smfnet_unet/rescore_val.py --name k9_full_run_e49_h100 --checkpoint $CKPT --out $S/rescored_e49_k9_full_run.json
$PIXI python $S/span_infill.py --checkpoint $CKPT --out $S/pilot_e49_k9_full_run.json 2>&1 | grep "^span"
