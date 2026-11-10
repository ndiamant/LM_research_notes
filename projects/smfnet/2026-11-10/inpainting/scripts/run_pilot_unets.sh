#!/bin/bash
# Inpainting pilot for two more U-Nets that isolate the norm and attention:
# the plain U-Net (per-position norm, no attention, AdamW) and the attention
# U-Net with per-read norm (AdamW, no EMA). Run with: bash run_pilot_unets.sh
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
M=/scratch/users/diamant/smf_models
S=$SCRATCH/smfnet_inpaint
P=$HOME/repos/SMFNet_patch/src
for spec in "unet_plain:$M/ONT_mESC_alphagenome_frozen_unet/epoch=377-step=171990.ckpt" \
            "unet_attn_adamw:$M/ONT_mESC_alphagenome_unet_attn/epoch=366-step=166985.ckpt"; do
  name=${spec%%:*}; echo "== $name $(date +%T)"
  PYTHONPATH=$P $PIXI python $S/span_infill.py --checkpoint "${spec#*:}" --out $S/pilot_$name.json 2>&1 \
      | grep -v "Warning\|it/s\]\|^batch\|weights_only" || echo "FAILED $name"
done
echo "all done $(date +%T)"
