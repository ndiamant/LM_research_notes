#!/bin/bash
# Span-inpainting pilot (span_infill.py) on 1,000 val active regions for
# several models, one after another on this session's GPU.
# Run with: bash run_pilot.sh
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
M=/scratch/users/diamant/smf_models
S=$SCRATCH/smfnet_inpaint
P=$HOME/repos/SMFNet_patch/src
A=$SCRATCH/src/SMFNet_atac/src
BW=/scratch/groups/btrippe/ndiamant/GSE184470_ESC_TKO_ATAC/tracks/ESC_TKO_pooled.tn5_counts.bw
run() {  # name, code, checkpoint, extra args...
  local name=$1 code=$2 ckpt=$3; shift 3
  echo "== $name $(date +%T)"
  PYTHONPATH=$code $PIXI python $S/span_infill.py --checkpoint "$ckpt" --out $S/pilot_$name.json "$@" 2>&1 \
      | grep -v "Warning\|it/s\]\|^batch" || echo "FAILED $name"
}
run k5 $P "$M/ONT_mESC_alphagenome_unet_attn_muon_ema/epoch=377-step=171990.ckpt"
run k9 $P "$M/ONT_mESC_alphagenome_unet_attn_muon_ema_k9/epoch=398-step=181545.ckpt"
run original $P "$M/SMFNet_full_run_2026-08-11/epoch=377-step=171990-converted.ckpt"
run ag_frozen $P "$M/ONT_mESC_alphagenome_frozen/epoch=377-step=171990.ckpt"
run p4_muon_ema $P "$M/ONT_mESC_unet_attn_p4_wide_h2_muon_ema/epoch=399-step=182000.ckpt"
run k9_1down $P "$M/ONT_mESC_alphagenome_unet_attn_muon_ema_r64_k9d1u5/epoch=199-step=182000.ckpt"
run atac $A "$M/ONT_mESC_ATAC_TKO_pooled/epoch=377-step=171990.ckpt" --atac-bigwig $BW
echo "all done $(date +%T)"
