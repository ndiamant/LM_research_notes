#!/bin/bash
# Train vs val memorization diagnostics (memorization_gap.py) for the main
# models, on 4,000 regions per h5 per split, plus the region-mean ceiling.
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
M=/scratch/users/diamant/smf_models
H=/scratch/groups/btrippe/ndiamant
S=$SCRATCH/smfnet_memorization
O=$S/results
P=$HOME/repos/SMFNet_patch/src
A=$SCRATCH/src/SMFNet_atac/src
BW=$H/GSE184470_ESC_TKO_ATAC/tracks/ESC_TKO_pooled.tn5_counts.bw
hostname; nvidia-smi --query-gpu=name --format=csv,noheader
run() {  # name, code, checkpoint ("" for the ceiling), h5 kind, split, extra args...
  local name=$1 code=$2 ckpt=$3 kind=$4 split=$5; shift 5
  local ck=(); [ -n "$ckpt" ] && ck=(--checkpoint "$ckpt")
  echo "== $name $kind $split $(date +%T)"
  PYTHONPATH=$code $PIXI python $S/memorization_gap.py ${ck[@]+"${ck[@]}"} --h5 $H/ONT_mESC_${kind}Regions.h5 \
      --split $split --out $O/${name}_${kind}_${split}.json "$@" || echo "FAILED $name $kind $split"
}
# Smoke test first: fail fast on a small subset.
run smoke_ceiling $P "" ac train --max-regions 64
run smoke_model $P "$M/ONT_mESC_alphagenome_frozen/epoch=377-step=171990.ckpt" ac train --max-regions 64
run smoke_atac $A "$M/ONT_mESC_ATAC_TKO_pooled/epoch=377-step=171990.ckpt" ac train --max-regions 64 --atac-bigwig $BW
for kind in ac bg; do
  for split in train val; do
    run ceiling $P "" $kind $split
    run baseline $P "$M/SMFNet_full_run_2026-08-11/epoch=377-step=171990-converted.ckpt" $kind $split
    run ag_frozen $P "$M/ONT_mESC_alphagenome_frozen/epoch=377-step=171990.ckpt" $kind $split
    run ag_ftall $P "$M/ONT_mESC_alphagenome_ftall_run2/epoch=377-step=171990.ckpt" $kind $split
    run attn_muon_ema $P "$M/ONT_mESC_alphagenome_unet_attn_muon_ema/epoch=377-step=171990.ckpt" $kind $split
    run atac $A "$M/ONT_mESC_ATAC_TKO_pooled/epoch=377-step=171990.ckpt" $kind $split --atac-bigwig $BW
  done
done
echo "all done $(date +%T)"
