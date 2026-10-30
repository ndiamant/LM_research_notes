#!/bin/bash
set -euo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
cd $HOME/repos/SMFNet_public
D=$SCRATCH/smfnet_ag_smoke/diagnostics
M=/scratch/users/diamant/smf_models
H=/scratch/groups/btrippe/ndiamant
BW=$H/GSE184470_ESC_TKO_ATAC/tracks/ESC_TKO_pooled.tn5_counts.bw
while pgrep -f run_all.sh > /dev/null; do sleep 60; done
for kind in ac bg; do
  h5=$H/ONT_mESC_${kind}Regions.h5
  R="--mask-rates 0.99 1.0"
  pixi run --frozen python $D/mask_rate_site_loss.py --checkpoint "$M/SMFNet_full_run_2026-08-11/epoch=377-step=171990-converted.ckpt" --h5 $h5 $R --out $D/full_baseline_$kind.json
  pixi run --frozen python $D/mask_rate_site_loss.py --checkpoint "$M/ONT_mESC_alphagenome_frozen/epoch=377-step=171990.ckpt" --h5 $h5 $R --out $D/full_alphagenome_$kind.json
  PYTHONPATH=$SCRATCH/src/SMFNet_atac/src pixi run --frozen python $D/mask_rate_site_loss.py --checkpoint "$M/ONT_mESC_ATAC_TKO_pooled/epoch=377-step=171990.ckpt" --h5 $h5 --atac-bigwig $BW $R --out $D/full_atac_$kind.json
done
