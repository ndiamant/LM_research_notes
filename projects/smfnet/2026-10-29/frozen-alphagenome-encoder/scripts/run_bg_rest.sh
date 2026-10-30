#!/bin/bash
set -euo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
cd $HOME/repos/SMFNet_public
D=$SCRATCH/smfnet_ag_smoke/diagnostics
M=/scratch/users/diamant/smf_models
H=/scratch/groups/btrippe/ndiamant
BW=$H/GSE184470_ESC_TKO_ATAC/tracks/ESC_TKO_pooled.tn5_counts.bw
h5=$H/ONT_mESC_bgRegions.h5
pixi run --frozen python $D/mask_rate_site_loss.py --checkpoint "$M/ONT_mESC_alphagenome_frozen/epoch=377-step=171990.ckpt" --h5 $h5 --out $D/alphagenome_bg.json
PYTHONPATH=$SCRATCH/src/SMFNet_atac/src pixi run --frozen python $D/mask_rate_site_loss.py --checkpoint "$M/ONT_mESC_ATAC_TKO_pooled/epoch=377-step=171990.ckpt" --h5 $h5 --atac-bigwig $BW --out $D/atac_bg.json
