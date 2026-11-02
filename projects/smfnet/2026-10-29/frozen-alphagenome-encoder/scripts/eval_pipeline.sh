#!/bin/bash
# Sample val active regions, compute metrics, and run mask-rate diagnostics
# for one model. Usage: eval_pipeline.sh <run name>
set -euo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
export PYTHONFAULTHANDLER=1
cd $HOME/repos/SMFNet_public
NAME=$1
R=/scratch/users/diamant/smf_models/$NAME
CKPT="$R/epoch=377-step=171990.ckpt"
H=/scratch/groups/btrippe/ndiamant
D=$SCRATCH/smfnet_ag_smoke/diagnostics
git log -1 --oneline
pixi run --frozen python scripts/sample.py hydra.run.dir=$R/samples_val_ac_epoch377 \
    "checkpoint='$CKPT'" regions.source=split regions.split=val \
    data.cfg.fasta_path=$H/mm10.fa "data.cfg.h5_path=[$H/ONT_mESC_acRegions.h5]"
pixi run --frozen python scripts/metrics.py hydra.run.dir=$R/metrics_val_ac_epoch377_plots \
    samples_path=$R/samples_val_ac_epoch377/samples.h5 plot.enabled=true
for kind in ac bg; do
  pixi run --frozen python $D/mask_rate_site_loss.py --checkpoint "$CKPT" --h5 $H/ONT_mESC_${kind}Regions.h5 --out $D/${NAME}_$kind.json
  pixi run --frozen python $D/mask_rate_site_loss.py --checkpoint "$CKPT" --h5 $H/ONT_mESC_${kind}Regions.h5 --mask-rates 0.99 1.0 --out $D/full_${NAME}_$kind.json
done
echo pipeline done
