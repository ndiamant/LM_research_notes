#!/bin/bash
# AdamW + EMA attention U-Net: fixed-mask val loss, full 32-step eval, and the
# 64-step subset, as for the Muon + EMA run. Run with: bash eval_adamw_ema.sh
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
M=/scratch/users/diamant/smf_models
H=/scratch/groups/btrippe/ndiamant
S=$SCRATCH/smfnet_unet
U=$HOME/repos/SMFNet_unet
N=ONT_mESC_alphagenome_unet_attn_adamw_ema
CKPT="$M/$N/epoch=396-step=180635.ckpt"
cd $U && git log -1 --oneline
for c in "attn_adamw_ema:$CKPT" "attn_adamw_ema_last:$M/$N/last.ckpt"; do
  PYTHONPATH=$U/src $PIXI python $S/rescore_val.py --name ${c%%:*} --checkpoint "${c#*:}" --out $S/rescored/${c%%:*}.json || echo "FAILED rescore ${c%%:*}"
done
bash -x $S/eval_run.sh $N "epoch=396-step=180635.ckpt" epoch396 || echo "FAILED eval_run"
out=$M/$N/samples_val_ac_subset1000_steps64
PYTHONPATH=$U/src $PIXI python scripts/sample.py hydra.run.dir=$out "checkpoint='$CKPT'" \
    regions.source=split regions.split=val sampling.sample_steps=64 \
    data.cfg.fasta_path=$H/mm10.fa "data.cfg.h5_path=[$H/ONT_mESC_acRegions.h5]" \
    data.cfg.max_regions_per_h5_per_split=1000 \
&& PYTHONPATH=$U/src $PIXI python scripts/metrics.py hydra.run.dir=${out/samples_/metrics_} \
    samples_path=$out/samples.h5 metrics.include_positional_baselines=false metrics.include_flank_baselines=false \
|| echo "FAILED 64-step"
echo all done
