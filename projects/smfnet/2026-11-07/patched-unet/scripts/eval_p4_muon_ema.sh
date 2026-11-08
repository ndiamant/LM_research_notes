#!/bin/bash
# Patched U-Net + Muon + EMA (400 epochs): fixed-mask val loss (best and
# epoch 399, which are the same checkpoint here), the full 32-step eval, the
# 64-step subset, and the train/val memorization diagnostics.
# Run with: bash eval_p4_muon_ema.sh
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
M=/scratch/users/diamant/smf_models
H=/scratch/groups/btrippe/ndiamant
S=$SCRATCH/smfnet_unet
export CODE=$HOME/repos/SMFNet_patch
N=ONT_mESC_unet_attn_p4_wide_h2_muon_ema
CKPT="$M/$N/epoch=399-step=182000.ckpt"
cd $CODE && git log -1 --oneline
export PYTHONPATH=$CODE/src
$PIXI python $S/rescore_val.py --name p4_muon_ema --checkpoint "$CKPT" --out $S/rescored/p4_muon_ema.json || echo "FAILED rescore"
bash -x $S/eval_run_code.sh $N "epoch=399-step=182000.ckpt" epoch399 || echo "FAILED eval_run"
out=$M/$N/samples_val_ac_subset1000_steps64
$PIXI python scripts/sample.py hydra.run.dir=$out "checkpoint='$CKPT'" \
    regions.source=split regions.split=val sampling.sample_steps=64 \
    data.cfg.fasta_path=$H/mm10.fa "data.cfg.h5_path=[$H/ONT_mESC_acRegions.h5]" \
    data.cfg.max_regions_per_h5_per_split=1000 \
&& $PIXI python scripts/metrics.py hydra.run.dir=${out/samples_/metrics_} \
    samples_path=$out/samples.h5 metrics.include_positional_baselines=false metrics.include_flank_baselines=false \
|| echo "FAILED 64-step"
G=$SCRATCH/smfnet_memorization
for kind in ac bg; do for split in train val; do
  $PIXI python $G/memorization_gap.py --checkpoint "$CKPT" --h5 $H/ONT_mESC_${kind}Regions.h5 --split $split \
      --out $G/results/p4_muon_ema_${kind}_${split}.json || echo "FAILED memorization $kind $split"
done; done
echo all done
