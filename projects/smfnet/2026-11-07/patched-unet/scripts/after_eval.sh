#!/bin/bash
# After the AdamW+EMA eval: rescore both EMA runs' final (epoch 399) weights
# on fixed masks, test and benchmark the head-block variant, then train it.
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
M=/scratch/users/diamant/smf_models
S=$SCRATCH/smfnet_unet
for n in muon_ema adamw_ema; do
  PYTHONPATH=$HOME/repos/SMFNet_unet/src $PIXI python $S/rescore_val.py --name attn_${n}_e399 \
      --checkpoint "$M/ONT_mESC_alphagenome_unet_attn_$n/periodic/epoch=399.ckpt" --out $S/rescored/attn_${n}_e399.json || echo "FAILED rescore $n"
done
cd $HOME/repos/SMFNet_patch
PYTHONPATH=$PWD/src $PIXI python -m pytest -q tests/test_fast_unet.py tests/test_unet.py 2>&1 | tail -3
PYTHONPATH=$PWD/src $PIXI python $SCRATCH/smfnet_patch/bench_variants.py p4_wide_h2 2>&1 | grep RESULT
bash $SCRATCH/smfnet_patch/run_variants2.sh p4_wide_h2
echo all done
