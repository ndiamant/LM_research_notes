#!/bin/bash
# Sample the same 1,000 GM12878 val regions (figwig loader, the H5's settings) with the
# new U-Net and the old CNN, then metrics.py on each (no positional baselines yet).
# Run with: sbatch eval_val_samples.sbatch
set -uo pipefail
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml"
export CODE=$HOME/repos/SMFNet_unet_figwig_runs  # detached at 6ce74dd
cd $CODE && git log -1 --oneline; hostname; nvidia-smi --query-gpu=name --format=csv,noheader
export PYTHONPATH=$CODE/src
M=/scratch/users/diamant/smf_models
F=$SCRATCH/fiberseq/GM12878
FASTA=/scratch/groups/btrippe/brian/projects/SMF/fiberseq/reference/hg38.fa
DATA="data=figwig data.cfg.fasta_path=$FASTA \
  data.cfg.figwig_dir=$F/figwig data.cfg.regions_bed=$F/kept_windows.bed \
  data.cfg.region_length=2114 data.cfg.model_input_length=2048 \
  data.cfg.min_overlap=500 data.cfg.min_spanning_molecules=32 \
  \"data.cfg.train_chrm=[chr2,chr4,chr5,chr7,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr21,chr22,chrX]\" \
  \"data.cfg.val_chrm=[chr8,chr20]\" \"data.cfg.test_chrm=[chr1,chr3,chr6]\" \
  data.cfg.max_regions_per_split=1000"
OUT=$SCRATCH/smfnet_figwig/eval
mkdir -p $OUT
for spec in \
  "new_unet_figwig=$M/GM12878_fiberseq_figwig_alphagenome_unet_attn_muon_ema_k9_dilated/epoch399-step188000.ckpt" \
  "old_cnn_h5=$SCRATCH/smfnet_figwig/GM12878_fiberseq_epoch369_ont_polarity_plainkeys.ckpt"; do
  name=${spec%%=*}; ckpt=${spec#*=}
  echo "=== $name $(date)"
  eval $PIXI python scripts/sample.py hydra.run.dir=$OUT/samples_val1000_$name "checkpoint='$ckpt'" \
      $DATA regions.source=split regions.split=val \
  && $PIXI python scripts/metrics.py hydra.run.dir=$OUT/metrics_val1000_$name \
      samples_path=$OUT/samples_val1000_$name/samples.h5 plot.enabled=true \
      metrics.include_positional_baselines=false metrics.include_flank_baselines=false \
  || echo "FAILED $name"
done
echo "=== all done $(date)"
