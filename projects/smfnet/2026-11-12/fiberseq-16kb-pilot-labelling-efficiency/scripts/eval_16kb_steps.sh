#!/bin/bash
# Does the generated molecule-wide co-accessibility (the long-range plateau) depend on the
# number of sampling steps? The same 300 val regions sampled at 32, 64 and 128 steps
# (512 reads each), then compare_samples.py on each.
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml"
export CODE=$HOME/repos/SMFNet_unet_figwig_runs  # detached at 961d8ba
cd $CODE && export PYTHONPATH=$CODE/src
CKPT=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_16kb_flank1024_alphagenome_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt
F=$SCRATCH/fiberseq/GM12878
OUT=$SCRATCH/smfnet_figwig/eval16/steps
mkdir -p $OUT
DATA="data=figwig data.cfg.fasta_path=/scratch/groups/btrippe/ndiamant/reference/hg38.fa \
  data.cfg.figwig_dir=$F/figwig data.cfg.regions_bed=$F/kept_windows.bed \
  data.cfg.region_length=18498 data.cfg.model_input_length=18432 data.cfg.dna_flank=1024 \
  data.cfg.min_read_fraction=0.75 data.cfg.min_overlap=500 data.cfg.min_spanning_molecules=32 \
  \"data.cfg.train_chrm=[chr2,chr4,chr5,chr7,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr21,chr22,chrX]\" \
  \"data.cfg.val_chrm=[chr8,chr20]\" \"data.cfg.test_chrm=[chr1,chr3,chr6]\" \
  data.cfg.max_regions_per_split=300"
for STEPS in 32 64 128; do
  echo "=== $STEPS steps $(date)"
  eval $PIXI python scripts/sample.py hydra.run.dir=$OUT/samples_val300_steps$STEPS "checkpoint='$CKPT'" $DATA \
      regions.source=split regions.split=val sampling.reads_per_region=512 sampling.sample_steps=$STEPS \
      sampling.reads_per_chunk=128 sampling.regions_per_batch=8|| { echo "FAILED sample $STEPS"; continue; }
  $PIXI python $SCRATCH/smfnet_figwig/compare_samples.py --samples $OUT/samples_val300_steps$STEPS/samples.h5 \
      --out $OUT/compare_val300_steps$STEPS --n-examples 0 || echo "FAILED compare $STEPS"
done
echo "=== all done $(date)"
