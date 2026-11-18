#!/bin/bash
#SBATCH --job-name=smfnet_eval_raw16_e50
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:30:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/long_range/eval_raw16_e50_%j.out
# The long-range suite on the raw (no emission head, no b_m) 16 kb jitter + short batches +
# background run trained 50 epochs (job 52625990): its final checkpoint (epoch 49) and the
# periodic one after 25 epochs (mid-schedule). For each: CRE pair scores with +/- 250 bp
# elements (cre_pairs_val_hw250) and +/- 100 bp (cre_pairs_val), and long_range_eval; then
# the short-range test with the 2 kb additive reference, the pilot, the additive background
# run and both new checkpoints. A step with output already is skipped, so a requeue resumes.
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_${SLURM_JOB_ID}
CODE=$HOME/repos/SMFNet_scatac_runs  # detached at 52518fb
export PYTHONPATH=$CODE/src
cd $CODE && git log -1 --oneline; hostname
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml"
S=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_
RAW=${S}16kb_flank1024_jitter_mixed_bg10_e50_alphagenome_unet_attn_muon_ema_k9_dilated
F=$SCRATCH/fiberseq/GM12878; L=$SCRATCH/smfnet_figwig/long_range
FASTA=/scratch/groups/btrippe/ndiamant/reference/hg38.fa
for pair in "raw16_e50:$RAW/epoch49-step182550.ckpt" "raw16_e50_ep24:$RAW/periodic/epoch24.ckpt"; do
  name=${pair%%:*}; ckpt=${pair#*:}
  if [ ! -f $L/scores_val_hw250_$name/pair_summary.json ]; then
    echo "=== pairs hw250 $name $(date)"
    $PIXI python scripts/score_cre_pairs.py hydra.run.dir=$L/scores_val_hw250_$name checkpoint=$ckpt \
      pairs_path=$L/cre_pairs_val_hw250/cre_pairs.tsv figwig_dir=$F/figwig efficiency_dir=$F/efficiency \
      fasta_path=$FASTA scoring.half_width=250
  fi
  if [ ! -f $L/scores_val_$name/pair_summary.json ]; then
    echo "=== pairs hw100 $name $(date)"
    $PIXI python scripts/score_cre_pairs.py hydra.run.dir=$L/scores_val_$name checkpoint=$ckpt \
      pairs_path=$L/cre_pairs_val/cre_pairs.tsv figwig_dir=$F/figwig efficiency_dir=$F/efficiency fasta_path=$FASTA
  fi
  if [ ! -f $L/eval_val_$name/long_range_summary.json ]; then
    echo "=== long_range_eval $name $(date)"
    $PIXI python scripts/long_range_eval.py hydra.run.dir=$L/eval_val_$name checkpoint=$ckpt \
      data.fasta_path=$FASTA data.figwig_dir=$F/figwig data.regions_bed=$F/kept_windows.bed \
      data.efficiency_dir=$F/efficiency data.region_length=18498 data.model_input_length=18432 \
      data.dna_flank=1024 data.min_read_fraction=0.75 data.min_overlap=500 data.min_spanning_molecules=32 \
      "data.val_chrm=[chr8,chr20]" spread.regions_per_batch=8
  fi
done
echo "=== same sites $(date)"
$PIXI python scripts/score_same_sites.py hydra.run.dir=$L/same_sites_val_raw16_e50 reference=add2 \
  "checkpoints=[add2:${S}additive_emission_alphagenome_unet_attn_muon_ema_k9_dilated/epoch398-step187530.ckpt,pilot16:${S}16kb_flank1024_alphagenome_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt,jmbg16:${S}16kb_flank1024_additive_emission_jitter_mixed_bg10_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt,raw50:$RAW/epoch49-step182550.ckpt,raw50_ep24:$RAW/periodic/epoch24.ckpt]" \
  data.fasta_path=$FASTA data.figwig_dir=$F/figwig data.regions_bed=$F/kept_windows.bed \
  data.efficiency_dir=$F/efficiency data.min_overlap=500 data.min_spanning_molecules=32 \
  "data.val_chrm=[chr8,chr20]" cre_bed=$F/inputs/ENCFF614SMH.bed.gz
echo "=== all done $(date)"
