#!/bin/bash
#SBATCH --job-name=smfnet_eval_new_runs
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/long_range/eval_new_runs_%j.out
# After the two 16 kb additive runs with strong jitter and short batches (jobs 52499650,
# 52499652): the long-range suite on their best checkpoints. Pair scores and
# long_range_eval for each; the short-range test with all six models (2 kb additive
# reference, original 2 kb, 16 kb pilot, 16 kb additive, and the two new runs).
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_${SLURM_JOB_ID}
CODE=$HOME/repos/SMFNet_unet_figwig_eval  # detached at e3f7032
export PYTHONPATH=$CODE/src
cd $CODE && git log -1 --oneline; hostname
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml"
S=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_
F=$SCRATCH/fiberseq/GM12878; L=$SCRATCH/smfnet_figwig/long_range
FASTA=/scratch/groups/btrippe/ndiamant/reference/hg38.fa
best() { ls $1/epoch*-step*.ckpt | head -1; }
JM=$(best ${S}16kb_flank1024_additive_emission_jitter_mixed_unet_attn_muon_ema_k9_dilated)
JB=$(best ${S}16kb_flank1024_additive_emission_jitter_mixed_bg10_unet_attn_muon_ema_k9_dilated)
echo "checkpoints: $JM $JB"
for pair in "jitter_mixed16:$JM" "jitter_mixed_bg10_16:$JB"; do
  name=${pair%%:*}; ckpt=${pair#*:}
  echo "=== pairs $name $(date)"
  $PIXI python scripts/score_cre_pairs.py hydra.run.dir=$L/scores_val_$name checkpoint=$ckpt \
    pairs_path=$L/cre_pairs_val/cre_pairs.tsv figwig_dir=$F/figwig efficiency_dir=$F/efficiency fasta_path=$FASTA
  echo "=== long_range_eval $name $(date)"
  $PIXI python scripts/long_range_eval.py hydra.run.dir=$L/eval_val_$name checkpoint=$ckpt \
    data.fasta_path=$FASTA data.figwig_dir=$F/figwig data.regions_bed=$F/kept_windows.bed \
    data.efficiency_dir=$F/efficiency data.region_length=18498 data.model_input_length=18432 \
    data.dna_flank=1024 data.min_read_fraction=0.75 data.min_overlap=500 data.min_spanning_molecules=32 \
    "data.val_chrm=[chr8,chr20]" spread.regions_per_batch=8
done
echo "=== same sites $(date)"
$PIXI python scripts/score_same_sites.py hydra.run.dir=$L/same_sites_val_6models reference=add2 \
  "checkpoints=[add2:${S}additive_emission_alphagenome_unet_attn_muon_ema_k9_dilated/epoch398-step187530.ckpt,orig2:${S}alphagenome_unet_attn_muon_ema_k9_dilated/epoch399-step188000.ckpt,pilot16:${S}16kb_flank1024_alphagenome_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt,add16:${S}16kb_flank1024_additive_emission_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt,jm16:$JM,jmbg16:$JB]" \
  data.fasta_path=$FASTA data.figwig_dir=$F/figwig data.regions_bed=$F/kept_windows.bed \
  data.efficiency_dir=$F/efficiency data.min_overlap=500 data.min_spanning_molecules=32 \
  "data.val_chrm=[chr8,chr20]" cre_bed=$F/inputs/ENCFF614SMH.bed.gz
echo "=== all done $(date)"
