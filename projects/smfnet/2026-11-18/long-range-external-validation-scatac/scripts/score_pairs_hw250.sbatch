#!/bin/bash
#SBATCH --job-name=smfnet_pairs_hw250
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/long_range/score_pairs_hw250_%j.out
# CRE pair scores (noise masking) of the four 16 kb models with +/- 250 bp elements, on the
# val pair table built with pairs.half_width=250 (cre_pairs_val_hw250), to match the
# scATAC +/- 250 bp pileups. A model already scored is skipped, so a requeue resumes.
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_${SLURM_JOB_ID}
CODE=$HOME/repos/SMFNet_scatac_runs  # detached at 52518fb
export PYTHONPATH=$CODE/src
cd $CODE && git log -1 --oneline; hostname
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml"
S=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_16kb_flank1024_
F=$SCRATCH/fiberseq/GM12878; L=$SCRATCH/smfnet_figwig/long_range
FASTA=/scratch/groups/btrippe/ndiamant/reference/hg38.fa
E=epoch24-step91275.ckpt
for pair in "pilot16:${S}alphagenome_unet_attn_muon_ema_k9_dilated/$E" \
  "additive16:${S}additive_emission_unet_attn_muon_ema_k9_dilated/$E" \
  "jitter_mixed16:${S}additive_emission_jitter_mixed_unet_attn_muon_ema_k9_dilated/$E" \
  "jitter_mixed_bg10_16:${S}additive_emission_jitter_mixed_bg10_unet_attn_muon_ema_k9_dilated/$E"; do
  name=${pair%%:*}; ckpt=${pair#*:}; out=$L/scores_val_hw250_$name
  if [ -f $out/pair_summary.json ]; then echo "=== $name done already"; continue; fi
  echo "=== pairs $name $(date)"
  $PIXI python scripts/score_cre_pairs.py hydra.run.dir=$out checkpoint=$ckpt \
    pairs_path=$L/cre_pairs_val_hw250/cre_pairs.tsv figwig_dir=$F/figwig efficiency_dir=$F/efficiency \
    fasta_path=$FASTA scoring.half_width=250
done
echo "=== all done $(date)"
