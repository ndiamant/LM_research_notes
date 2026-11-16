#!/bin/bash
#SBATCH --job-name=smfnet_same_sites
#SBATCH --partition=btrippe,hns,owners
#SBATCH --requeue
#SBATCH -C GPU_SKU:H100_SXM5
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/long_range/same_sites_%j.out
# Short-range paired test (score_same_sites.py): central 1 kb of the 2 kb val reads, the
# 2 kb models at 2,048 bp of DNA and the 16 kb models at 3,072 bp (outside their training
# distribution: they never trained on short inputs). Reference: the 2 kb additive model.
set -uo pipefail
export PATH=/share/software/user/open/pixi/0.53.0/bin:$PATH PIXI_CACHE_DIR=$SCRATCH/pixi-cache PYTHONFAULTHANDLER=1
export NUMBA_CACHE_DIR=${L_SCRATCH:-/tmp}/numba_cache_${SLURM_JOB_ID}
CODE=$HOME/repos/SMFNet_unet_figwig_eval  # detached at 515ea75
export PYTHONPATH=$CODE/src
cd $CODE && git log -1 --oneline; hostname
M=/scratch/users/diamant/smf_models/GM12878_fiberseq_figwig_
F=$SCRATCH/fiberseq/GM12878
pixi run --frozen --manifest-path $HOME/repos/SMFNet_unet_figwig/pyproject.toml python scripts/score_same_sites.py \
  hydra.run.dir=$SCRATCH/smfnet_figwig/long_range/same_sites_val reference=add2 \
  "checkpoints=[add2:${M}additive_emission_alphagenome_unet_attn_muon_ema_k9_dilated/epoch398-step187530.ckpt,orig2:${M}alphagenome_unet_attn_muon_ema_k9_dilated/epoch399-step188000.ckpt,pilot16:${M}16kb_flank1024_alphagenome_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt,add16:${M}16kb_flank1024_additive_emission_unet_attn_muon_ema_k9_dilated/epoch24-step91275.ckpt]" \
  data.fasta_path=/scratch/groups/btrippe/ndiamant/reference/hg38.fa data.figwig_dir=$F/figwig \
  data.regions_bed=$F/kept_windows.bed data.efficiency_dir=$F/efficiency data.min_overlap=500 \
  data.min_spanning_molecules=32 "data.val_chrm=[chr8,chr20]" cre_bed=$F/inputs/ENCFF614SMH.bed.gz
