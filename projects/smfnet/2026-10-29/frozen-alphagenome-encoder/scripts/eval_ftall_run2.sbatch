#!/bin/bash
#SBATCH --job-name=smfnet_eval_ftall
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_ag_smoke/eval_ONT_mESC_alphagenome_ftall_run2_%j.out
# Sampling, metrics, and mask-rate diagnostics for the fully fine-tuned model.
$SCRATCH/smfnet_ag_smoke/eval_pipeline.sh ONT_mESC_alphagenome_ftall_run2
