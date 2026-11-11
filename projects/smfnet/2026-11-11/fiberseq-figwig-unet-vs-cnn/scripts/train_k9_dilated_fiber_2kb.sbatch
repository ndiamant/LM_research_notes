#!/bin/bash
#SBATCH --job-name=smfnet_fiber_figwig_k9
#SBATCH --partition=btrippe
#SBATCH -G 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --time=1-12:00:00
#SBATCH --output=/scratch/users/diamant/smfnet_figwig/train_k9_dilated_fiber_2kb_%j.out

# GM12878 Fiber-seq via the figwig loader on the H5's settings: k=9 dilated attention
# U-Net, Muon + EMA, 400 epochs (~14 h on an H100, as the mESC run's ~460 steps per epoch),
# then a val sample; see train_k9_dilated_fiber_2kb.sh.
bash $SCRATCH/smfnet_figwig/train_k9_dilated_fiber_2kb.sh
