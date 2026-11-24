#!/bin/bash
# Run a command inside the Rocky Linux 9 container with the GPU and Sherlock filesystems.
# Usage: ctr.sh <command...>   (run from the jax-container-bench worktree)
exec apptainer exec --nv \
  -B /home/groups,/scratch,/home/users,/share/software,/lscratch \
  --env PIXI_CACHE_DIR=/scratch/users/diamant/pixi-cache \
  --env XDG_CACHE_HOME=/scratch/users/diamant/xdg-cache \
  --env PATH=/share/software/user/open/pixi/0.53.0/bin:/usr/local/bin:/usr/bin:/bin \
  /home/groups/btrippe/diamant/containers/rockylinux9.sif "$@"
