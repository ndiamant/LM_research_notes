#!/bin/bash
# Three alternating baseline/container benchmark runs, sequential on the one GPU.
set -uo pipefail
B=/scratch/users/diamant/2026-11-23_jax_container_bench
RUN=/scratch/users/diamant/2026-11-17_fiber_hist190_tf60_ol61_kl300_alphagenome_bf16
ARGS="--run-dir $RUN --init-params $B/runs/init.npz --region-subset 40000 --steps 280 --pure-steps 200"
module load devel pixi/0.53.0
export PIXI_CACHE_DIR=$SCRATCH/pixi-cache JAX_LOG_COMPILES=1 WANDB_MODE=disabled
for rep in 1 2 3; do
  echo "[$(date +%T)] baseline rep$rep"
  (cd /home/users/diamant/repos/smf_hmm/.claude/worktrees/pixi-migration && \
   pixi run --frozen python $B/bench_train_step.py $ARGS --out $B/runs/baseline_rep$rep.json) \
   > $B/runs/baseline_rep$rep.log 2>&1
  tail -1 $B/runs/baseline_rep$rep.log
  echo "[$(date +%T)] container rep$rep"
  (cd /home/users/diamant/repos/smf_hmm/.claude/worktrees/jax-container-bench && \
   $B/ctr.sh env JAX_LOG_COMPILES=1 WANDB_MODE=disabled \
   pixi run --frozen python $B/bench_train_step.py $ARGS --out $B/runs/container_rep$rep.json) \
   > $B/runs/container_rep$rep.log 2>&1
  tail -1 $B/runs/container_rep$rep.log
done
echo "[$(date +%T)] all done"
