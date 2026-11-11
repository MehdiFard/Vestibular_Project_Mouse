#!/bin/bash
#SBATCH --job-name=ieqtl_pool
#SBATCH --partition=<gpu_partition>
#SBATCH --gres=gpu:a100:1
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96g
#SBATCH --time=04:00:00
#SBATCH --output=logs/ieqtl_pool_%j.log
# Pooled-ages interaction eQTL on GPU, both covariate sets (sensitivity).
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PY="$(conda info --base)/envs/Env1/bin/python"
S=scripts/08_interaction_eqtl_pooled.py
for tis in saccule utricle; do
  for cov in core full; do
    echo ">>> $tis $cov"; $PY $S --tissue "$tis" --covar "$cov" --celltypes all
  done
done
echo ">>> ALL DONE"
