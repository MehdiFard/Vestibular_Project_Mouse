#!/bin/bash
#SBATCH --job-name=ieqtl_core
#SBATCH --partition=<gpu_partition>
#SBATCH --gres=gpu:a100:1
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96g
#SBATCH --time=04:00:00
#SBATCH --output=logs/ieqtl_core_%j.log
# Per-group interaction eQTL with the CORE covariate model (geno PCs only, no RNA PCs).
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PY="$(conda info --base)/envs/Env1/bin/python"
S=scripts/05_interaction_eqtl.py
for g in saccule_w8 saccule_w22 utricle_w8 utricle_w22; do
  echo ">>> $g (core)"
  $PY $S --group "$g" --covar core --outroot results/ieqtl_core --celltypes all
done
echo ">>> ALL DONE"
