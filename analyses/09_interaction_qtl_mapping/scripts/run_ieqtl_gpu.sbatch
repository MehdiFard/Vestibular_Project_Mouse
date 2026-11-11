#!/bin/bash
#SBATCH --job-name=ieqtl
#SBATCH --partition=<gpu_partition>
#SBATCH --gres=gpu:a100:1
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96g
#SBATCH --time=04:00:00
#SBATCH --output=logs/ieqtl_%j.log
# ---------------------------------------------------------------------------
# Cell-type interaction-eQTL mapping on GPU (tensorQTL), all groups x cell types.
# Interaction covariate = ssGSEA marker-enrichment (inverse-normal transformed).
#   sbatch scripts/run_ieqtl_gpu.sbatch
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PY="$(conda info --base)/envs/Env1/bin/python"
S=scripts/05_interaction_eqtl.py

for g in saccule_w8 saccule_w22 utricle_w8 utricle_w22; do
    echo ">>> $g"
    $PY $S --group "$g" --celltypes all
done
echo ">>> ALL DONE"
