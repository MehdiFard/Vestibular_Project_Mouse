#!/bin/bash
#SBATCH --job-name=sac_utr_pca
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=<cpu_partition>
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2:00:00

set -euo pipefail
BASE="$PWD/analyses/05_dge"
ENV="$(conda info --base)/envs/Env1"
cd "$BASE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "PCA starting on $(hostname) at $(date)"
"$ENV/bin/Rscript" scripts/02_pca.R
echo "PCA done at $(date)"
