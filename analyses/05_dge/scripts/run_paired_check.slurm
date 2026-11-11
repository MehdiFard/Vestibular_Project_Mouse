#!/bin/bash
#SBATCH --job-name=sac_utr_pu
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=<cpu_partition>
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
set -euo pipefail
BASE="$PWD/analyses/05_dge"
ENV="$(conda info --base)/envs/Env1"
cd "$BASE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "start $(date)"
"$ENV/bin/Rscript" scripts/15_paired_vs_unpaired_check.R
echo "end $(date)"
