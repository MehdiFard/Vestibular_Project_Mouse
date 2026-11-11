#!/bin/bash
### Differential expression for the saccule-vs-utricle analysis.
### The paired model carries ~70 strain coefficients, so this is far too heavy
### for a login node - the login-node reaper kills it after a few minutes.

#SBATCH --job-name=sac_utr_de
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=<cpu_partition>
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail

BASE="$PWD/analyses/05_dge"
ENV="$(conda info --base)/envs/Env1"

cd "$BASE"

# One BLAS thread per process; parallelism comes from BiocParallel over genes.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Differential expression starting on $(hostname) at $(date)"
"$ENV/bin/Rscript" scripts/03_differential_expression.R
echo "Differential expression finished at $(date)"

echo "Figures"
"$ENV/bin/Rscript" scripts/04_figures.R

echo "Report"
"$ENV/bin/python" scripts/05_build_report.py

echo "All done at $(date)"
