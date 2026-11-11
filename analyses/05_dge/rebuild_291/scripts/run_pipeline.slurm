#!/bin/bash
### Everything downstream of quantification, for the 291-library rebuild.
#SBATCH --job-name=v2_pipeline
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=<cpu_partition>
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
set -euo pipefail
RB="$PWD/analyses/05_dge/rebuild_291"
ENV1="$(conda info --base)/envs/Env1"
ENRICH="$(conda info --base)/envs/enrich_r"
cd "$RB"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="$RB/logs"

echo "=== 01 inputs $(date) ==="   && "$ENV1/bin/Rscript" scripts/01_build_inputs.R
echo "=== 02 pca $(date) ==="      && "$ENV1/bin/Rscript" scripts/02_pca.R
echo "=== 03 DE $(date) ==="       && "$ENV1/bin/Rscript" scripts/03_differential_expression.R
echo "=== 04 figures $(date) ===" && "$ENV1/bin/Rscript" scripts/04_figures_and_compare.R
echo "=== 05 gsea $(date) ==="     && "$ENRICH/bin/Rscript" scripts/05_gsea.R
echo "=== 06 report $(date) ==="   && "$ENV1/bin/python"  scripts/06_build_report.py
echo "=== all done $(date) ==="
