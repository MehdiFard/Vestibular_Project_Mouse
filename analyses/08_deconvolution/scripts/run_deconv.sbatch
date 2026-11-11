#!/bin/bash
#SBATCH --job-name=deconv
#SBATCH --partition=<cpu_partition>
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128g
#SBATCH --time=03:00:00
#SBATCH --output=logs/deconv_pipeline_%j.log
# ---------------------------------------------------------------------------
# End-to-end deconvolution pipeline (run from the project root):
#   sbatch scripts/run_deconv.sbatch
# Each step uses the env that has its packages:
#   01 bulk TPM      -> Env1      (python/pandas)
#   02 sc reference  -> seurat_env(Seurat)
#   03 deconvolution -> deconv_env(Bisque/MuSiC/GSVA)
#   03b validation   -> deconv_env
#   04 compare/plots -> deconv_env
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PY_ENV1="$(conda info --base)/envs/Env1/bin/python"
R_SEURAT="$(conda info --base)/envs/seurat_env/bin/Rscript"
R_DECONV="$(conda info --base)/envs/deconv_env/bin/Rscript"
D=scripts

echo ">>> [01] bulk TPM matrices";       $PY_ENV1 $D/01_build_bulk_tpm.py
echo ">>> [02] single-cell reference";   $R_SEURAT $D/02_build_sc_reference.R
echo ">>> [03] deconvolution (3 methods)"; $R_DECONV $D/03_run_deconv.R
echo ">>> [03b] pseudobulk validation";  $R_DECONV $D/03b_pseudobulk_check.R
echo ">>> [04] compare + plots";         $R_DECONV $D/04_validate_compare.R
echo ">>> DONE"
