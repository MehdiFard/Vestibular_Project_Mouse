#!/bin/bash
# Run the whole saccule-vs-utricle analysis end to end.
#
#   bash scripts/run_all.sh
#
# Env1 is the only environment on this machine with DESeq2, tximport and
# python-pptx together, so everything runs from there.

set -euo pipefail

BASE="$PWD/analyses/05_dge"
ENV="$(conda info --base)/envs/Env1"

cd "$BASE"
mkdir -p logs

echo "[1/3] Building count matrix and sample table"
"$ENV/bin/Rscript" scripts/01_build_inputs.R 2>&1 | tee logs/01_build_inputs.log

echo "[2/3] PCA"
"$ENV/bin/Rscript" scripts/02_pca.R 2>&1 | tee logs/02_pca.log

# Steps 3-5 run on a compute node. The paired model carries ~70 strain
# coefficients and the login-node reaper kills it part way through, silently.
# On a compute node the whole thing takes about three minutes.
echo "[3/3] Submitting differential expression, figures and report to SLURM"
JOB=$(sbatch --parsable scripts/run_de.slurm)
echo "  job $JOB - follow it with: tail -f logs/sac_utr_de-$JOB.log"

echo
echo "When that job finishes: $BASE/report/saccule_vs_utricle_report.pptx"
echo
cat <<'NOTE'

NOTE: 05_build_report.py also reads tables written by the supplementary scripts,
so on a completely fresh run those have to exist before the deck can be built:

  06_batch_correction_check.R      (batch-correction demonstration)
  07_biotype_separation_check.R    biotype_separation_comparison.tsv
  08_limma_crosscheck.R            limma_vs_deseq2_summary.tsv
  09_vst_and_ntop_check.R          (VST / ntop figures)
  10_pairing_and_paired_age.R      paired_age_vs_batch_adjusted.tsv
  11_model_diagram.R               (model diagram figure)
  12_build_master_metadata.R       master_metadata.tsv
  13_infer_sex.R                   (sex validation figure)
  14_sex_effect.R                  sex_effect_summary.tsv
  15_paired_vs_unpaired_check.R    paired_vs_unpaired_summary.tsv   <- heavy, use
                                   sbatch scripts/run_paired_check.slurm
  16_interaction_illustration.R    interaction_examples.tsv

All but 15 run fine on a login node in a minute or two. Then rerun
  "$ENV/bin/python" scripts/05_build_report.py
NOTE
