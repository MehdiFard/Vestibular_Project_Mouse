#!/bin/bash
# Running Pantry - Pheast on GPU for one group (QTL mapping, heritability, TWAS weights) in the pantry_pheast env.
# Called by ../Snakefile (rule pheast), inside the Apptainer image when config pheast_container is set.
#   bash run_pheast_gpu.sh <group> [extra snakemake arguments]
# MAMBA_SH: miniforge mamba.sh (config mamba_sh); the number of jobs follows the CPUs given by Slurm.
DATASET=$1; shift
JOBS=${SLURM_CPUS_PER_TASK:-34}

# Activating conda env (in the container, if any)
MAMBA=${MAMBA_SH:?set MAMBA_SH (miniforge mamba.sh; config mamba_sh)}
source $MAMBA
mamba activate pantry_pheast

# Start time
echo "Starting workflow at: $(date)"

# Check
echo "==== GPU Info ===="
nvidia-smi
echo ""

echo "==== CUDA Compiler Version ===="
nvcc --version
echo ""

echo "==== Running Snakemake workflow ===="
echo ""

# Running snakemake workflow
snakemake   --configfile config.yml \
            --config dataset=$DATASET \
            --jobs $JOBS \
            --default-resources mem_mb=256000 cpus=$JOBS \
            --resources gpu=1 \
            --printshellcmds \
            --latency-wait 60 \
            --keep-going "$@"


# --keep-going
# -R tensorqtl_trans
# --rerun-incomplete


# Finish time
echo "Finished workflow at: $(date)"


### END
