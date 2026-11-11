#!/bin/bash
### Pantry - phenotyping - wrapper

#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=gpu-core-sponsored
#SBATCH --account=gpu-core-sponsored
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --gres=gpu:1
#SBATCH --mem=256G
#SBATCH --time=24:00:00  
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<email>

# Activating env to use apptainer
MAMBA="$(conda info --base)/etc/profile.d/mamba.sh"
source $MAMBA
mamba activate apptainer_env

# Activating docker container and running the scripts
CONTAINER="/path/to/posit-base-20250410-ppm-testing.sif"
apptainer exec --nv --bind /path/to/data $CONTAINER run_phenotyping_gpu.sh $DATASET

### END