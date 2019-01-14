#!/bin/bash
## Project:
#SBATCH --account=NN9526K
#SBATCH --partition=bigmem
#SBATCH --job-name=cosmic_t
#SBATCH --time=2-00:00:00
#SBATCH --output=/cluster/work/users/dsarria/slurm_logs/R-%x.%j.out
#SBATCH --error=/cluster/work/users/dsarria/slurm_logs/R-%x.%j.err
#SBATCH --nodes=1 --ntasks-per-node=32 --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## Software modules
module restore system   # Restore loaded modules to the default

module load foss/2020b
module load CMake/3.9.1
module load OpenMPI/2.1.1-GCC-6.4.0-2.28
module load gompi/2017b
module load Python/3.6.2-foss-2017b
module load Boost/1.66.0-foss-2018a-Python-2.7.14
module load GCC/10.2.0
module load GCCcore/10.2.0

module list             # List loaded modules, for easier debugging

## Run the application
echo "Running code..."
cd ${SLURM_SUBMIT_DIR}

srun python python_job_fram_BM.py

wait
echo "Done."
