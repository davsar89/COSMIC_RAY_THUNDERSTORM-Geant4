#!/bin/bash
## Project:
#SBATCH --account=NN9526K
#SBATCH --partition=bigmem
#SBATCH --job-name=MAKE_BDF
#SBATCH --time=0-8:0:0
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH -c 32
#SBATCH --mem-per-cpu=4G

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

module restore system
module load MATLAB/2019a
#module load Python/3.6.2-foss-2017b
#module load GCC/10.2.0
#module load GCCcore/10.2.0
#module load CMake/3.9.1

#exit 1 

echo "starting..."
cd ${SLURM_SUBMIT_DIR}
matlab -nodisplay -nodesktop -nojvm -r "global RECORD_PDG_TO_PROCESS; RECORD_PDG_TO_PROCESS = 22; MAKE_BIG_DATAFILES"
matlab -nodisplay -nodesktop -nojvm -r "global RECORD_PDG_TO_PROCESS; RECORD_PDG_TO_PROCESS = 11; MAKE_BIG_DATAFILES"
matlab -nodisplay -nodesktop -nojvm -r "global RECORD_PDG_TO_PROCESS; RECORD_PDG_TO_PROCESS = -11; MAKE_BIG_DATAFILES"

echo "done"
## Note: if you are using the Parallel Computing Toolbox, remove -nojvm
