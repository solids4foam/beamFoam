#!/bin/bash
#SBATCH --job-name=qub_mooredSemiSub
#SBATCH --output=log_%x_%j.out
#SBATCH --error=log_%x_%j.err
#SBATCH --partition=main
#SBATCH --ntasks=60
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00

# Good practice to navigate to the submission directory
cd $SLURM_SUBMIT_DIR

# Load OpenFOAM
#source /usr/lib/openfoam/openfoam2412/etc/bashrc
source /usr/lib/openfoam/openfoam2312/etc/bashrc
#source /usr/lib/openfoam/openfoam2212/etc/bashrc
#source /usr/lib/openfoam/openfoam2112/etc/bashrc

# Set number of threads for programs that use OpenMP parallelisation
# Note that OMP_NUM_THREADS defaults to 1 if not set
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=10

# Execute the run script
./Allrun_postProcess &> log.Allrun
#mpirun -np 60 interFoam -parallel &> log.interFoam
