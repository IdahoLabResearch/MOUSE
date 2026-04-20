#!/bin/bash
#SBATCH --job-name=mpi_example
#SBATCH --partition=short
#SBATCH --time=5:58:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=25
#SBATCH --wckey=edu_res
#SBATCH --error=serial.err.%J
#SBATCH --output=serial.txt.%J

source ~/.bashrc
conda activate snap

cd examples/
mpirun -np 2 python -u watts_exec_LTMR_mpi.py