#!/bin/bash

# Submit this script with: sbatch thefilename

#SBATCH --time=40:00:00 # walltime
#SBATCH --ntasks-per-node=56 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1 # number of nodes
#SBATCH --wckey edu_class # Project Code
#SBATCH -J "2d_dep_700k" # job name
#SBATCH --mail-user=reham.abdelnasser@inl.gov # email address

#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
conda activate openmc-env
python 2d_dep_600k.py