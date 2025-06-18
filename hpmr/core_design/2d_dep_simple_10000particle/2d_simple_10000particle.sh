#!/bin/bash

# Submit this script with: sbatch thefilename

#SBATCH --time=5:00:00 # walltime
#SBATCH --ntasks-per-node=56 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1 # number of nodes
#SBATCH --wckey edu_class # Project Code
#SBATCH -J "2d_simple_10000p" # job name
#SBATCH --mail-user=reham.abdelnasser@inl.gov # email address

#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#export OPENMC_CROSS_SECTIONS='/home/abderi/cross_sections/endfb-viii.0-hdf5/cross_sections.xml'
conda activate openmc-env
python 2d_simple_10000particle.py