#!/bin/bash

# Submit this script with: qsub thefilename

#PBS -l select=1:ncpus=16:mem=12gb
#PBS -l walltime=00:20:00
#PBS -j oe
#PBS -N openmc_eco
#PBS -P edu_class
#PBS -m be
#PBS -M botros.hanna@inl.gov

cd /home/hannbn/projects/MARVEL_MRP/Github_repos/watts/examples/openmc_economics

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source activate /home/hannbn/.conda/envs/openmc-env/envs/watts_openmc
# Your job commands go here
python watts_exec.py > output.txt