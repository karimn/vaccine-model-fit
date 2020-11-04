#!/usr/bin/env bash

#SBATCH -n 24               # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-03:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared           # serial_requeue   # Partition to submit to
#SBATCH --mem=100000 # MB    # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o temp/log/cgd_fit_%j.log  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e temp/log/cgd_fit_%j.log  # File to which STDERR will be written, %j inserts jobid

module load gcc R_core 

#Rscript cgd_fit.R --num-runs=10 --parallel=$SLURM_NTASKS 
Rscript cgd_fit.R --num-runs=10 --vcov-moments-weight=0.02 --output=cgd_optim_vcov_weighted.rds --parallel=$SLURM_NTASKS  
# Rscript cgd_fit.R --num-runs=10 --output=cgd_optim_novcov --no-vcov-moments
#Rscript cgd_fit.R --append-stage-2=data/cgd_optim.rds

