#!/usr/bin/env bash

#SBATCH -n 12               # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-06:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared           # serial_requeue   # Partition to submit to
#SBATCH --mem=50000 # MB    # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o temp/log/cgd_fit_%j.log  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e temp/log/cgd_fit_%j.log  # File to which STDERR will be written, %j inserts jobid

module load gcc R # gcc/9.2.0-fasrc01 R_core/3.6.3-fasrc01

Rscript cgd_fit.R --num-runs=10 --output=cgd_optim2
#Rscript cgd_fit.R --append-stage-2=data/cgd_optim.rds

