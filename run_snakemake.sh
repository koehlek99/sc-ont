#!/bin/bash
#SBATCH --job-name=snakem_scl
#SBATCH --partition=compute     # the partition to submit into
#SBATCH --account=sc-users      # the SLURM account to charge
#SBATCH --nodes=1               # number of nodes to allocate
#SBATCH --ntasks=1              # number of processes (tasks) the job will start
#SBATCH --time=02-00:00         # how long the job is permitted to run, here 30 minutes


snakemake --slurm --default-resources --snakefile workflow/Snakefile -j 1

