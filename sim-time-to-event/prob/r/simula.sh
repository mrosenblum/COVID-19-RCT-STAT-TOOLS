#!/bin/bash
#SBATCH --job-name=array_job
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #not required if only 1 thread as 1 thread per task is the default
#SBATCH --mem=4G
#SBATCH --partition=curie-cpu
#SBATCH --array=1-500
echo "$SLURM_ARRAY_TASK_ID"

source ~/.bashrc
spack load /mjrrusu
R CMD BATCH simula.r simula.rout
