#!/bin/bash
#
#SBATCH --job-name=2runRplots
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --ntasks=2
#SBATCH --mem=6GB
#SBATCH --partition=normal
#
module load r-mass r-proto

srun --ntasks 1 --exclusive --mem-per-cpu=1GB R --vanilla < ./Results_plotMKT.R&
srun --ntasks 1 --exclusive --mem-per-cpu=1GB R --vanilla < ./Results_plotMKT_Theta.R&
wait
