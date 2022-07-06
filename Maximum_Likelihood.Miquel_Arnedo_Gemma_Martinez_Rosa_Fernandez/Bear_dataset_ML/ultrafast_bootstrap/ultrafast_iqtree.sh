#!/bin/sh
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=10G
#SBATCH -o slurm.%j.%x.out
#SBATCH -e slurm.%j.%x.err
#SBATCH --time=0-5:0

#module load
module load iq-tree

iqtree -s 50_genes.fa -m MODEL -bb 1000 -nt 8
