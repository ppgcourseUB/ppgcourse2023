#!/bin/sh
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=10G
#SBATCH -o slurm.%j.%x.out
#SBATCH -e slurm.%j.%x.err
#SBATCH --time=0-5:0

module load cesga/2020
module load iq-tree/2.1.2

iqtree2 -s 50_genes.fa -m MODEL -alrt 1000 -T AUTO
