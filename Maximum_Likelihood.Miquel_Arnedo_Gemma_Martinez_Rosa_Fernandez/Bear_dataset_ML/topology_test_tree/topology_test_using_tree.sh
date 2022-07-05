#!/bin/sh
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -o slurm.%j.%x.out
#SBATCH -e slurm.%j.%x.err
#SBATCH --time=0-1:0

module load cesga/2020
module load iq-tree/2.1.2

iqtree2 -s 50_genes.fa -z ../bear_species_trees_topologies.tre -te PREVIOUS_TREE
