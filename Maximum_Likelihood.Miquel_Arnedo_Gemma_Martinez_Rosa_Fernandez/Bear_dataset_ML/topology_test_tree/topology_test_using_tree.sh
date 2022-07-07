#!/bin/sh
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -o slurm.%j.%x.out
#SBATCH -e slurm.%j.%x.err
#SBATCH --time=0-1:0

#module load
module load iq-tree

iqtree -s 50_genes.fa -z ../bear_species_trees_topologies.tre -te ../ultrafast_bootstrap/50_genes.fa.treefile -redo
