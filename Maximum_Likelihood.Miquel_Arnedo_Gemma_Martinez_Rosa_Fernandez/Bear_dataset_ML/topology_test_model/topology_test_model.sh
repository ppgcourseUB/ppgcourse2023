#!/bin/sh
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -o slurm.%j.%x.out
#SBATCH -e slurm.%j.%x.err
#SBATCH --time=0-2:0

#module load
module load iq-tree

iqtree -s 50_genes.fa -z ../bear_species_trees_topologies.tre -m LG+C20+F+G -redo
