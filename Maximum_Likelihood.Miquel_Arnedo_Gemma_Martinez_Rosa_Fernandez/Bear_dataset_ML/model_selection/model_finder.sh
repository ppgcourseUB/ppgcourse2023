#!/bin/sh
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=10G
#SBATCH -o slurm.%j.%x.out
#SBATCH -e slurm.%j.%x.err
#SBATCH --time=0-5:0

module load cesga/2020
module load iq-tree/2.1.2

iqtree2 -s 50_genes.fa -m MF -mset LG+F+G,WAG+F+G,JTT+F+G,GTR20 -madd LG+C20+F+G,LG+C10+F+G,LG+C30+F+G,LG+C40+F+G,LG+C50+F+G,LG+C60+F+G,C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,LG4M,LG4X -T AUTO
