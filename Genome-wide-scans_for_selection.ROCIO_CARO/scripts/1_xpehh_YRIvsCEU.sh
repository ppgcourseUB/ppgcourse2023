#!/bin/bash     

# define names
#SBATCH --job-name=selscan
#SBATCH --error selscan-%j.err
#SBATCH --output selsan-%j.out 

# memory and CPUs request
#SBATCH --mem=6G
#SBATCH --cpus-per-task=8
module load selscan/v1.2.0a

HAP1=data/YRI/YRI.recode.vcf.hap
HAP2=data/CEU/CEU.recode.vcf.hap
MAP=data/CEU/CEU_filled.map
XPEHH_FILE=YRIvsCEU

NORM_XPEHH_FILE=YRIvsCEU.xpehh.out

# XPEHH computation
selscan --xpehh --hap $HAP1 --ref $HAP2 --map $MAP --out $XPEHH_FILE --threads 8

# Normalization
norm --xpehh --files $NORM_XPEHH_FILE

