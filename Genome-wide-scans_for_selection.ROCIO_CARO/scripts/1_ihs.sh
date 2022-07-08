#!/bin/bash     
                                                                                                 
# define names                                                                                                                                         
#SBATCH --job-name=selscan
#SBATCH --error selscan-%j.err
#SBATCH --output selsan-%j.out 

# memory and CPUs request
#SBATCH --mem=6G
#SBATCH --cpus-per-task=8


#### For CEU #### 
HAP_FILE=data/CEU/CEU.recode.vcf.hap
MAP_FILE=data/CEU/CEU_filled.map
IHS_FILE=CEU
NORM_IHS_FILE=CEU.ihs.out

# iHS computation
selscan --ihs --hap $HAP_FILE --map $MAP_FILE --out $IHS_FILE --threads 8

# Normalization
norm --ihs --files $NORM_IHS_FILE


#### For YRI ####
HAP_FILE=data/YRI/YRI.recode.vcf.hap
MAP_FILE=data/YRI/YRI_filled.map
IHS_FILE=YRI
NORM_IHS_FILE=YRI.ihs.out

# iHS computation
selscan --ihs --hap $HAP_FILE --map $MAP_FILE --out $IHS_FILE --threads 8

# Normalization
norm --ihs --files $NORM_IHS_FILE
