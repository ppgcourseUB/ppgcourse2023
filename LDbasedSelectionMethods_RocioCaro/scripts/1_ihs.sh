#!/bin/bash 

#SBATCH -J iHS_tutorial # A single job name for the array 

#SBATCH -n 1 # Number of cores (task per node)

#SBATCH -N 1 # All cores on one machine (nodes)

#SBATCH --cpus-per-task=8

#SBATCH -p normal # Partition 

#SBATCH --mem=6GB # Memory request  

#SBATCH -o iHS_tutorial.out # Standard output 

#SBATCH -e iHS_tutorial.err # Standard error  

#SBATCH --mail-type=ALL # mail alert at start, end and abortion of execution

#SBATCH --mail-user=YOURMAIL # send mail to this address

module load selscan/v1.2.0a


#### For CEU ####

HAP_FILE=CEU.recode.vcf.hap

MAP_FILE=CEU_filled.map

IHS_FILE=CEU

NORM_IHS_FILE=CEU.ihs.out

# iHS computation
selscan --ihs --hap $HAP_FILE --map $MAP_FILE --out $IHS_FILE --threads 8

# Normalization
norm --ihs --files $NORM_IHS_FILE


#### For YRI ####

HAP_FILE=YRI.recode.vcf.hap

MAP_FILE=YRI_filled.map

IHS_FILE=YRI

NORM_IHS_FILE=YRI.ihs.out

# iHS computation
selscan --ihs --hap $HAP_FILE --map $MAP_FILE --out $IHS_FILE --threads 8

# Normalization
norm --ihs --files $NORM_IHS_FILE