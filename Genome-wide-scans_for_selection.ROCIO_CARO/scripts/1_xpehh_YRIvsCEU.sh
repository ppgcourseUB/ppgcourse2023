#!/bin/bash 

#SBATCH -J XPEHH_tutorial # A single job name for the array 

#SBATCH -n 1 # Number of cores (task per node)

#SBATCH -N 1 # All cores on one machine (nodes)

#SBATCH --cpus-per-task=8

#SBATCH -p normal # Partition 

#SBATCH --mem=6GB # Memory request  

#SBATCH -o XPEHH_tutorial.out # Standard output 

#SBATCH -e XPEHH_tutorial.err # Standard error  

#SBATCH --mail-type=ALL # mail alert at start, end and abortion of execution

#SBATCH --mail-user=rocio.caro@upf.edu # send mail to this address

module load selscan/v1.2.0a

HAP1=YRI.recode.vcf.hap
HAP2=CEU.recode.vcf.hap

MAP=CEU_filled.map

XPEHH_FILE=YRIvsCEU

NORM_XPEHH_FILE=YRIvsCEU.xpehh.out

# XPEHH computation
selscan --xpehh --hap $HAP1 --ref $HAP2 --map $MAP --out $XPEHH_FILE --threads 8

# Normalization
norm --xpehh --files $NORM_XPEHH_FILE

