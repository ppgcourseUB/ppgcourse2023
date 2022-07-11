#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=bp_core                                                                                         
#SBATCH --error bp_core-%j.err                                                                                     
#SBATCH --output bp_core-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8 

# directories
INPUT=../input/hgdp.geno
cd $INPUT

# module load                                                                                                           
module load BayPass   

# run BayPass (CORE Model) with different seeds
./g_baypass -npop 52 -gfile hgdp.geno -nthreads 8 -seed 15263 -outprefix hgdp_core_s1
./g_baypass -npop 52 -gfile hgdp.geno -nthreads 8 -seed 26847 -outprefix hgdp_core_s2
./g_baypass -npop 52 -gfile hgdp.geno -nthreads 8 -seed 94875 -outprefix hgdp_core_s3
