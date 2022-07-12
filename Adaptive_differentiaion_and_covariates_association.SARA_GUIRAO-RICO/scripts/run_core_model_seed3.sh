#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=bp_core_seed3                                                                                         
#SBATCH --error bp_core_seed3-%j.err                                                                                     
#SBATCH --output bp_core_seed3-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8 

# directories
INPUT=../input
cd $INPUT

# module load                                                                                                           
module load BayPass   

# run BayPass (CORE Model) with seed3
g_baypass -npop 52 -gfile hgdp.geno -nthreads 8 -seed 94875 -outprefix hgdp_core_s3
