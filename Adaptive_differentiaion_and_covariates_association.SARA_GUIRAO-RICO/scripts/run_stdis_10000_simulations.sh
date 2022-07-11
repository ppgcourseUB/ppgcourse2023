#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=bp_stdis_10000_sim                                                                                         
#SBATCH --error bp_stdis_10000-%j.err                                                                                     
#SBATCH --output bp_stdis_10000-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8 

# directories
INPUT=../input/hgdp.geno
cd $INPUT

# module load                                                                                                           
module load BayPass   

# run BayPass (CORE Model) with the 10000 PODs as input
./g_baypass -npop 52 -gfile G.hgdp_pods_10000 -efile covariates -scalecov -nthreads 8 -outprefix hgdp_stdis_10000_pods
 
