#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=bp_core_1000_sim                                                                                         
#SBATCH --error bp_core_1000_sim-%j.err                                                                                     
#SBATCH --output bp_core_1000_sim-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8 

# directories
INPUT=../input
cd $INPUT

# module load                                                                                                           
module load BayPass   

# run BayPass (CORE Model) with the 1000 PODs as input
g_baypass -npop 52 -gfile G.hgdp_pods_1000 -nthreads 8 -outprefix hgdp_pod_1000 
