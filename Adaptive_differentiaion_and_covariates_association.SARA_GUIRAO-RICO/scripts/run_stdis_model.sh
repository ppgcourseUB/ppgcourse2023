#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=bp_stdis                                                                                         
#SBATCH --error bp_stdis-%j.err                                                                                     
#SBATCH --output bp_stdis-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8 

# directories
INPUT=../input
cd $INPUT

# module load                                                                                                           
module load BayPass   

# run BayPass (STDis Model)
g_baypass -npop 52 -gfile hgdp.geno -efile covariates -scalecov -nthreads 8 -outprefix hgdp_stdis
