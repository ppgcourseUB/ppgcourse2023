#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=bp_stdis_contr                                                                                         
#SBATCH --error bp_stdis_contr-%j.err                                                                                     
#SBATCH --output bp_stdis_contr-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8 

# directories
INPUT=../input/hgdp.geno
cd $INPUT

# module load                                                                                                           
module load BayPass   

# run BayPass (STDis and Contrast Model)
./g_baypass -npop 52 -gfile hgdp.geno -contrastfile covariates_eu -efile covariates_eu -nthreads 8 -d0yij 20 -outprefix hgdp_contrast
