#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=bp_1000_sim                                                                                         
#SBATCH --error bp_1000_sim-%j.err                                                                                     
#SBATCH --output bp_1000_sim-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8 

# directories
INPUT=../input/hgdp.geno
cd $INPUT

# module load                                                                                                           
module load BayPass   

# get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution obtained when running the CORE Model
pi.beta.coef=read.table("hgdp_core_s1_summary_beta_params.out",h=T)$Mean

# upload the original data to obtain total allele count (sample size for each population). Do this by using the geno2YN() function in baypass_utils.R script
hgdp.data<-geno2YN("hgdp.geno")

# read the omega matrix from seed1 obtained when running the CORE Model:
omega_s1=as.matrix(read.table(file="hgdp_core_s1_mat_omega.out", header=F))

# simulated 1000 PODs
simu.hgdp_1000 <- simulate.baypass(omega.mat=omega_s1, nsnp=1000, 
    sample.size=hgdp.data$NN, beta.pi=pi.beta.coef, pi.maf=0, suffix="hgdp_pods_1000")
