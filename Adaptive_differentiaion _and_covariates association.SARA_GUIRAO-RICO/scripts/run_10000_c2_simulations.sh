#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=bp_10000_c2_sim                                                                                         
#SBATCH --error bp_10000_c2_sim-%j.err                                                                                     
#SBATCH --output bp_10000_c2_sim-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8 

# directories
INPUT=../input/hgdp.geno
cd $INPUT

# module load                                                                                                           
module load BayPass   

# get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
c2.pi.beta.coef=read.table("hgdp_contrast_summary_beta_params.out",h=T)$Mean

# upload the original data to obtain total allele count
hgdp.data<-geno2YN("hgdp.geno")

# read the omega matrix:
omega_contrast=as.matrix(read.table(file="hgdp_contrast_mat_omega.out", header=F))

# generate 10000 PODs
simu.C2.10000 <- simulate.baypass(omega.mat=omega_contrast,nsnp=10000,
	sample.size=hgdp.data$NN, beta.pi=c2.pi.beta.coef, pi.maf=0, 
	suffix="hgdp_C2_10000_pods")
