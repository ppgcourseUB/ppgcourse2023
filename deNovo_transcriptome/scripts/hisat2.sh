#!/bin/bash

# define names
#SBATCH --job-name=hisat2
#SBATCH --error hisat2-%j.err
#SBATCH --output hisat2-%j.out

# memory and CPUs request
#SBATCH --mem=50G
#SBATCH --cpus-per-task=5

# module load
module load hisat/2.2.1

# jobs to launch
srun hisat2-build /xxx/session2/4trinity/Trinity.fasta /xxx/session2/4trinity/Trinity

srun hisat2 -p 10 -x /xxx/session2/4trinity/Trinity \
	-1 /xxx/session2/2trimmed_data/reads_1.P.fq.gz -2 /xxx/session2/2trimmed_data/reads_2.P.fq.gz \
	-S /xxx/session2/6QC_trinity/reads.sam &> /xxx/session2/6QC_trinity/reads.sam.info
