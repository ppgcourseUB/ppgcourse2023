#!/bin/bash

# define names
#SBATCH --job-name=FastQC
#SBATCH --error FastQC-%j.err
#SBATCH --output FastQC-%j.out

# memory and CPUs request
#SBATCH --mem=10G
#SBATCH --cpus-per-task=5

# module load
module load fastqc/0.11.9

# jobs to launch
srun fastqc -t 5 /xxx/session2/0data/reads.left.fq.gz -o /xxx/session2/1QC
srun fastqc -t 5 /xxx/session2/0data/reads.right.fq.gz -o /xxx/session2/1QC

