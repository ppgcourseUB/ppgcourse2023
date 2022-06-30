#!/bin/bash

# define names
#SBATCH --job-name=trinityQC
#SBATCH --error trinityQC-%j.err
#SBATCH --output trinityQC-%j.out

# memory and CPUs request
#SBATCH --mem=10G
#SBATCH --cpus-per-task=5

# module load
module load trinity

# jobs to launch

srun /cluster/software/trinity/trinityrnaseq-v2.11.0/util/TrinityStats.pl /xxx/session2/4trinity/Trinity.fasta > /xxx/session2/6QC_trinity/Trinity_assembly.metrics
srun /cluster/software/trinity/trinityrnaseq-v2.11.0/util/TrinityStats.pl /xxx/session2/5cdhit/Trinity_cdhit.fasta > /xxx/session2/6QC_trinity/Trinity_cdhit_assembly.metrics
