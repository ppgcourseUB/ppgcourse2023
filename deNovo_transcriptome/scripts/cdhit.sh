#!/bin/bash

# define names
#SBATCH --job-name=cdhit
#SBATCH --error cdhit-%j.err
#SBATCH --output cdhit-%j.out

# memory and CPUs request
#SBATCH --mem=50G
#SBATCH --cpus-per-task=5

# module load
module load cdhit

cd-hit-est -i /xxx/session2/4trinity/Trinity.fasta -o /xxx/session2/5cdhit/Trinity_cdhit.fasta -c 0.9 -M 50 -T 5 > /xxx/session2/5cdhit/Trinity_cdhit.err
