#!/bin/bash

# define names
#SBATCH --job-name=transdecoder
#SBATCH --error transdecoder-%j.err
#SBATCH --output transdecoder-%j.out

# memory and CPUs request
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10

# module load
module load transdecoder

# jobs to launch

srun TransDecoder.LongOrfs -t /xxx/session2/4trinity/Trinity.fasta 
