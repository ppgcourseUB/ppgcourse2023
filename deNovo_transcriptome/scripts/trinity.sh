#!/bin/bash

# define names
#SBATCH --job-name=trinity
#SBATCH --error trinity-%j.err
#SBATCH --output trinity-%j.out

# memory and CPUs request
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10

# module load
module load trinity

# jobs to launch
srun Trinity --left /xxx/session2/2trimmed_data/reads_1.P.fq.gz \
	--right /xxx/session2/2trimmed_data/reads_2.P.fq.gz --seqType fq \
	--normalize_reads --normalize_max_read_cov 30 \
	--max_memory 100G \
	--CPU 10 \
	--output /xxx/session2/4trinity
