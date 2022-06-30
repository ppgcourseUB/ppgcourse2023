#!/bin/bash

# define names
#SBATCH --job-name=trimmomatic
#SBATCH --error trimmomatic-%j.err
#SBATCH --output trimmomatic-%j.out

# memory and CPUs request
#SBATCH --mem=10G
#SBATCH --cpus-per-task=5

# module load
module load java
module load trimmomatic

# jobs to launch
srun java -jar /cluster/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 5 \
	/xxx/session2/0data/reads.left.fq.gz /xxx/session2/0data/reads.right.fq.gz \
	/xxx/session2/2trimmed_data/reads_1.P.fq.gz /xxx/session2/2trimmed_data/reads_1.U.fq.gz \
	/xxx/session2/2trimmed_data/reads_2.P.fq.gz /xxx/session2/2trimmed_data/reads_2.U.fq.gz \
	ILLUMINACLIP:/xxx/session2/TruSeq3-PE.fa:2:30:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

