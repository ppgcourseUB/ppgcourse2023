#!/bin/bash

# define names
#SBATCH --job-name=FastQC_loop
#SBATCH --error FastQC_loop-%j.err
#SBATCH --output FastQC_loop-%j.out

# memory and CPUs request
#SBATCH --mem=10G
#SBATCH --cpus-per-task=5

# module load
module load fastqc/0.11.9

# jobs to launch
for file in /xxx/session2/2trimmed_data/*fq.gz;
	do
		srun fastqc -t 5 $file -o /xxx/session2/3QC_trimmed
	done

