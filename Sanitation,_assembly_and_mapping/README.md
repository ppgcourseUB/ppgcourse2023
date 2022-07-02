# From raw reads to processed alignments

Instructor: **Cinta Pegueroles Queralt**

The goal of this practical session is to learn how to obtain a <em>de novo</em> transcriptome from RNAseq data. The practical is organized in four main steps:

1. **Sanitation**: pre-processing RNAseq data  
2. **<em>de novo</em> assembly**: in the absence of a reference genome
3. **QC**: assessing the quality of de novo transcriptome
4. **Post-processing analyses**


### 0. Downloading data:
Due to time and computer resources we will work with a test data set from [TRINITY GITHUB][https://github.com/trinityrnaseq/trinityrnaseq/tree/devel/sample_data/test_Trinity_Assembly]

```
cd session2
mkdir 0data
cd 0data
wget https://github.com/trinityrnaseq/trinityrnaseq/raw/devel/sample_data/test_Trinity_Assembly/reads.left.fq.gz
wget https://github.com/trinityrnaseq/trinityrnaseq/raw/devel/sample_data/test_Trinity_Assembly/reads.right.fq.gz
```

### 1.1. Quality control of RNAseq data
Now we will perform the quality control (QC) of the RNAseq data using the software FastQC. First of all, we will create a folder for the QC output files

```
cd session2
mkdir 1QC
```
First we will create/edit a bash script to run FastQC software (fastqc.sh). IMPORTANT: mind to replace xxx by the current path.

```
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
```
To launch FastQC script in the cluster write on the terminal:
```
sbatch fastqc.sh
```
Use "s_jobs" command to check the state of the process. To check the results you need to download the files on a local computer. WARNING: mind to change the paths accordingly
```
scp -r user@cluster:/my_CLUSTER_path_to_session2/* my_LOCAL_path_to_session2/
```

Alternatively we could run the FastQC software directly on the terminal:
```
cd session2
fastqc -o 1QC 0data/reads.left.fq.gz 0data/reads.right.fq.gz
```

### 1.2. Trimming sequences

If we detect a drop in the quality of the sequences and/or the presence of adapters we need to trim our sequences. We will do it by using the Trimmomatic software.

First, we will create a folder for the trimmed sequences:

```
cd session2
mkdir 2trimmed_data
cd 2trimmed_data

```
Then we will create/edit a bash script to run Trimmomatic software (trimseq.sh). IMPORTANT: mind to replace xxx by the current path.

```
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
        /xxx/session2/0data/reads.left.fq.gz \ 
        /xxx/session2/0data/reads.right.fq.gz \
        /xxx/session2/2trimmed_data/reads_1.P.fq.gz \ 
        /xxx/session2/2trimmed_data/reads_1.U.fq.gz \
        /xxx/session2/2trimmed_data/reads_2.P.fq.gz \
        /xxx/session2/2trimmed_data/reads_2.U.fq.gz \
        ILLUMINACLIP:/xxx/session2/TruSeq3-PE.fa:2:30:1 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

To launch the Trimmomatic script in the cluster write on the terminal:
```
sbatch trimSeq.sh
```
Use "s_jobs" command to check the state of the process. 

Alternatively we could run the Trimmomatic software directly on the terminal:

```
trimmomatic PE ../0data/reads.left.fq.gz ../0data/reads.right.fq.gz \
    ../2trimmed_data/reads_1.P.fq.gz ../2trimmed_data/reads_1.U.fq.gz \
    ../2trimmed_data/reads_2.P.fq.gz ../2trimmed_data/reads_2.U.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 1.3. Quality control of the trimmed sequences
Now we will perform the quality control (QC) of the trimmed sequences using again the software FastQC. First of all, we will create a folder for the QC output files

```
cd session2
mkdir 3QC_trimmed
```
Since now we have 4 files to check (2 paired, 2 unpaired) we will use a loop to avoid repeating commands. We will create/edit a bash script to run FastQC software in loop over all trimmed sequences (fastqc_loop.sh). IMPORTANT: mind to replace xxx by the current path.

```
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
```

To launch fastQC_loop.sh script in the cluster write on the terminal:
```
sbatch fastqc_loop.sh
```
Use "s_jobs" command to check the state of the process. To check the results you need to download the files on a local computer. WARNING: mind to change the paths accordingly
```
scp -r user@cluster:/my_CLUSTER_path_to_session2/* my_LOCAL_path_to_session2/
```
Alternatively we could run the FastQC software directly on the terminal:
```
cd session2
for file in /xxx/session2/2trimmed_data/*fq.gz; do fastqc $file -o /xxx/session2/3QC_trimmed; done
```

### 2. <em>De novo</em> transcriptome assembly
### 2.1. Running trinity software

Now we will obtain a <em>de novo</em> transcriptome using the trimmed paired end sequences. In the same script we will also calculate some stats on the obtained transcriptome. First, we will create two folders to store the analysis:

```
cd session2
mkdir 4trinity #will contain the de novo transcriptome
```

First, we will create/edit a bash script to run trinity software (trinity.sh) to obtain a <em>de novo</em> transcriptome. IMPORTANT: mind to replace xxx by the current path.
```
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

```
To launch the trinity script in the cluster write on the terminal:
```
sbatch trinity.sh
```
Use "s_jobs" command to check the state of the process. 

Alternatively we could run the Trinity software directly on the terminal:
```
cd session2/4trinity

#launch trinity
Trinity --left /xxx/session2/2trimmed_data/reads_1.P.fq.gz \
	--right /xxx/session2/2trimmed_data/reads_2.P.fq.gz --seqType fq \
	--normalize_reads --normalize_max_read_cov 30 \
	--max_memory 100G \
	--CPU 10 \
	--output /cluster/home/cpegueroles/phylo/session2/4trinity
```
### 2.1. cd-hit: redundance reduction

Now we will avoid reduncamce in the de novo transcriptome by using the cdhit software. First we will create a folder to drop the analysis:

```
cd session2
mkdir 5cdhit
```

First we will create/edit a bash script to run cdhit software (cdhit.sh). IMPORTANT: mind to replace xxx by the current path.
```
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

```
To launch the script in the cluster write on the terminal:
```
sbatch cdhit.sh
```
Alternatively we could run the cdhit software directly on the terminal:
```
cd session2/5cdhit

cd-hit-est -i /xxx/session2/4trinity/Trinity.fasta -o /xxx/session2/5cdhit/Trinity_cdhit.fasta -c 0.9 -n 8 -M 0 -T 0 -d 40 > /xxx/session2/5cdhit/Trinity_cdhit.err


```
### 3. Post-processing

In this example we are analyszing the de novio transcriptoem obtained with trinity, but we could apply the same analysis to the reduced fasta file obtained after running cdhit software

### 3.1. Basic stats
We will check basic statistics of both trinity assembly and CDhit reduced transcriptome

```
cd session2
mkdir 6QC_trinity 
```

First we will create/edit a bash script to run the stats (trinityQC.sh) to obtain a de novo transcriptome. IMPORTANT: mind to replace xxx by the current path.

```
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
```
To launch the script in the cluster write on the terminal:
```
sbatch trinityQC.sh
```

### 3.2. Representation of reads
In this step we will backmap the RNAseq reads to the de novo assembly to check the completeness of the transcriptome. The higher the precentge of reads that backmap, the better.

To do so, we will use the paired trimmed reads (which are those that we used for the de novo assembly) and the software hisat2 for mapping.

First we will create/edit a bash script to run hisat2 software (hisat2.sh). IMPORTANT: mind to replace xxx by the current path.

```
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
srun hisat2-build /xxx/session2/trinity/Trinity.fasta /xxx/session2/4trinity/Trinity

srun hisat2 -p 10 -x /xxx/session2/4trinity/Trinity \
	-1 /xxx/session2/2trimmed_data/reads_1.P.fq.gz -2 /xxx/session2/2trimmed_data/reads_2.P.fq.gz \
	-S /xxx/session2/6QC_trinity/reads.sam &> /xxx/session2/6QC_trinity/reads.sam.info
```

### 3.3. From ncl to aa: TransDECODER

We will use transdecoder to get our predicted proteome. First, we will create a folder to drop the analysis, we will move to this folder and we will launch the analysis.

First we will create/edit a bash script to run transdecoder software (transdecoder.sh). IMPORTANT: mind to replace xxx by the current path.

```
cd session2
mkdir 6proteome
cd  6proteome
```
```
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
```

To launch the script in the cluster write on the terminal:
```
sbatch transdecoder.sh
```

#### Bibliography

Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-seq data without a reference genome. Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. PubMed PMID: 21572440.


Haas BJ, Papanicolaou A, Yassour M, Grabherr M, Blood PD, Bowden J, Couger MB, Eccles D, Li B, Lieber M, Macmanes MD, Ott M, Orvis J, Pochet N, Strozzi F, Weeks N, Westerman R, William T, Dewey CN, Henschel R, Leduc RD, Friedman N, Regev A. De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nat Protoc. 2013 Aug;8(8):1494-512. Open Access in PMC doi: 10.1038/nprot.2013.084. Epub 2013 Jul 11. PubMed PMID:23845962.
