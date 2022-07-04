# MULTIPLE SEQUENCE ANALYSES 

**Instructor**: Ana Rojas

The goal of the following exercises is to get familiarised with Multiple Sequence Alignments (MSA) general procedures.

Any question email me at a.rojas.m@csic.es

***

## Software description and hints

==> to work in your personal computers [*local work*]: 
- Visualisers: Jalview (https://www.jalview.org/) or Belvu (https://www.sanger.ac.uk/tool/seqtools/)
- Taxonomy ids of the species to work: Ursus (taxid:9639) species and Ailuropod (taxid:9645) for BLAST searches 
- BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi, we will use ALWAYS BLAST PROTEIN.
- To retrieve the sequences: https://www.uniprot.org/
- To randomise and extract proteins: bbmap package (https://sourceforge.net/projects/bbmap/)

==> installed in the cluster: 
- To make non redundant datasets: cd-hit http://weizhong-lab.ucsd.edu/cd-hit/ or http://weizhong-lab.ucsd.edu/cdhit-web-server/cgi-bin/index.cgi?cmd=cd-hit
- Alignment method, mafft: https://mafft.cbrc.jp/alignment/software/source.html 
- Alignment method,clustal omega: https://www.ebi.ac.uk/Tools/msa/clustalo/
- Profile-based, hmmer https://www.ebi.ac.uk/Tools/hmmer/
***

## Data description and access

`MSA/BlockA/groups_to_align/*.ufasta`

`MSA/BlockB/queries/GTP_query2launch.txt`

`MSA/BlockB/queries/NR_query2launch.txt`

`MSA/BlockB/templates/GTP_template_aln.hmm`

`MSA/BlockB/templates/GTP_template_aln.ufasta`

`MSA/BlockB/templates/NR_template_aln.hmm`

`MSA/BlockB/templates/NR_template_aln.ufasta`

***

## Exercise  1: BlockA, ALIGNING bear sequences.

To get familiar with methods and visualisers.

0.- Select All groups from `MSA/BlockA/groups_to_align/*.ufasta`

1.- Run CD-HIT via command line (*you can prepare a job script to submit with `sbatch` like in the previous practical session*) or via web [*local work*] to generate non redundant sets (set -c 0.7). 

> **Remember that to work locally with these and other files during the course you must download the data FROM your personal computers first. This can be done by typing `scp -r user@host:/folder_to_download` in your terminal or using a FTP client like [FileZilla](https://filezilla-project.org/).**

2.- Write a small code to run an alignment method for each group. 

Here you have an example of job script to submit MAFFT works (*mafft.run*) with default parameters (using multiple threads):

```
#!/bin/bash

# define names
#SBATCH --job-name=mafft
#SBATCH --error mafft-%j.err
#SBATCH --output mafft-%j.out

# memory and CPUs request
#SBATCH --mem=6G
#SBATCH --cpus-per-task=8

# module load
module load mafft

# jobs to launch
(shopt -s failglob; for A in ./BlockA/groups_to_align/*.ufasta; do mafft --thread 8 ${A} > ${A}.mfft; done)

#post-processing
(shopt -s failglob; for i in ./BlockA/groups_to_align/*.mfft; do  mv $i ${i%.ufasta.mfft}.mfft.fasta; done)
```

To submit the job:

```
sbatch mafft.run
```

3.- Make a multiple sequence alignment of each of each alignment using clustal omega. 
Here you have an example job script to submit with `sbatch`:

```
#!/bin/bash

# define names
#SBATCH --job-name=clustalo
#SBATCH --error clustalo-%j.err
#SBATCH --output clustalo-%j.out

# memory and CPUs request
#SBATCH --mem=6G
#SBATCH --cpus-per-task=8

# module load
module load clustal-omega

# jobs to launch
(shopt -s failglob; for A in ./BlockA/groups_to_align/*.ufasta; do clustalo --threads 8 -i ${A} -o ${A}.clustalo  -v; done)

#post-processing
(shopt -s failglob; for i in ./BlockA/groups_to_align/*.clustalo; do  mv $i ${i%.ufasta.clustalo}.clustalo.fasta; done)
```

4.- Select 5 alignments of your choice (from most similar to more divergent), and randomise the order of your sequences within the group and repeat step number 3. 

To randomise files you can: `shuffle.sh in=file.fa out=shuffled.fa`  
To extract the 1% of files you can `reformat.sh in=file.fa out=sampled.fa samplerate=0.01`


5.- Visualise all the alignments using Jalview (https://www.jalview.org/) or Belvu (https://www.sanger.ac.uk/tool/seqtools/) [*local work*]. **Remember that to do that you must download the data from your personal computers first. In your local folder, type `scp -r user@host:/folder_to_download.`**

8.- Select any protein of the five alignments and IDENTIFY the proteins included in the alignment via BLAST or HMMER search (local work).

## Exercise 2: Block B, PROTEIN FAMILIES

In this Block you are going to analyse TWO big superfamilies of proteins (GTPasas and Nuclear Receptors). These families are essential signalling proteins. We want to analyse the presence of these proteins in our Bears. These proteins are usually multidomain, and also the subfamilies are distant, therefore posing a challenge for MSA.

0.- Create an individual workinf directory for each query file, i.e.: type `mkdir GTPases/` and `mkdir NRs/`

1.- Retrieve the fasta sequences using UNIPROT of [*local work*]: 

`MSA/BlockB/queries/GTP_query2launch.txt`

`MSA/BlockB/queries/NR_query2launch.txt`

2.- Run BLAST of protein of each query file restricting the taxonomy to Ursus (taxid:9639) species and Ailuropod (taxid:9645) [*local work*].

3.- select 10-20 diverse sequences (<80% seq, id) and save the FULL sequence to a file `bear_GTPases.ufasta` or `bear_NR.ufasta`. **Again, you can use either an `scp` command to upload the files to the cloud or an FTP client.**

4.- In the cloud, concatenate the file with the existing `MSA/BlockB/templates/GTP_template_aln.ufasta` or `MSA/BlockB/templates/NR_template_aln.ufasta`.

`cat  PATH/bear_GTPases.ufasta GTP_template_aln.ufasta > GTPall.ufasta`.

`cat  PATH/bear_NR.ufasta NR_template_aln.ufasta > NRall.ufasta`.

5.- run MAFFT or clustal with the BIG `GTPall.ufasta` file  (IMPORTANT HERE: those will be BIG alignments, select your parameters accordingly. TIP: you can remove the sequence redundancy using cd-hit http://weizhong-lab.ucsd.edu/cd-hit/.). **Create a new bash script to submit the job using sbatch`  accordingly**.

6.- Identify regions well-aligned, trim by hand regions not suitable for the alignment. Cut them visualy and save the short alignment. 

7.- Alternatively you could fit your sequences to the `hmmer` profiles available in the templates to check for regions which do not belong to our sequences of interest. This option can be only made via command line. Here you have a template script (**finish filling the running line with the correct file names and paths**) with: 

```
#!/bin/bash

# define names
#SBATCH --job-name=hmmer
#SBATCH --error hmmer-%j.err
#SBATCH --output hmmer-%j.out

# memory and CPUs request
#SBATCH --mem=6G
#SBATCH -c 8

# module load
module load hmmer

# jobs to launch
hmmalign -o OUTfile.fasta  <*.hmm> <YourBIG unaligned file>
```

8.-  Visualise all the alignments using Jalview (https://www.jalview.org/) or Belvu (https://www.sanger.ac.uk/tool/seqtools/), and check for blosks of unalignable regions [*local work*].


***

<br>**QUESTIONS (or other actions)**</br>
- Are the sequences of your alignments in both blocks single/multidomain?

- Which part of the sequence you have in your alignments?

- How diverse are your alignments overall on each block?

- which method works best for less similar sequences?

- Does randomising sequence order affect the final alignment?

- Why do I need to use a profile?

- Have I recovered all the Bear sequences belonging to this protein families?

