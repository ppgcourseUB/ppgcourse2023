# MULTIPLE SEQUENCE ANALYSES 

The goal of the following exercises is to get familiarised with Multiple Sequence Alignments general procedures.
Ana Rojas. Any question email me at a.rojas.m@csic.es

***

## Software description and hints

- Visualisers: Jalview (https://www.jalview.org/) or Belvu (https://www.sanger.ac.uk/tool/seqtools/)
- Taxonomy ids of the species to work: Ursus (taxid:9639) species and Ailuropod (taxid:9645) for BLAST searches
- BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi, we will use ALWAYS BLAST PROTEIN.
- To retrieve the sequences: https://www.uniprot.org/

== check if you have it istalled in the cluster: 
- To make non redundant datasets: cd-hit http://weizhong-lab.ucsd.edu/cd-hit/ or http://weizhong-lab.ucsd.edu/cdhit-web-server/cgi-bin/index.cgi?cmd=cd-hit
- Alignment method, mafft: https://mafft.cbrc.jp/alignment/software/source.html 
- Alignment method,clustal omega: https://www.ebi.ac.uk/Tools/msa/clustalo/
- Alignment method,t-coffee: https://tcoffee.crg.eu/
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

## Exercise  1, BlockA

To get familiar with methods and visualisers.

0.- Select All groups from `MSA/BlockA/groups_to_align/*.ufasta`

1.- Run CD-HIT via command line or web to generate non redundant sets (-c 0.7)

2.- Write a small code to run an alignment method for each group. 


A shell example will be this using MAFFT with default parameters:
`(shopt -s failglob; for A in *.ufasta; do mafft ${A} > ${A}.mfft; done)`

`(shopt -s failglob; for i in *.mfft; do  mv $i ${i%.ufasta.mfft}.mfft.fasta; done)`


3.- Make a multiple sequence alignment of each of each alignment using an alternative method of your choice: clustal, t-coffee.

4.- Select 5 alignments of your choice (from most similar to more divergent), and randomise the order of your sequences within the group and repeat step number 3. 

5.- Visualise all the alignments using Jalview (https://www.jalview.org/) or Belvu (https://www.sanger.ac.uk/tool/seqtools/).

8.- Select any protein of the five alignments and IDENTIFY the proteins included in the alignment via BLAST or HMMER.

## Exercise (or example, etc..) 2

In this Block you are going to analyse TWO big superfamilies of proteins (GTPasas and Nuclear Receptors). These families are essential signalling proteins. We want to analyse the presence of these proteins in our Bears. 
These proteins are usually multidomain, and also the subfamilies are distant, therefore posing a challenge for MSA.

0.- Create an individual workinf directory for each query file, i.e.: type `mkdir GTPases/` and `mkdir NRs/`

1.- Retrieve the fasta sequences using UNIPROT  of: 

`MSA/BlockB/queries/GTP_query2launch.txt`

`MSA/BlockB/queries/NR_query2launch.txt`

2.- Run BLAST of protein of each query file restricting the taxonomy to Ursus (taxid:9639) species and Ailuropod (taxid:9645).

3.- select 10-20 diverse sequences (<80% seq, id) and save the FULL sequence to a file `bear_GTPases.ufasta` or `bear_NR.ufasta`.

4.- Concatenate the file with the existing `MSA/BlockB/templates/GTP_template_aln.ufasta` or `MSA/BlockB/templates/NR_template_aln.ufasta`

`cat  PATH/bear_GTPases.ufasta GTP_template_aln.ufasta > GTPall.ufasta`.

`cat  PATH/bear_NR.ufasta NR_template_aln.ufasta > NRall.ufasta`.

5.- run MAFFT, Clustal, or T-Coffee with the BIG `GTPall.ufasta` file  (IMPORTANT HERE: those will be BIG alignments, select your parameters accordingly. TIP: you can remove the sequence redundancy using cd-hit http://weizhong-lab.ucsd.edu/cd-hit/.)

6.- Identify regions well-aligned, trim by hand regions not suitable for the alignment. Cut them visualy and save the short alignment. 

7.- Alternativey you could fit your sequences to the hmmer profiles available in the templates to check for regions which do not belong to our sequences of interest. This option can be only made via command line typing: 

`hmmalign -o OUTfile.fasta  <*.hmm> <YourBIG unaligned file>`.



***

<br>**QUESTIONS (or other actions)**</br>
- Are the sequences of your alignments in both blocks single/multidomain?
- 
- Which part of the sequence you have in your alignments?
- 
- How diverse are your alignments overall on each block?
- 
- which method works best for less similar sequences?
- 
- Does randomising sequence order affect the final alignment?
- 
- Why do I need to use a profile?

