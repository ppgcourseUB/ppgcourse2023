# ORTHOLOGY PREDICTION FOR PHYLOGENOMIC ANALYSES 

The objective of this practice is to learn how to use [Orthofinder](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y) to infer orthogroups that can then be used to do a phylogenomics analysis.

***

## Software desrciption and hints

For this practice we will use [Orthofinder](https://github.com/davidemms/OrthoFinder). You can find a full description on the program in the provided link. It is important to understand that OrthoFinder is a pipeline that uses different programs and that some of the steps can be performed by different programs. For instance, the first step is an all against all homology search that can be done using either BlastP or Diamond.

The input of the program will be a folder containing one file per proteome. These files need to contain fasta formated sequences of proteins and there needs to be a single file per species included in the analysis. Files need to have .fa, .faa, .fasta or similar for their termination, else orthofinder will not recognise them.

***

## Data description and access

We will work with a subset of four bear proteomes. The original proteomes can be downloaded from UniProt:

* Ursus arctos (brown bear): https://www.uniprot.org/UP000286642
* Ursus maritimus (polar bear): https://www.uniprot.org/UP000261680
* Ursus americanus (black bear): https://www.uniprot.org/UP000291022
* Ailuropoda melanoleuca (giant panda, outgroup): https://www.uniprot.org/UP000008912

You can find the subset of proteomes we are going to use in the folder data/starting/

***

## Exercise (or example, etc..) 1

### 1.- Download the data into your Desktop folder. Once there uncompress the file by typing:

`tar -zxvf proteomes.tar.gz`

Have a look at the files:

`head proteomes/*`

### 2.- Headers in uniprot contain a lot of information, and while that can be useful, it can become a problem
farther down the pipeline. As such it is better to deal with it from the start. For each of the files execute the
following command:

2.1.- Create a new folder where we will put the parsed files:

`mkdir proteomes_parsed/`

2.2.- Now cut the headers for each proteome:

`cat proteomes/fileName | cut -f1 -d “ “ >proteomes_parsed/fileName`

Where you have to substitute fileName for each of the file names in the proteomes folder. The cat option
reads a file, the | symbol is used to concatenate two different commands so that one thing is done after the
other, usually it used the output of one command to do something. The cut option will cut every line in the
file by a given term, in this case it is a space given by the -d option. Finally the -f 1 option indicates you
only want to keep the first field.

For instance:

`cat proteomes/Amelanoleuca.fasta | cut -f1 -d ‘ ‘ >proteomes_parsed/Amelanoleuca.fasta`

2.3.- Compare now the headers:

`head -n1 proteomes/Amelanoleuca.fasta`

`head -n1 proteomes_parsed/Amelanoleuca.fasta`

Repeat that for each file and make sure you have four files in the proteomes_parsed folder.


## Exercise (or example, etc..) 2

Now, you will...

***

>**QUESTIONS (or other actions)**</br>
How ...?   
What is...?   
If you...?  

