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

## Exercise 1

1.- Download the data into your Desktop folder. Once there uncompress the file by typing:

`tar -zxvf proteomes.tar.gz`

Have a look at the files:

`head proteomes/*`

2.- Headers in uniprot contain a lot of information, and while that can be useful, it can become a problem
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


## Exercise 2

We will run orthoFinder with the provided dataset and default values. To make sure it works quickly we
provided only a subset of the data.

1.- To execute orthofinder, type this command:

`orthofinder -f proteomes_parsed -t 2`

Where -f indicates the folder where the proteomes are found and -t indicates the number of threads
available to you.

2.- While the program works, lets have a look to the options that orthofinder has. Open a new terminal
and type orthofinder -h

## Exercise 3

We provided you with a dataset of four bear proteomes. Orthofinder automatically predicts which is the
correct species tree and uses it for the prediction of orthologs. The first step, always, is to compare it to
your current knowledge of the species tree of your species.

1.- Check the species tree that has been automatically build by orthofinder. You can find it in the folder
Species_tree within your orthofinder folder. You can visualize the tree using this website:

http://phylo.io/

Is there anything wrong with it?

2.- In the folder Species_tree, there is an additional folder called Potential_Rooted_Species_Trees/ where
you can find a list of different trees with different rooting points. Select the properly rooted species tree.
OrthoFinder allows you to re-run only parts of the analysis. In this case, we do not need to re-compute the
orthogroups therefore we can start running orthofinder from there:

`orthofinder -fg FolderRun -s speciesTree_file -t 2`

* -fg indicates you want to re run an analysis from the orthogroups on.

* -s indicates you will provide a user defined species tree

* -t indicates the number of threads

Running this will generate a second folder, which will contain the new results.

## Exercise 4

Until now orthofinder has been using distance matrices to build the gene trees, but, while fast, this is not the
most accurate way to build them.

1.- In order to change the way trees are build we can use the following command:

`orthofinder -fg FolderRun -s speciesTree_file -t 2 -M msa`

-M indicates you want to build trees using mafft as a multiple sequence aligner and fasttree as a tree
reconstruction program.

## Exercise 5

The default options for multiple sequence alignment reconstruction and species tree reconstruction can be
changed. The orthofinder manual explains how to do that.

We have run Orthofinder in three different ways: the initial run in which we were using default parameters
(initial run), a second run in which we corrected the rooting of the species tree and the final run where we
used a more accurate tree reconstruction method. Now it is time to look through the results. Try to answer
the following questions using the information found inside the different orthofinder folders.

1.- Lets study the evolutionary events that affect protein G1LY09_AILME.

1.1.- Search for the orthologs to this protein in the different runs. Have they been affected by the
differences in the analysis?

![image](https://user-images.githubusercontent.com/9434530/176634125-90b05063-77e1-4d3b-991f-a460978293e5.png)

1.2.- Repeat the same for the number of paralogs. Try to fill in this table:

![image](https://user-images.githubusercontent.com/9434530/176634263-6fb8a552-0c72-4fee-92a8-8a290e53a089.png)

1.3.- To which orthogroup does this protein belong to?

1.4.- phylo.io offers the possibility to see two trees side by side so that you can compare them. Go to the
compare window in phylo.io. Open on one side the Resolved Tree belonging to the orthogroup where your protein can be found for the first run and on the other side open the same tree for the second run. Can you
see the differences that caused the different predictions?

1.5.- Repeat the comparison between the tree in the second run and in the third run.

3.- We’re interested in studying the evolution of a given orthogroup (OG0000039). Scan the files to obtain
the following information:

3.1.- How many sequences are in each of the species within this orthogroup?

3.2.- Are there any single copy orthologs in there?

3.3.- How many duplication events were there in this orthogroup? Are you able to find where they
happened?

4.- Which two species have the highest number of one-to-one orthologs?

5.- One of the questions we had when beginning this study was whether the polar bear was more closely
related to the brown bear or to the american bear. Can you answer that question now?



