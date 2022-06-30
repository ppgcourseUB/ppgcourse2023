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

You can find the subset of proteomes we are going to use in the folder data/proteomes/

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

`awk '{if ($0 ~ /^>/)  {split($0,a,"|"); split(a[3],b," "); print ">"b[1]} else { print;}}' proteomes/fileName >proteomes_parsed/fileName`

This small piece of code will parse your header and keep only one part of it which is unique for each protein. Make sure to name the output fileName with something short and descriptive of your species as OrthoFinder will be using it to identify each species.

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

1.- OrthoFinder has few parameters, but the most important one of them is the inflation parameter. This parameter indicates whether the orthogroups are going to be smaller or bigger. By default it is set to 1.5. We are now going to run orthoFinder with a bigger inflation parameter. There is no need to re-calculate the diamond search, in order to re-use previously calculated results, run the following command:

`orthofinder -b folderName/OrthoFinder/Results_XXXX/ -I 3.0 -og`

Note that we are using -b instead of -f and we are prividing previously calculated results. Also we are changing the inflation parameter using -I and setting it to 3.0. At this point we are only interested in comparing the orthogroups, the -og parameter will stop the run of orthoFinder after it calculates orthogroups. This is a time-saving trick if you want to assess different inflation parameters.

This will generate a second folder which will be called Results_XXX_1 where the new results of orthoFinder can be found.

2.- Focus on the folder called Orthogroups. In this folder you will find several files of interest:

* Orthogroups.tsv: Will print all the orthogroups detected during the analysis
* Orthogroups_UnassignedGenes.tsv: Genes that have not been assigned to an orthogroup will go to this file
* Orthogroups_SingleCopyOrthologues.txt: This will give you a list of orthogroups that did not have duplications and in which all species are present

3.- Now compare the results obtained by the two runs of orthoFinder and try and answer the following questions:

3.1.- Are the number of orthogroups the same in both runs? How many orthogroups are in each of them?

3.2.- Did the change in inflation parameter affect the detection of single copy genes?

3.3.- Note that in order to do a good phylogenomics analysis you need to find groups of orthologous genes that are present in all species without duplications. Which file would give you this information and how many orthogroups can you use in each run?

3.4.- Nowadays many reviewers are asking how the inflation parameter can affect your results. Can you think on a way to show which is the correct inflation parameter?

## Exercise 4

Orthogroups can contain duplications which means we can have a mix of orthologs and paralogs. OrthoFinder implements a method to distinguish between them, but to do that it needs a species tree as reference. OrthoFinder tries to calculate the species tree on its own. We will go back to using the first run of orthoFinder. In that folder you should see a folder called Species_tree and in there is a file called SpeciesTree_rooted.txt

1.- Check the species tree that has been automatically build by orthofinder.

http://phylo.io/

Based on your knowledge on how bears have evolved, do you see anything wrong with this species tree?

2.- If the tree is not rooted correctly, re-root the tree by clicking on the branch and pressing on re-root (make sure you press on the branch and not on the species name!). Now export the newick (find the button on the upper right part of the image, press on export nwk). Save the file into your folder.

3.- As before we are not going to re-run the whole pipeline. This time we will start from the pre-calculated orthogroups and will just change the species tree:

`orthofinder -fg FolderRun -s speciesTree_file`

* -fg indicates you want to re run an analysis from the orthogroups on.

* -s indicates you will provide a user defined species tree

Running this will generate a new folder, which will contain the new results. 

Note that if you have a species tree before running orthoFinder it is more convenient to just provide it from the start using the -s option. Just make sure that the name of the proteome files are the same as the ones found in the species tree. 

Additionally it is possible that in this toy example with very few species the change in species tree will have little effect on the prediction of orthologs, but this can change in more complicated scenarios)

## Exercise 4

Until now orthofinder has been using distance matrices to build the gene trees, but, while fast, this is not the
most accurate way to build them.

1.- In order to change the way trees are build we can use the following command:

`orthofinder -fg FolderRun -s speciesTree_file -t 2 -M msa`

-M indicates you want to build trees using mafft as a multiple sequence aligner and fasttree as a tree
reconstruction program. Due to time constrains we will not run orthoFinder with iqtree, though it is advisable to do so if you have enough computational power. You would need to add -T iqtree to the command line.

## Exercise 5

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



