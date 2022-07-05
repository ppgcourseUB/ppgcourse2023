# ORTHOLOGY PREDICTION FOR PHYLOGENOMIC ANALYSES 

**Instructor**: Marina Marcet-Houben

The objective of this practice is to learn how to use [OrthoFinder](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y) to infer orthogroups that can then be used to do a phylogenomics analysis.

***

## Software desrciption and hints

For this practice we will use [OrthoFinder](https://github.com/davidemms/OrthoFinder). You can find a full description on the program in the provided link. It is important to understand that OrthoFinder is a pipeline that uses different programs and that some of the steps can be performed by different programs. For instance, the first step is an all against all homology search that can be done using either BlastP or Diamond.

The input of the program will be a folder containing one file per proteome. These files need to contain fasta formated sequences of proteins and there needs to be a single file per species included in the analysis. Files need to have .fa, .faa, .fasta or similar for their termination, else `OrthoFinder` will not recognise them.

***

## Data description and access

We will work with a subset of four bear proteomes. The original proteomes can be downloaded from UniProt:

* Ursus arctos (brown bear): https://www.uniprot.org/proteomes/UP000286642
* Ursus maritimus (polar bear): https://www.uniprot.org/proteomes/UP000261680
* Ursus americanus (black bear): https://www.uniprot.org/proteomes/UP000291022
* Ailuropoda melanoleuca (giant panda, outgroup): https://www.uniprot.org/proteomes/UP000008912

You can find the subset of proteomes we are going to use in the folder proteomes/

***

## Exercise 1

We will run `OrthoFinder` with the provided dataset and default values. To make sure it works quickly we
provided only a subset of the data.

1.- To execute `OrthoFinder`, submit the job script *orthofinder.run* using the command `sbatch`.
Here you can see the contents of the script *orthofinder.run*:

```bash
#!/bin/bash                                                                                                             

##This is a script to run orthofinder                                                        

#SBATCH -p normal                                                                                      
#SBATCH -c 8                                                                                                     
#SBATCH --mem=6GB                                                                                       
#SBATCH --job-name orthofinder-job01                                                                            
#SBATCH -o %j.out                                                          
#SBATCH -e %j.err                                                      

#module loadding. Check available modules with `module avail` 
module load orthofinder

#running orthofinder
orthofinder -f proteomes -og -t 8 -a 2
```

> The flag -f indicates the folder where the proteomes are found. The -og option will tell `OrthoFinder` to stop after inferring the orthogroups. You can adapt the program for multi-threading using the -t and -a options.

2.- While the program works, lets have a look to the options that `OrthoFinder` has. Type `orthofinder -h`

> Note that one of the options is `-S` which determines how the homology search will be done. `Diamond` is very fast, much faster than `Blast`, but is less sensitive when running with distantly related species. Consider this when running `OrthoFinder`.

3.- Once `OrthoFinder` has finished, you will find that in the folder proteomes there will be a new folder called OrthoFinder and within it will be a file called Results_Jul06. The results of `OrthoFinder` can be found there sorted in different folders. New runs of `OrthoFinder` will generate new folders called Results_Jul06_1, Results_Jul06_2, ... 

Go to the folder and focus on the folder called *Orthogroups*. In this folder you will find several files of interest:

* *Orthogroups.tsv*: Will print all the orthogroups detected during the analysis
* *Orthogroups_UnassignedGenes.tsv*: Genes that have not been assigned to an orthogroup will go to this file
* *Orthogroups_SingleCopyOrthologues.txt*: This will give you a list of orthogroups that did not have duplications and in which all species are present.

```diff
- 3.1.- How many orthologous groups has orthoFinder been able to build?

- 3.2.- How many genes have been left outside the orthologous groups?
```

4.- Now go to the folder *Comparative_Genomics_Statistics*. In this folder you will find files summarizing the statistics of orthogroups for the different species.

```diff
- 4.1.- How many genes have been placed within the orthogroups? And which percentage do they represent from the total?

- 4.2.- Which is the species with the highest percentage of genes that are species specific? Do you think this is due to a biological reason?
```

## Exercise 2

1.- `OrthoFinder` has few parameters, but the most important one of them is the inflation parameter. This parameter indicates whether the orthogroups are going to be smaller or bigger. By default it is set to 1.5. We are now going to run `OrthoFinder` with a bigger inflation parameter. 

Replace the running line in the script *orthofinder.run* with (**remember to also change job name!**):

`orthofinder -b proteomes/OrthoFinder/Results_Jul06/ -I 3.0 -og`

> Note that we are using `-b` instead of `-f` and we are providing previously calculated results, this will avoid having to re-calculate the all-vs-all comparison. Also we are changing the inflation parameter using `-I` and setting it to 3.0. At this point we are only interested in comparing the orthogroups, the `-og` parameter will stop the run of orthoFinder after it calculates orthogroups. This is a time-saving trick if you want to assess different inflation parameters and how they affect your orthogroups.

This will generate a second folder which will be called *Results_Jul06_1* where the new results of `OrthoFinder` can be found.

2.- Repeat this step with an inflation parameter of 5.0

3.- Now compare the results obtained by the three runs of `OrthoFinder` and try and answer the following questions:

```diff
- 3.1.- Are the number of orthogroups the same in all three runs? How many orthogroups are in each of them?

- 3.2.- Did the change in inflation parameter affect the detection of single copy genes?

- 3.3.- Note that in order to do a good phylogenomics analysis you need to find groups of orthologous genes that are present in all species without duplications. Which file would give you this information and how many orthogroups can you use in each run?

- 3.4.- Orthogroups have a direct impact on orthology prediction. Depending on the analysis you are running, a reviewer may ask how the inflation parameter can affect your results. Can you think on a way to show which is the correct inflation parameter?
```

## Exercise 3

Orthogroups can contain duplications which means we can have a mix of orthologs and paralogs. `OrthoFinder` implements a method to distinguish between them. This time we will ask `OrthoFinder` to start from the orthogroups and to finish the process till the end.

Replace the running line in the script *orthofinder.run* with (**remember to also change job name!**):

`orthofinder -fg proteomes/OrthoFinder/Results_Jul06/`

> In this case the -fg option is telling `OrthoFinder` to use the orthogroups previously calculated to run the analaysis for orthology

Now, to correctly predict orthologous relationships a species tree is needed as reference. `OrthoFinder` tries to calculate the species tree on its own but you should always make sure the tree it has inferred is correct and rooted properly. Go to the results folder (proteomes/OrthoFinder/Results_Jul06_3) where you should see a folder called *Species_tree* and in there is a file called *SpeciesTree_rooted.txt*

1.- Check the species tree that has been automatically build by `OrthoFinder`. You can download the file to your personal computer and use `phylo.io` (
http://phylo.io/) to visualize the tree:

```diff
- Based on your knowledge on how bears have evolved, is the tree correct?
```

2.- If the tree is not rooted correctly, re-root the tree by clicking on the branch and pressing on re-root (make sure you press on the branch and not on the species name!). Now export the newick (find the button on the upper right part of the image, press on export nwk). Save the file into your folder and run orthofinder again with the option -s speciesTree_file.

Replace the running line in the script *orthofinder.run* with (**remember to also change job name!**):

`orthofinder -fg proteomes/OrthoFinder/Results_Jul06/ -s speciesTree_file`

> Note that if you have a species tree before running `OrthoFinder` it is more convenient to just provide it from the start using the `-s` option. Just make sure that the name of the proteome files are the same as the ones found in the species tree. 

> Other things that will affect the prediction of orthologs is how the gene trees are build. By default orthoFinder uses distance matrices and fastME to build the gene trees. This, while fast, can give faulty gene trees when dealing with more complex datasets. Once you are sure that your species tree is correct and that you are satisfied with the orthogroups, it is recommended that you do the orthology prediction using multiple sequence alignments (`mafft` is implemented in `OrthoFinder`) and make the gene trees using `IQ-TREE` (also implemented in `OrthoFinder`). Due to time constraints we are not going to run these, but you can find the pre-calculated results in the results folder.

3.- Go to the *Comparative_Genomics_Statistics*, here you can find main statistics for the analysis you have run. Search for the following information:

```diff
- 3.1.- Which two species have the highest number of orthologs?

- 3.2.- Which kind of orthologous relationships are most common between bears (one-to-one? many-to-one?). Can you think of a scenario in which this could be different?
```

4. We will analyse the orthogroup OG0000001 which should contain protein A0A7N5K5T5_AILME. If for some reason it does not search for the orthogroup that contains this protein and analyse that one.

```diff
- 4.1.- Go to the first set of results and check how many members this family has. Can you tell only from the orthogroup information which proteins are orthologs and which ones are paralogs?

- 4.2.- Now go to the last results. Go to the Orthologues folder and search for the orthologs to A0A7N5K5T5_AILME. How many do you find?

- 4.3.- Are they all one-to-one orthologs? ```

- 4.4.- Search for the gene tree of this family, copy the newick and visualize it in phylo.io. With the tree next to you, search your results for information on duplication events. How many duplication events can you find? Are any of them specific for a single species.

- 4.5.- Most of the duplications observed in the previous exercise were ancient, why do you think orthoFinder did not separate them? Were they separated in the analysis run with -I 3.0? If they were, and looking again to the gene tree, did the split make sense?
```


