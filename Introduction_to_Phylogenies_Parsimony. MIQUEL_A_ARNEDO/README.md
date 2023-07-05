# PARSIMONY: exercises for TNT

**Instructor:** Miquel A. Arnedo

## Software

<ins>**TNT**</ins> is available at http://www.lillo.org.ar/phylogeny/tnt/

**OTHER PROGRAMS:** 
+ <ins>**Mesquite:</ins>** Matrix and tree manipulation - http://www.mesquiteproject.org/Installation.html
+ <ins>**FigTree:</ins>** Tree visualization and manipulation - https://github.com/rambaut/figtree/releases
+ <ins>**Sequence Matrix:</ins>** Assemble multigene/multitaxon phylogenetic datasets easily - http://www.ggvaidya.com/taxondna/

**CODING GAPS AS ABSENCE/PRESENCE DATA:**
+ <ins>**SeqState:</ins>** http://bioinfweb.info/Software/SeqState

> **NOTICE that TNT comes in two flavours, command line with executables for Widows, Linux and Mac OSX and with a Shell with windows only for Windows. If you want to run the version with windows on Mac OSX or Linux you will need a Windows emulator (Parallels, VirtualBOX) or WINE (highly recommended)**

For this exercise we will use the command line executables
Download TNT form the website link above.

These are some wikis and websites that can help you:
+ http://phylo.wikidot.com/tntwiki 
+ http://phylobotanist.blogspot.com/2015/03/parsimony-analysis-in-tnt-using-command.html

---

You can find additional information about TNT in its web page (see above). There is a general description for running the program in Windows in a powerpoint presentation (QuickTutorial). You can get online help from all instructions with the help command. Examples of scripts come with the download zip. 
Finally, there is a TNT wiki (http://tnt.insectmuseum.org) where you can find additional information and wiki tutorials and a command list at http://phylo.wikidot.com/tntwiki#TNT_scripts

## Data files

We will investigate the phylogenetic relationships of bears—Family [Ursidae](http://tolweb.org/Ursidae/16015) by using a concatenated matrix of 6 mtDNA genes, namely cox1, cytb, nad1, nad5, 16S and 12S and and 5 nuclear genes, MC$R, RNASE1, CHRNA1, FES, IRBP, obtained from Genbank.

TNT uses it own format. Supposedly, also accepts NEXUS format, although it is a little bit picky! (e.g. use the simple NEXUS format in Mesquite). Having the sequences in fasta format as well will help you to recode the gaps as absence/presence characters (with Seqedit).

You can download (to local or cluster) the dataset in different formats using the following instructions (alternatively download it from the Campus Virtual at (https://campusvirtual.ub.edu/mod/folder/view.php?id=4480052)
```
svn export https://github.com/ppgcourseUB/ppgcourse2023/trunk/Introduction_to_Phylogenies_Parsimony.%20MIQUEL_A_ARNEDO/bears_2023_data
```
Below you will find the instruction to conduct some basic analyses using parsimony as an optimality criterion for phylogenetic inference. We have included the instructions and commands for running analyses using the terminal (Linux, MAC and windows).

You can find complete information on the commands used by TNT at http://phylo.wikidot.com/tntwiki#TNT_scripts

**Some interesting commands:**

```
RDIR set path for "run" files
STATS.RUN: Script that calculates RI and CI. 
XWIPE: Remove data set from memory (allows changing settings) 
BEST: filter trees, discarding suboptimal 
N   keep trees up to N steps (weighted) worse than best 
-   invert the selection criterion 
[   discard trees not fulfilling constraints of monophyly 
]   discard trees fulfilling them 
*   collapse the trees before comparing 
KEEP N number of trees in memory
CDIR: changes local directory
TAXNAME +N: increases the length of the taxa labels to N characters
```


For the terminal commands, we used some of the explanations provided by [Alexander Schmidt-Lebuhn](https://www.anbg.gov.au/cpbr/staff/schmidt-lebuhn-alexander_staff.html) available at http://phylobotanist.blogspot.com/2015/03/parsimony-analysis-in-tnt-using-command.html
 
**Although you can run any analyses in TNT using sequential commands, the easiest and fastest way is to wrap up a series of commands in small scripts using a text editor**. Then, name them using the extension .run. You can run the scripts within TNT using ```run file_name.run```. 

As an example, TNT includes a “one-shot analysis” that can be run with the script "aquickie.run". It includes a search strategy and node support using resampling. 

To run it from the command line type

```
proc filename.tnt; 
taxname= aquickie;
```

## A. Search for the shortest tree with TNT

TNT includes three different search strategies, namely **implicit enumeration** (exhaustive search, it does guarantee the shortest tree, but there is a taxa limit, ~15), **heuristic** (does not guarantee to find the shortest tree), and **new search strategy**, optimized algorithm for heuristic searchers of large data matrices (e.g. >100 taxa).
	
1. A useful feature in TNT is the possibility of creating a log file where all commands, analyses, and results can be saved for future review. Start logging the analysis output into a text file with a log filename.

```log file_name_log.out;```, and ```log /;```to close it


2. TNT is set to ridiculously low memory usage and will thus generally throw an error if you attempt to import a large matrix

Enter mxram 200 to set memory usage to 200 MB or whatever is realistic and necessary ```mxram 2000;```


3. Read the data matrix: 

```procedure input_file_name.tnt; ``` you can abbreviate the commands to the maximum unambiguous letter, here  ```p input_file_name.tnt; ```


+ For <ins>**exhaustive searches</ins>**

```IENUM;```search algorithm, results guaranteed to be optimal but feasible only for small data sets (<15 terminals)


+ For <ins>**heuristic searches</ins>**

1. Set the maximum number of trees to be stored in RAM (at least, the product of the number of RAS by the number of trees kept by replication) ```hold 1000;```
2. Select the search strategy (ej: RAS+TBR) ```mult 1000 =tbr h 10;```

TNT reports the number of trees found and the steps (Best Score), and the number of replicates (repls) that found the shortest tree. If the repls that found the Best score tree are less than 20% of the total repls. increase the number of repls. (e.g. by an order of magnitude). If any of the rpls. has reached the maximum number of Tree to save per replication (overflow) then, conduct a new round of searching but selecting the trees in RAM

To visualize trees:```tplot;```
```tplot N;``` show tree(s) N
```tplot *N;``` show tree(s) N in parenthetical notation
```naked=;``` show tree diagrams without numbers (default)
```naked-;``` show tree diagrams with numbers 

Repeat the operation to return to the previous screen or press ESC

To select an outgroup (in case it is not the first taxon in the matrix):

```outgroup N;``` where N is the order of the taxon
```taxname=;```turns on the name labels, so you can use the terminal name instead
```outgroup Ailuropoda_melanoleuca;```
```reroot;``` enforces the new root

If you found more than one shortest tree, make the strict consensus

```nelsen*;``` Calculate strict consensus tree, * keep consensus as last tree in memory

To save the trees in a file in TNT format

```tsave file_name.tre;``` open tree files (with + at the end, append)
```save;```save trees to file (previously opened with "tsave"), 
```save N;``` save N tree(s)
```save /``` close tree file

Alternatively, export the tree in NEXUS format to open them with FIGTREE or any other tree visualizing tools. We will use the command line. Write the instruction  

```taxname =;``` To record the names of the taxa instead of numbers
```export - file_name.tre``` Export the tree in NEXUS format. Nnote the use of "-", otherwise it exports the data matrix

If you want to export the trees with branch lengths in nexus format 
```
ttags = ;
blength *;
ttags );
export> consensus_tree.tre ;
ttags -;
```
Note that the "ttag" command controls the information that is added as labels to the tree nodes/branches

If you have generated a consensus tree and you only one to export the consensus tree, use the following

```tchoose /;``` select the last memory tree, discard the rest, and then proceed as above



+ For <ins>**New technology Search</ins>**

The heuristic strategy is highly inefficient and only justified when you have a small number of terminals (>100). Alternatively, the best strategy is to combine a selection of tree search shortcuts collectively known as New Technology Searches. Below you can find an example of a combination of different search strategies optimized for large datasets (>100 taxa). It conducts a new technology tree search combining ratchet, fuse, drift, and automatic stopping criteria, runs an additional round of branch swapping on shortest trees, makes the consensus, and saves it in NEXUS with branch lengths

```
xmult=hits 10 noupdate nocss replic 10 ratchet 10 fuse 1 drift 5 hold 100 noautoconst keepall;
bbreak = tbr ;
nelsen * ;
tchoose / ;
ttags = ;
blength *;
ttags );
ttags ;
export> consensus_tree.tre ;
ttags -;
```
 
## B. Implementing alternative gap treatments

TNT considers the **GAPS as 5th state by default**. If you want to treat them as missing data (for example if you have recorded them as absence/presence data), **before** opening the matrix select:

```
nstates NOGAPS;
```

• To get back to gaps as 5th state

```
nstates GAPS;
```

• To tret gaps as **absence/presence characters**, you have to recode them first. We can use the program **SeqState** (graphic interface in Java, mutliplatform)

<ins>**SeqState</ins>**

1. Double-click on the SeqState Java executable
2. Go to the top menu

*File> Load FASTA File*

3. Navigate to the file that you want to open (NOTICE: SeqState opens either NEXUS or FASTA format, however out of experience FASTA works better)
4. On the top menu, go to

*IndelCoder> Simmons and Ochoterena (2000) – Simple coding*

It will automatically generate an output file in NEXUS format with the name “input_file_name**_sic.nex**”

> **ATTENTION, this format is not directly readable by TNT, we must edit it using a standard text editor, to convert it into TNT format. Please, not that note DNA, protein, and numeric coding can coexist in the same matrix, but they have to be entered sequentially and the format specified with the following notation: &[dna], &[num], &[prot] for DNA sequence data, morphology/binary and aminoacid, respectively.**

**In our example we will make the following changes:**

```xread
number_of_chars number_of_taxa
&[dna]
‘DNA matrices’
&[num]
‘gaps scores as absence/presence’
;
cc-.;
proc/;
```

> **ATTENTION, you must choose the option Read gaps as missing if you will not be counting the gaps twice, as a 5th state and as an absence/presence!**

Proceed as explained for tree searching to find the shortest tree 

## C. Implied weights (Goloboff 1993)

If you want to change the value of K, first select the value (weighting function, any value from 1 to 100, the higher the lesser homoplasy is weighted against. Then check using the implied weights box. 

```
piwe=50;
```
Search for the fittest tree using the abovementioned search strategies for parsimony
```
fit;
```
reports the fit of the tree(s)
If you want to know what the actual length of the fittest tree is found use
```
length;
```
To deactivate implied weights
```
piwe-; 
```
 
## D. Node support: resampling

You can use different resampling methods: boot=bootstrap, jak=jackknife, sym=symetrical resampling (recommended when using differential weighting)

An example of assessing node support using Jackknife
```
ttags -;
ttags =;
resample jak replications 1000 [ mult 15 =tbr hold 20 ];
keep 1 ;
ttags );
ttags ;
export - boots.tre; proc/;
ttags -;
```

Note: by default, TNT provides two support indexes, the absolute frequencies and the ‘Group present/Contradicted’, which measures the support as the difference in frequency between the group
and its most frequent contradictory group. GC values of -1, 0, and 1 indicate (respectively) maximum contradiction, indifference, and maximum support.

## E.  Node support: Bremer support

1. First, we will search for suboptimal trees, that is, with a number of steps above the minimum, to define the upper limit where the clade is not found anymore
2. Define the maximum number of suboptimal trees to keep (eg = 10,000)
3. Define the maximum extra steps to consider (e.g. = 50)
4. Start searching for suboptimal trees
5. Calculate Bremer supports

An example for assessing Bremer support to  tree nodes. Perform branch-swapping, using pre-existing trees as starting point. Alternatively, start by running a new, simple run, eg: mult 1000 =tbr h 10;
```
ttags-;
ttags=;
hold 10000 sub 50 bbreak = tbr ;
collapse 3;
bsupport;
keep 1;
ttags );
ttags ;
export – bremer.tre;
ttags-;
```
Alternatively, to avoid collapsing the RAM with very suboptimal trees, conduct the suboptimal search step by step, as follows:
```
hold 1000 ; sub 1 ; find * ; <enter>*</br>
hold 2000 ; sub 3 ; find * ; <enter>*</br>
hold 4000 ; sub 5 ; find * ; <enter>*</br>
collapse 3;
bsupport;
```

To conduct Partition Bremer Support you can use a TNT script (pbsup.run) written by Carlos Peña and available at https://github.com/carlosp420/pbsup.run. Please refer to that GitHub site to conduct the analysis. You have all the files required to run the script available in the bears_2023_data/PBS folder  

 
## F. Partial analyses

To analyse only part of the characters use the command ```blocks```;

To display the defined partitions
```blocks;``` 
To define blocks to start at character(s) J K L
```blocks J K L;```
To save partitions
```blocks*; ```
To disable all characters that are not listed in Nn (if the list is headed by the "&", only the shared partitions are kept active)
```blocks= Nn;``` 

In our bear example, partitions correspond to the following genes and positions:
nuclear genes<enter>
 	CHRNA1 = 1-362
	FES = 363-831
 	IRBP = 832-2111
  	MC4R = 2112-3110
   	RNASE1 = 3111-3530
mt genes
	mt12S = 3531-4517
 	mt16S = 4518-6107
  	mtCOX1 = 6108-7652
   	mtCYTB = 7653-8787
    	mtNAD1 = 8788-9743
     	mtNAD5 = 9744-11573

If I wish to analyze the nuclear data only:
```
blocks 3531;
blocks*;
blocks= 0;
``` 

If, instead, I just one to analyze the mitochondrial genes, simply:
```blocks= 1;```

## G. ILD incongruence test

The Incongruence Length Difference (ILD) test can detect cases of inconsistency between data partitions and can be implemented in TNT through the ILDtntk.run script (http://phylo.wikidot.com/tntwiki#TNT_scripts). Please, note that the script can only compare two partitions simultaneously.

```
run ildtntk.run XX YY ZZ;
```
runs the script, if the script is not in the same folder as the TNT executable, then write the full path
XX: defines the number of characters of the first partition. 
YY: the number of replications to do, the more the better, minimum 1000, but you can use 99 to speed up the analysis
ZZ the search algorithms to use when calculating length for each partition in each replication (default: 1 random addition sequence + SPR), e.g:

```mult 15 =tbr hold 20```

For instance, to test possible incongruence between the nuclear and the mitochondrial genes use the following command line (1000 permutations)

```run ildtntk.run 3530 8043 1000 mult 15 =tbr hold 20;```
