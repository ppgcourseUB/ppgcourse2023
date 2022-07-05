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

TNTA uses it own format. Supposedly, also accepts NEXUS format, although it is a little bit picky! (e.g. use the simple NEXUS format in Mesquite). Having the sequences in fasta format as well will help to recode the gaps as absence/presence characters.

Below you will find the instruction to conduct some basic analyses using parsimony as an optimality criterion for phylogenetic inference. We have included the intstructions and commands for running analyses using the graphic interface (only available for Windows) and the terminal (Linux, MAC and windows).

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
KEEP N number of trees in memeory
CDIR: changes local directory
TAXNAME +N: increases the length of the taxa lables to N characters
```


For the terminal commands we used some of the explanations provided by by [Alexander Schmidt-Lebuhn](https://www.anbg.gov.au/cpbr/staff/schmidt-lebuhn-alexander_staff.html) available at http://phylobotanist.blogspot.com/2015/03/parsimony-analysis-in-tnt-using-command.html
 

## A. Search for the shortest tree with TNT

TNT includes three different search strategies, namely **implicit enumeration** (exhaustive search, it does guarantee the shortest tree, but there is a taxa limit, ~15), **heuristic** (do not guarantee finding the shortest tree) and **new search strategy**, optimised algorithm for heuristic searchers of large data matrices (e.g. >100 taxa).
	
1. A useful feature in TNT is the possibility of create a log file where all commands, analyses and results can be saved for future review. Start logging the analysis output into a text file with log filename.

```
> log file_name_log.out;
> log /;
#to close it
```

2. TNT is set to a ridiculously low memory usage and will thus generally throw an error if you attempt to import a large matrix

```
>mxram 2000;
##Enter mxram 200 to set memory usage to 200 MB or whatever is realistic and necessary
```

3. Read data matrix: 

```
>procedure input_file_name.tnt; 
##you can abbreviate the commands to the mimum unambiguos letter , here p input_file_name.tnt
```

+ For <ins>**exhaustive searches</ins>**

```
>IENUM;
##search algorithm, results guaranteed to be optimal but feasible only for small data sets (<15 terminals)
```

+ For <ins>**heuristic searches</ins>**

4. Set the maximum number of trees to be stored in RAM (at least, the product of the number of RAS by the number of trees kept by replication)
5. Select the search strategy (ej: RAS+TBR)
6. Start searching

```
>hold 1000;
>mult 1000 =tbr h 10;
```

7. TNT reports the number of trees found and the steps (Best Score), and the number pf replicates (rpls) that found the shortest tree. 
8. If the repls that found the Best score tree are less than 20% of the total repls. increase the number of repls. (e.g. by an order of magnitude) 
9. If any of the rpls. has reach the maximum number of Tree to save per replication (overflow) then, conduct a new round of searching but selecting the trees in RAM
10. To visualize trees:

```
>tplot;
##N    show tree(s) N, *N   show tree(s) N in parenthetical notation. 
>naked=;
## =   show tree diagrams without numbers (default),   -   with
```

11. Repeat the operation to return to the previous screen or press esc

12. To select outgroup (in case it is not the first taxon in the matrix):

```
>outgroup N;
##where N is the order of the taxon
>taxname=;
##turns on the name labes, so you can use the terminal name instead
>outgroup NC_009492_Ailuropoda_melanoleuc;
>reroot;
## enforces the new root
```

13. To save the trees in a file
14. The file remains open and will store any other tree obtained.
15. If you want to stop saving trees, close the file

```
>tsave file_name.tre;
##open tree files (with + at the end, append)
>save;
##save trees to file (previously opened with "tsave"), 
>save N;
##save N tree(s) 
>save /
##close tree file
```

16. If you found more than one shortest tree, make the consensus

```
>nelsen*;
##Calculate strict consensus tree, * keep consensus as last tree in memory
```

17. Alternatively, export the tree in NEXUS format to open them with FIGTREE. We will use the command line. Write the instruction taxname =; To record the names of the taxa instead of numbers, and then export file_name;

```
> tchoose /;
##select last memory tree, discard the rest 
ttags = ;
blength *;
ttags );
export> consensus_tree.tre ;
ttags -;
##save the consensus tree with branch lengs in nexus format 
```
+ For <ins>**New technology Search</ins>**

This is a collection of strategies that can be combined in different ways. A possible example:
• Check driven search
• Combine options (check box): sectorial search, tree fusing and tree drifting
• Initial addseqs= 50; 
• Initial level=15; 
• Stabilize consensus 5 times (with default factor of 75)

TNT includes a “one-shot analysis” that can be run with the script "aquickie.run". It includes a serach strategy and node support using resampling. 

To run it from the command line type

*>proc filename.tnt ; aquickie ; [enter].*

```
> xmult=hits 10 noupdate nocss replic 10 ratchet 10 fuse 1 drift 5 hold 100 noautoconst keepall;
>bbreak = tbr ;
>nelsen * ;
>tchoose / ;
ttags = ;
blength *;
ttags );
ttags ;
export> consensus_tree.tre ;
ttags -;
##conducts a new technology tree search combining ratchet, fues, drift and automatic stopping criteria, runs an additonal round of branch swapping on shortest trees, makes the consensus and save it in NEXUS with branch lenthgs 
```
 
## B. Implementing alternative gap treatments

TNT considers the **GAPS as 5th state by default**. If you want to treat them as missing data (for example if you have recoded them as absence/presence data), **before** opening the matrix select:

```
>nstates NOGAPS;
```

• To get back to gaps as 5th state

```
> nstates GAPS;
```

• To tret gaps as **absence/presence characters**, you have to recode them first. We can use the program **SeqState** (graphic interface in Java, mutliplatform)

<ins>**SeqState</ins>**

1. Double click on the SeqState java executable
2. Go to the top menu

*File> Load FASTA File*

3. Navigate to the file that you want to open (NOTICE: SeqState opens either NEXUS or FASTA format, however out of experience FASTA works better)
4. On the top menu, go to

*IndelCoder> Simmons and Ochoterena (2000) – Simple coding*

It will automatically generate an output file in NEXUS format with name “input_file_name**_sic.nex**”

> **ATTENTION, this format is not directly readable by TNT, we must edit it using a standard text editor, to convert it into TNT format. Please, not that note DNA, protein and numeric coding can coexist in the same matrix, but they have to be entered sequentially and the format specify with the following notation: &[dna], &[num], &[prot] for DNA sequence data, morphology/binary and aminoacid, respectively.**

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

> **ATTENTION, you must choose the option Read gaps as missing if you will not be counting the gaps twice, as a 5th state and as an absence / presence!**

Proceed as explained for tree searching to find the shortest tree 

## C. Implied weights (Goloboff 1993)

1. If you want to change the value of K, first select the value (weighting function, any value from 1 to 100, the higher the lesser weighted against is homoplasy. Then check using implied weights box. 
2. Search for fittest tree using the abovementioned searh strategies for parsimony
3. If you want to know what the actual length of the fittest tree is found use

```
> piwe=10;
## implements maximum implied weight as optimality criterion, with K=10;
> fit;
##reports the fit of the tree(s)
>length;
##reports the actual number of steps of the tree(s)
piwe-; 
##desactivates
```
 

## D. Node support: Bootstrap

```
>ttags -;
>ttags =;
>resample boot replications 1000 [ mult 15 =tbr hold 20 ];
##you can use different resampling methods: boot=bootstrap, jak=jackknife, sym=symetrical resampling
keep 1 ;
ttags );
ttags ;
export - boots.tre; proc/;
>ttags -;
```

## E.  Node support: Bremer support

1. First, we will search for suboptimal trees, that is, with a number of steps above the minimum, to define the upper limit where the clade is not found any more
2. Define the maximum number of suboptimal trees to keep (eg = 10,000)
3. Define the maximum extra steps to consider (e.g. = 50)
4. Start searching for suboptimal trees
5. Calculate Bremer supports

Alternatively, to avoid collapsing the RAM with very suboptimal trres, conduct the suboptimal search step by step, as follows:

*hold 1000 ; sub 1 ; find * ; <enter>*</br>
*hold 2000 ; sub 3 ; find * ; <enter>*</br>
*hold 4000 ; sub 5 ; find * ; <enter>*</br>
*bsupport;*

```
>ttags-;
>ttags=;
>hold 10000;
> sub 50;
> bbreak = tbr ;
##Perform branch-swapping, using pre-existing trees as starting point. Alternatively make a new, simple run, eg: mult 1000 =tbr h 10;
> bsupport;
>keep 1;
>ttags );
>ttags ;
>export – bremer.tre;
>ttags-;
```

> TIP: TNT allows you to write small scripts to automate certain tasks. As an example, there is a script that calculates Partition Bremer Support  (pbsup.run) written by Carlos Peña and available at https://github.com/carp420/pbsup.run.

Please refer to https://github.com/carlosp420/pbsup.run for any additional information about this routine.
 
## F. Partial analyses

1. To analyse only part of the characters

```
>blocks;
##; Displays the defined partitions
##*; Saves partitions
##= Nn Disables all characters that are not listed in Nn (if the list is headed by the "&", only the shared partitions are kept active)
```

## G. ILD incongruence test

1. The Incongruence Length Difference (ILD) test can detect cases of inconsistency between data partitions and can be implemented in TNT through the ILD.run script (http://phylo.wikidot.com/tntwiki#TNT_scripts). Please, note that the script can only compare two partitions simultaneously.

```
> run ildtntk.run XX YY ZZ;
##runs the script, if the script is not in the same folder as the TNT executable, then write the full path
##XX: defines the number of characters of the first partition. 
##YY: the number of replications to do, the more the better, minimum 1000, but you can use 99 to speed up the analysis
##ZZ the search algorithms to use when calculating length for each partition in each replication (default: 1 random addition sequence + SPR), e.g: mult 15 =tbr hold 20
```
