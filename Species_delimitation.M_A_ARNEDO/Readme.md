# Species delimitation

**Instructor:** Miquel A. Arnedo \& Vanina Tonzo

## Programs

+ ABGD, available at http://wwwabi.snv.jussieu.fr/public/abgd/
Automatic identification of the barcoding gap (Puillandre et al., 2012). ABGD requires users to provide several parameters: the genetic distance (for example, Jukes-Cantor or p-distances), a prior limit on intraspecific diversity (P), and a proxy for minimum width of barcoding gaps (X). For each set of user parameters, it returns two constraints: a "primary partition" and a "recursive" one obtained after applying its algorithm recursively.
ATTENTION: This approach has recently been superseded by ASAP (Assemble species by automatic partitioning; Puillandre et al., 2021), available at https://bioinfo.mnhn.fr/abi/public/asap/.
+ GMYC can be run in R environment with SPLITS package package v 1.0-19 avialbale at https://rdrr.io/rforge/splits/. 
It can also be run on-line at: https://species.h-its.org/gmyc/
There is a Bayesian version bGMYC that allows uncertainty measures to be assigned to the constraint: R package available at: https://nreid.github.io/software/
+ (m)PTP , available at https://github.com/Pas-Kapli/mptp
It uses a dynamic programming implementation that estimates the ML delimitation faster and more accurately than the original PTP. It assumes a different exponential distribution for branching events for each of the bounded species, allowing it to fit a wider range of empirical data sets. It also incorporates the option of determining minimum branch length for your own data sets.
The PTP-like methods family can be run online at:
PTP/mPTP: https://mptp.h-its.org/#/tree
bPTP: https://cme.h-its.org/exelixis/web/software/PTP/index.html
+ BPP v. 4, available at http://abacus.gene.ucl.ac.uk/software.html
It incorporates the multiple species coalescent model (MSC) and allows (1) the estimation of the parameters of divergence times of species and population sizes when the phylogeny of the species is given (Rannala and Yang, 2003), (2) the Species tree inference using user-supplied individual-to-species assignments (Rannala and Yang, 2017), (3) species delimitation using a user-specified guide tree (Yang and Rannala, 2010; Rannala and Yang, 2013 ), (4) the joint delimitation of species and estimation of species trees (Yang and Rannala 2014). 
The last version of BPP now also implements the multispecies-coalescing-with-introgression (MSci) model (see Flouri et al, 2020), an extension of the multispecies coalescing model to incorporate introgression/hybridization.

## Data

+ Single marker discovery
    1. bears_c1.fasta
mtDNA cox1 gene from 90 taxa (1,545 bp) in FASTA format of specimens belonging to the Ursus genus, using the Sloth bear as outgroup (NC009970_Melursus_ursinus)
    2. bears_c1_root.treefile
Preferred ML tree inferred with `IQ-TREE` (full codon partition)

+ Multi-locus validation
    1. bears_bpp.txt
12 protein-coding genes from 89 taxa
    2. bears_map.txt
An individual-to-population map file (Imap file) that assigns every individual to a population
    3. bears_A10.bpp.ctl 
Control file for BPP 
    4. bears_out_2.txt, bears_mcmc_2.txt 
BPP output files (outfile, mcmcfile)

## Aim

Become familiar with the use of different types of programs for species delimitation. Perform different types of quantitative species delimitation analyses in bears (_Ursus_) using single and multilocus approaches from unique markers. 

## Procedure

## 1. Single marker discovery

**1.1. Automated Barcoding Gap Discovery**

We will use the online tool ABDG, available at http://wwwabi.snv.jussieu.fr/public/abgd/abgdweb.html

Load the file (bears_c1.fasta) by clicking on the Choose file option. The program requires the matrix in fasta format.

You can either select an alignment in FASTA format or a distance matrix in MEGA format. ABDG only includes the JC69 and K2P distance models if you want to use another model (eg. p-distance) then use MEGA to calculate the distance matrix and load it in ABDG.

Use the parameter values by default and click on the “GO” button. You will be referred to a new page that shows the calculations. Click on here to see the results, which includes the Histogram of distances, the Ranked distances and the number of groups depending on the prior intraspecific distances, both primary and recursive. 

If you click on the dots of the last graphic, it will show you the actual groupings. You can download and visualize the inferred groups in a tree like format. Copy&paste the resulting file into FigTree to visualize the groups.

How many clusters have you found between 1 to 3% of prior interspecific pairwise divergence?


**1.2. Evolutionary criteria: the GMYC model**

To use this approach, we first need to infer a phylogenetic tree form our data, then transform it to ultrametric and finally apply the GMYC model. There are several ways to obtain an ultrametric tree:

1. Inferring first a tree with branches proportional to substitutions using either max. likelihood (eg. `iqtree`) and then transforming the tree to ultrametric using programs such as `PATHd8` or `R8S` (penalized likelihood).

2. Using `BEAST`, which incorporates the time estimations as a parameter of the analysis. If you choose this option, it is better to use a coalescent prior, to minimize type 1 error.

We will use the first option and will run all the analyses using a **R script (*bears_gmyc.R*)**. This script has the following dependencies: R packages `splits`, `paran`, `mass` and `rncl`. You can either run this analysis in your personal computers, and in that case you must install the packages in your laptops before running the script, or, if the analysis is not computationally demanding, you can run it exceptionally in the control node using these commands:

```bash
module load r-ape r-mass r-rncl
```

Then, in the R console, type:

```R
install.packages("paran", dependencies=TRUE)
install.packages("splits", repos="http://R-Forge.R-project.org", dependencies=TRUE)
quit()
```

After doing that, you can now run the script bears_gmyc.R:

```bash
Rscript ./sp_data/bears_gmyc.R
```

`splits` incorporates support values from the GMYC clusters, based on the AIC. The function `gmyc.support` calculates support values of the GMYC-delimited species by using the multimodel comparison approach described by Burnham & Anderson (2002). The support value of a node is defined as the sum of Akaike weights of candidate delimitation models in which the node is included. Only models included in the p% confidence set obtained by `confset.gmyc` are used for calculation. For further details, see *Fujisawa, T., & Barraclough, T. G. (2013). Delimiting species using single-locus data and the Generalized Mixed Yule Coalescent approach: a revised method and evaluation on simulated data sets. Systematic Biology, 62(5), 707–724. http://doi.org/10.1093/sysbio/syt033.*


**1.3. Evolutionary criteria: the mPTP model**

Unlike GMYC, `(m)PTP` methods do not require an ultrametric tree. Therefore, you can proceed with your ML or BI preferred tree(s). `mPTP` implements two flavours of the point-estimate solution. First, it implements the original method from Zhang et al. (2013) where all within-species processes are modeled with a single exponential distribution. `mPTP` uses a dynamic programming implementation which estimates the ML delimitation faster and more accurately than the original `PTP`. The dynamic programming implementation has similar properties as (Gulek et al. 2010). See the github of the program at https://github.com/Pas-Kapli/mptp/wiki for more information. The second method assumes a distinct exponential distribution for the branching events of each of the delimited species allowing it to fit to a wider range of empirical datasets.

MCMC method `mPTP` generates support values for each clade. They represent the ratio of the number of samples for which a particular node was in the between-species process, to the total number of samples.

You have different alternatives to run `(m)PTP`. There is an online service available at. It is friendly but does not include all the options. We will use a command line version instead.

+ Here you have the commands to implement mPTP with MCMC support in a batch file to submit the jobs in the cluster (**Please do not launch the `mptp` command from the control node**):

First run to get the best_minbr value:

```bash
#!/bin/bash

##Script to submit bpp jobs

#SBATCH -p normal
#SBATCH -c 8
#SBATCH --mem=6GB
#SBATCH --job-name mptp-1
#SBATCH -o %j.out
#SBATCH -e %j.err

# modules
module load mptp

# running the program
mptp --tree_file ./sp_data/bears_c1_root.treefile --minbr_auto ./sp_data/bears_c1.fasta --output_file bears_c1_best_minbr.out
```

> --mcmc Support values for the delimitation    
> --mcmc_runs Perform multiple MCMC runs   
> --mcmc_log Log samples and create SVG plot of log-likelihoods      
> --mcmc_burnin Ignore all MCMC steps below threshold    

Additional options
--mcmc_startnull	Start each run with the null model (one single species)   
--mcmc_startrandom	Start each run with a random delimitation (default)  
--mcmc_startml	Start each run with the delimitation obtained by the Maximum-likelihood heuristic

Main run:

```bash
#!/bin/bash

##Script to submit bpp jobs

#SBATCH -p normal
#SBATCH -c 8
#SBATCH --mem=6GB
#SBATCH --job-name mptp-1
#SBATCH -o %j.out
#SBATCH -e %j.err

# modules
module load mptp

# running the program
mptp --seed 767 --multi --tree_file ./sp_data/bears_c1_root.treefile --outgroup NC009970_Melursus_ursinus --outgroup_crop --minbr 0.0006459066 --mcmc 50000000 --mcmc_startnull --mcmc_runs 3 --mcmc_log 1000000 --mcmc_burnin 2000000 --output_file bears_c1_nout_mptpt.out
```

> mptp: runs the executable  
> --seeds: random number generator  
> --tree_file: define tree file   
> --multi: use one lambda per coalescent (mPTP, this is default. Alternatively, one lambda for all coalescents –single = PTP)    
> --outgroup: define out group     
> --outgroup_crop: and remove it    
> --minbr: Set minimum branch length, you can select the value by running a preliminary analysis using the following command line:    
> --output_file: Name of the output file    

+ **\#Output Files For maximum likelihood delimitation:**   
  1. **output_filename.txt**: contains information about the run:     
  2. **output_filename.svg**: is a vector graphic of the phylogenetic. The ML delimitation scheme is illustrated in the graphic; with black are illustrated the branches of the speciation process and with red the branches of the coalescent process.    

+ **\#Output Files For a single mcmc run, four files will be created:**   
  1. **filename.run_seed.stats**: This file reports the frequency of all possible number of species for the input phylogeny (ie. `1 -n`, where n is the number of tips in the phylogeny).   
  2. **filename.run_seed.svg**: This file corresponds to a graphical representation of the phylogenetic input tree. The support values for each node to be part of the speciation process is also provided in this file. The branches of the tree are colored with a gradient from black to red, scaled by the corresponding support values. If a branch is certainly part of the speciation process (i.e., the support value of its ascending node is 100%) it will be colored black and if branch is certainly part of the coalescence (i.e., the support value of its ascending node is 0) it is colored red.   
  3. **filename.run_seed.tree**. It contains the tree in format and the support values for the each node being part of the speciation process.   
  4. **filename.run_seed.txt**. It contains information about the run and the ML delimitation similar to the ML output.     


+ If the `--mcmc_log`option is activated, two additional files will be generated: 

  1. **filename.run_seed.log**, contains the likelihood values in a text format. 
  2. **filename.run_seed.logl.svg**, plots the likelihood for every MCMC sample stored in memory (set by the `--mcmc_log`option). 
    
+ Finally, if multiple mcmc runs are executed (i.e.,`--mcmc_runs`> 1), then all of the output files described above will be created for each of the independent runs. Two additional files will be created:

  1. **filename.run_seed.combined.svg**. A graphical representation of the phylogenetic input tree with the support values derived from all independent MCMC runs.
  2. **filename.run_seed.combined.tree**. The phylogenetic input tree in *newick* format with the support values derived from all independent MCMC runs.   

### Monitoring for chain convergence:

The ASDDSV is inspired by the standard deviation of split frequencies (Ronquist et al., 2012) and it is used for quantifying the similarity among independent MCMC runs. To calculate it we average the standard deviation of per-node delimitation support values across the independent runs. ASDDSV approaches zero as runs converge to the same distribution of delimitation.   
     
In the end of the MCMC runs the ASDDSV will be printed at the end of the screen as follows:

Average standard deviation of support values among runs: 0.001385

## 2. Multi-locus species delimitation: `BPP` program

`BPP` requires 3 files to run:

* **The sequence file**:    
It uses a sequential phylip format with one gene (locus) after the other. Locus may differ in the number of sequences. Labelling is important. Each sequence name must be followed by ^ and an ID corresponding to the individual the sequence comes from.
* **The IMAP file**:   
An individual-to-population map file (Imap file) that assigns each individual to a population.
* **The control file**:    
Specifies the model and the parameters and priors used, and effectively “drives” the analysis. Lines beginning with an asterisk are comments. BPP implements several types of analyses other than the species delimitation.

   + A00 (speciesdelimitation = 0, speciestree = 0): Estimation of the parameters of species divergence times and population sizes under the MSC model when the species phylogeny is given (Rannala and Yang, 2003).
   + A01 (speciesdelimitation = 0, speciestree = 1): Inference of the species tree when the assignments are given by the user (Rannala and Yang, 2017).
   + A10 (speciesdelimitation = 1, speciestree = 0): Species delimitation using a user-specified guide tree (Yang and Rannala, 2010; Rannala and Yang, 2013);
   + A11 (speciesdelimitation = 1, speciestree = 1): Joint species delimitation and species tree inference, or unguided species delimitation (Yang and Rannala, 2014).

Additionally, a file specifying the heredity scalar (eg. 0.25 for mtDNA or 1 for diploid nucDNA) may be included.

To run the program in the cluster, submit this batch file using `sbatch` (**Please do not launch the `mptp` command from the control node**):

```bash
#!/bin/bash

##Script to submit bpp jobs

#SBATCH -p normal
#SBATCH -c 8
#SBATCH --mem=6GB
#SBATCH --job-name bpp-1
#SBATCH -o %j.out
#SBATCH -e %j.err

# modules
module load bpp

# running the program
bpp --cfile sp_data/BPP/inputs/bears_A10.bpp.ctl 
```
See BPP_tutorial.pdf for further details.
