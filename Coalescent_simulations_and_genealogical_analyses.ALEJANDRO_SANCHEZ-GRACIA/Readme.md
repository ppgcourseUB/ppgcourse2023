# Coalescent simulations and genealogical analyses

**Instructor:** Alejandro Sánchez-Gracia

In this session you will use the program *fastsimcoal2* to:
* Simulate samples for **statistical model-based inference**
* Simulate genome-wide SNP data for **validating statistical methods**
* Simulate DNA sequences for **predictive uses**

You will use the program [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2) (fast sequential Markov coalescent simulation of genomic data under complex evolutionary models), developed by Excoffier, L. and Foll, M  (and subsequently expanded with other contributors) at the University of Berne. More information about the rationale of this software and the implemented methods can be found [here](https://www.ncbi.nlm.nih.gov/pubmed/21398675) and [here](https://www.ncbi.nlm.nih.gov/pubmed/24204310).

According to the website, *fastsimcoal2* is a simulation-based framework that "...*preserves all the simulation flexibility of simcoal2 but that is now implemented under a faster continuous-time sequential Markovian coalescent approximation, allowing it to efficiently generate genetic diversity for different types of markers along large genomic regions, for both present or ancient samples. It includes a parameter sampler allowing its integration into Bayesian or likelihood parameter estimation procedure.*

---------------------------------------------------

## *fastsimcoal2* general options

To see a full list of the options availalbe in `fastsimcoal2` type:

```bash
fsc2709 -h
```
> Remember that if you are in your server account, you should load the module first.

<details><summary>CLICK HERE to see the result of this command</summary>

<p>

```
fastSimcoal3 (ver 3.7.9 - 11.04.22)

Usage:

 -h  --help              : prints this help 
 -i  --ifile test.par    : name of parameter file
 -n  --numsims 1000      : number of simulations to perform
                           Also applies for parameter estimation
 -t  --tplfile test.tpl  : name of template parameter file (optional)
 -f  --dfile test.def    : name of parameter definition file (optional)
 -F  --dFile test.def    : same as -f, but only uses simple parameters defined
                           in template file. Complex params are recomputed
 -e  --estfile test.est  : parameter prior definition file (optional)
                           Parameters drawn from specified distributions are 
                           substituted into template file.
 -E  --numest 10         : number of draws from parameter priors (optional)
                           Listed parameter values are substituted in template file
 -g  --genotypic         : generates arlequin projects with genotypic data
 -p  --phased            : specifies that phase is known in arlequin output
                           default: phase is unknown
 -s  --dnatosnp 2000     : output DNA as SNP data, and specify maximum no. of SNPs
                           to output (use 0 to output all SNPs).
 -S  --allsites          : output the whole DNA sequence, incl. monomorphic sites 
 -I  --inf               : generates DNA mutations according to an
                           infinite site (IS) mutation model
 -d  --dsfs              : computes derived site frequency spectrum
                           (for SNP or DNA as SNP (-s) data only).
 -m  --msfs              : computes minor site frequency spectrum
                           (for SNP or DNA as SNP (-s) data only)
 -j  --jobs              : output one simulated or bootstrapped SFS per file
                           in a separate directory for easier analysis
                           (requires -d or -m and -s0 options)
 -b  --numboot 10        : number of bootstraps to perform on polymorphic sites to extract SFS
                           (should be used in addition to -s0 and -j options)
 -H  --header            : generates header in site frequency spectrum files
 -q  --quiet             : minimal message output to console
 -T  --tree              : outputs coalescent tree in nexus format
 -k  --keep  10000       : number of simulated polymorphic sites kept in memory
                           If the simulated no. is larger, then temporary files
                           are created. Default value is 10000
 -r  --seed              : seed for random number generator (positive integer <= 1E6)
 -x  --noarloutput       : does not generate Arlequin output
 -G  --indgenot          : generates an individual genotype table
 -M  --maxlhood          : perform parameter estimation by max lhood from SFS
                           values between iterations
 -L  --numloops 20       : number of loops (ECM cycles) to perform during
                           lhood maximization. Default is 20
 -l  --minnumloops 2     : number of  loops (ECM cycles) for which the lhood is 
                           computed on both monomorphic and polymorphic sites
                           if REFERENCE parameter is defined
 -C  --minSFSCount 1     : minimum observed SFS entry count taken into account in
                           likelihood computation (default = 1, but value can be < 1. e.g  0.5)
 -0  --removeZeroSFS     : do not take into account monomorphic sites for SFS
                           likelihood computation
 -a  --ascDeme 0         : This is the deme id where ascertainment is performed
                           when simulating SNPs. Default: no ascertainment.
 -A  --ascSize 2         : number of ascertained chromosomes used to define SNPs in
                           a given deme. Optional parameter. Default value is 2
 -u   --multiSFS         : generate or use multidimensional SFS
 -w   --brentol 0.01     : tolerance for Brent optimization
                           Default = 0.01. Smaller value imply more precise estimations
                           but require more computation time (min;max) = (1e-1;1e-5)
 -c   --cores 1          : number of openMP threads for parameter estimation
                           (default=1, max=numBatches, use 0 to let openMP choose optimal value)
 -B   --numBatches 12    : max. no. of batches for multi-threaded runs
                           (default=12)
 -P   --pooledsfs        : computes pooled SFS over all samples.
                           Assumes -d or -m, but not -u flag activated
      --recordMRCA       : records tMRCAs for each non recombining segment and outputs
                           results in <generic name>_mrca.txt. Beware: huge slow down of computing time
      --foldedSFS        : computes the 1D and 2D MAF SFS by simply folding the DAF SFS
      --logprecision 23  : precision for computation of logs of random numbers. Max value is 23
                           Default value is 23 (full precision). Recommended lower value is 18
      --initValues my.pv : specifies a file (*.pv) containing initial parameter values
                           for parameter optimization
      --nosingleton      : ignores singletons in likelihood computation
 -R   --repetition  1    : number of repetitions of current computations
                           default = 1
      --ASM              : estimates probability of Ancestral State Misidentification (ASM)
                           and corrects expected sfs accordingly
      --pureSMC          : use the pure SMC algorithm for recombinations
                           instead of the default SMC' algorithm
      --SMC2             : use a modified SMC' algorithm for recombinations
                           instead of the default SMC' algorithm
 -4   --F4GT             : computes distribution of failure of 4 gamete test
      --testValue  xxx   : some value that I can use for testing things
                            Alternative to Brent optimization (still under devt)
                           Lower values imply faster computations, but less precisions on logs
      --randParam 100    : computes SFS-based likelihoods for a given number of
                           random parameter combinations
 -Z                       : performs a single likelihood estimation with current tpl and est file
```
	

</p>

</details>

## Observed data

The folder "observed_data" contains the joint and multidimensional SFS calculated from 3000 independent SNPs in 32 samples of three Human and one Neanderthal populations, see [Figure1](figure1.pdf).

---------------------------------------------------

## Submitting SLURM job scripts

All job scripts available for this practical session MUST be submitted from the scripts folder.

## Example 1: Demographic inference based on SFS

One of the most powerful utilities of computer simulations based on the coalescent is **statistical inference of demographic histories and model parameters**. Here you will run a toy example in which you will compare two demographic histories (vaguely) inspired in Human evolution, and estimate the underlying parameters from the multidimensional frequency spectrum (multiSFS). You will use *fastsimcoal2* simulations to compute the expected SFS under each model and to estimate model parameters. To determine the best model (based on the composite likelihood estimated for each model), you can apply an [AIC](https://es.wikipedia.org/wiki/Criterio_de_informaci%C3%B3n_de_Akaike) (Akaike information criterium) approach.

In a real expermient, you sould create first the input files for *fastsimcoal2*, one per each competing model ([Figure1](figure1.pdf)). However, for purely practical reasons and due the limited duration of the hands-on sesion, we will provide us with the correct input files, ready to be used with the program. **Some very convenient homework to understand the rationale behind coalescent samplers settings would be to try to fill in the incomplete template files that are in the "input/template" folder based on the models in figure 1.**

To run the first analysis, submit the job script "sim_for_inference.run" with the command `sbatch`. This script launches 5 independent runs of each model (which will be used to find the best run, i.e. the run with the highest likelihood).

*sim_for_inference.run* content:

```bash
#!/bin/bash

##Script to simulate samples for statistical inference

#SBATCH -p normal                
#SBATCH -c 8                     
#SBATCH --mem=6GB
#SBATCH --job-name fsc2709-ex01                
#SBATCH -o %j.out              
#SBATCH -e %j.err

# modules
module load fastsimcoal2

# setting variables and creating folders
INPUT=../input_files/SFS_for_inference
cd $INPUT

# running fastsimcoal2
for i in {1..5}
do
	mkdir run0$i
	cp ppgcourse_model0.tpl ppgcourse_model0.est ppgcourse_model0_DSFS.obs run0$i
	cd run0$i
	fsc2709 -t ppgcourse_model0.tpl -e ppgcourse_model0.est -d -C8 -B8 -n20000 -L20 -M -q -u
	cd ..
done

for i in {1..5}
do
	mkdir run1$i
	cp ppgcourse_model1.tpl ppgcourse_model1.est ppgcourse_model1_DSFS.obs run1$i
	cd run1$i
	fsc2709 -t ppgcourse_model1.tpl -e ppgcourse_model1.est -d -C8 -B8 -n20000 -L20 -M -q -u
	cd ..
done
```

> -t: name of the template parameter file  
> -e: name of the prior definition file  
> -d: runs fastsimcoal2 with DAF 
> -C and -B: are multithreading parameters  
> -n20000: 20000 simulations to estimate de SFS in each cycle    
> -L20: 20 optimization cycles   
> -M: indicates parameter estimation  
> -q: quite run  
> -u: uses the multidimensional SFS  

### Finding the best model:

To find the best model, you should select the runs with the highest likelihood for each model (MaxEstLhood) and use them to compute the AIC. Here you have some shell commands to facilitate you this task. **Still form the "scripts"" folder** type:

model 0:

```bash
 grep -v "MaxEstLhood" ../input_files/SFS_for_inference/run0*/ppgcourse_model0/ppgcourse_model0.bestlhoods | sort -k 3 | head -1 | awk '{print "For model 0, replicate "NR " has the highest MaxEstLhood:",$3}'
```
model 1:

```bash
 grep -v "MaxEstLhood" ../input_files/SFS_for_inference/run1*/ppgcourse_model1/ppgcourse_model1.bestlhoods | sort -k 3 | head -1 | awk '{print "For model 1, replicate "NR " has the highest MaxEstLhood:",$6}'
```
You can now calculate the AIC for each model.

*AIC = 2k - 2lnL*

> Remember that k is the number of estimated parameters. **Notice that `fastsimcoal2` "MaxEstLhood" values are in log<sub>10</sub> units!..you need to convert them to ln before computing AIC**
</br>

<details><summary>CLICK ME IF YOU HAVE FINISHED THE EXAMPLE 1</summary>

### QUESTIONS

<p>

```ruby
1. Which is the best model based on AIC?  
2. Which are the maximum likelihood estimates of model parameters under this best model? HINT: Check the file "ppgcourse_model*_maxL.par" of the best run from the best model
```

</p>

</details>

## Example 2: Simulations for validating statistical methods

A second potential application of computer simulations under the coalescent is the validation of statistical methods used for demographic inference. Using simulation programs, you can generate **pseudo-observed data sets/samples (PODS)** that can be used to check the performance of different approaches and software (since you know the real model and model parameters that generated the data). You can for instance simulate genomic sequence data under the neutral model to test the properties of a particular set of summary statistics for detecting deviations of neutrality, or simulate genome-wide SNP data to study the performance of the composite likelihood approximation used in programs such as *fastsimcoal2* (example 1) or the one of software than you will use tomorrow (momi2) based on, for instance, simulated regions with partially linked SNPs and variable mutation rates. 

Here you will use *fastsimcoal2* to **simulate SNP data** under the best model (and model parameters) from example 1.  
To do that, you should copy the file with the maximum likelihood estimates under the best model from previous runs ("input_files/SFS_for_inference/run*/ppgcourse_model*/ppgcourse_model*_maxL.par") into the folder "/input_files/SNP_for_validation". **If you prefer, you can use the input file obtained from a prevous run and that is already available in the last folder "ppgcourse_bestmodel_maxL.par".**

> Notice that since you have done very few runs and optimization cycles to maximize likelihoods in the example runs, each of you may have different ML estimates of model parameters under the best model (due to local maxima).

Here you have the chromosome structure and linkage block configuration that you will use for SNP data simulation (**if you use your own "*_maxL.par" file from example 1, you must replace the last part of the file with these lines**). You will simulate 5 independent genomic regions of different sizes and with different recombination rates (within each region!...regions are independent) and mutations rates combinations:

(last part of the file "ppgcourse_bestmodel_maxL.par"):

```
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
5
//per Block: data type, num loci, rec. rate + optional parameters
DNA 10000 0.000000016 0.000000016
DNA 30000 0.0000002 0.0000002
DNA 20000 0.000000006 0.0000001
DNA 10000 0.0000001 0.00000005
DNA 30000 0.000000000 0.000000016
```
> We are now specifying DNA as the marker to simulate

So, the job script to run this example 2 looks like this:

```bash
#!/bin/bash

##Script to simulate samples for validation

#SBATCH -p normal                
#SBATCH -c 8                     
#SBATCH --mem=6GB
#SBATCH --job-name fsc2709-ex02                
#SBATCH -o %j.out              
#SBATCH -e %j.err

# modules
module load fastsimcoal2

# setting variables and creating folders
INPUT=../input_files/SNP_for_validation
cd $INPUT

# running fastsimcoal2
fsc2709 -i ppgcourse_bestmodel_maxL.par -n100 -s0 -C8 -B8
```

> -i: name of the parameter file   
> -n100: 1000 simulations (replicates)  
> -C and -B are multithreading parameters 
> -L20: 20 optimization cycles   
> -s0: output SNP data (although DNA is specified as a marker in the file “ppgcourse.par”, this flag causes the program to print only SNP) 

This script launches a run that generates SNP data (100 replicates) under the specified demographic model and the chromosomal and block configuration set in the “par” file. 

**NOTICE** that in this exercise you have simulated an example with SNP data under a given demographic model and chromosome configuration. In order to use these data to validate software like *fastimcoal2* or *momi2*, you should calculate the SFS of each of these simulated SNPs, using for instance the software [Arlequin](http://cmpg.unibe.ch/software/arlequin35). 

You can use these data (and any other simulation scheme) to run and validate an approximate Bayesian computation ([ABC]( https://www.annualreviews.org/doi/abs/10.1146/annurev-statistics-030718-105212)) approach. In this case, you should apply PODS the same summary statistics and transformations made to the simulated and observed data in the ABC analysis. 

## Example 3: Simulations for predictive uses

Finally, you might also use computer simulations for **predictive purposes and general understanding of genome evolution**. In this third example, you will investigate the effect that an introgression event from a ghost population has on the history of a group of sampled populations. Newly, you will use *fastsimcoal2* to simulate genomic data under the two models in [Figure1](figure1.pdf). In this example, however, you will simulate DNA data, specifically a region of 5Kb, which will be *sampled ONLY in Europe, Asia and Africa*. Here, we will consider that this region was **introgressed from Neanderthals to Eurasians and maintained in extant populations**. 

In the folder "/input_files/DNA_for_predictive" you have the input files to run this example.
Here you can see the chromosome structure and linkage block configuration set in these file:

```
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
DNA 5000 0 0.00000016
```

..and here you have the content of job script to run this example 3:

```bash
#!/bin/bash

###Script to simulate samples for predictive uses

#SBATCH -p normal                
#SBATCH -c 8                     
#SBATCH --mem=6GB
#SBATCH --job-name fsc2-ex03                
#SBATCH -o %j.out              
#SBATCH -e %j.err

# modules
module load fastsimcoal2

# setting variables and creating folders
INPUT=../input_files/DNA_for_predictive
cd $INPUT

# executing fastsimcoal2 for model0
fsc2709 -i ppgcourse_model0.par -n5 -c8 -B8

# executing fastsimcoal2 for model1
fsc2709 -i ppgcourse_model1.par -n5 -c8 -B8
```

This script launches a run that generates 5 replicates of a DNA fragment under the specified demographic models and chromosomal configuration set in the input filee. 

Now, you can for example estimate the phylogenetic trees with the alleles of the three sampled populations using the program IQTree. To do that, you need to convert the “\*arp” files into FASTA files and run `IQTree`. Use the script "iqtree.run" in the folder "scripts" to submit these jobs using `sbatch`.

The content of the "iqtree.run" script is:

```bash
#!/bin/bash

##Script to generate trees

#SBATCH -p normal                
#SBATCH -c 8                     
#SBATCH --mem=6GB
#SBATCH --job-name fsc2709-ex04                
#SBATCH -o %j.out              
#SBATCH -e %j.err

# modules
module load iq-tree

# iqtree model 0
cd ../input_files/DNA_for_predictive/ppgcourse_model0
FILES="*.arp"

for f in $FILES
do
        grep "^[0-9]" $f | perl -p -e  's/^/>/g' |  perl -p -e 's/\t[0-9]+//g' | perl -p -e 's/\t/\r/g' > model0.fas
        iqtree -s model0.fas -pre $f -nt 8 -redo
done


# iqtree model 1                                                                                                                                  
cd ../ppgcourse_model1
FILES="*.arp"

for f in $FILES
do
        grep "^[0-9]" $f | perl -p -e  's/^/>/g' |  perl -p -e 's/\t[0-9]+//g' | perl -p -e 's/\t/\r/g' > model1.fas
        iqtree -s model1.fas -pre $f -nt 8 -redo
done
```

To be able to visualize the generated trees (one for each model) you need to execute the program `Figtree` in your personal computers.<br/><br/>

</br>

<details><summary>CLICK ME IF YOU HAVE FINISHED THE EXAMPLE 3</summary>

### QUESTIONS

<p>

```ruby
1. How could the effect of introgression be made more evident?  
2. How should the simulations be changed? what parameters and to what extent? 
```

</p>

</details>

 
