# Divergence Time Estimation based on Molecular clocks
# Tutorial in a nutshell


Instructor: **Jesus Lozano-Fernandez**

We are going to implement a **NODE CALIBRATED** divergence time estimation and a **TOTAL EVIDENCE ANALYSIS** in the context of bear evolution; and to compare both outputs to infer the age of the MRCA of Ursinae, and assess the impact of these different methods, models and priors. For coherence with the Bayesian class, we are going to use also [**RevBayes**](https://revbayes.github.io), which is a flexible software with his own programming language similar to R, named Rev. Is out of the scope of this short tutorial to be able to fully understand the language, the models or the parameters. At the end of the class, you are expected to have a grasp of the main difference between methods, the utility of fossil and stratigraphic calibrations and the relevance of the prior assumptions. This tutorial is based on material available for learning divergence time estimatin using RevBayes, available at [https://revbayes.github.io/tutorials/dating/](https://revbayes.github.io/tutorials/dating/). All the scripting to run the analyss and the datasets has been already built and are in place, so you won't need to modify them, just running the analyses. Morevoer, it has been uploaded also the output of the analyses which were previously ran, on case you face some difficulties. We will, though, examine the inpt files and run the analyes, and learn about the interpretation of results (using the softwares [**Tracer**](http://tree.bio.ed.ac.uk/software/tracer/), and [**FigTree**](http://tree.bio.ed.ac.uk/software/figtree/)).

 ![question](img/MolecularDating.jpg)

***

## NODE DATING

* *Based on the 3rd exercise ‘Estimating speciation times using node dating’*
https://revbayes.github.io/tutorials/dating/nodedate
* Additionally, RevBayes uses a scripting language to call functions and provide arguments: the `Rev` language. It is similar to `R` but more heavily scripted.


## Background

In this exercise we will use information from the fossil record to calibrate the molecular substitution rate to absolute time using node dating. This approach involves assigning probability densities that incorporate temporal information from the fossil record to particular nodes in the tree.

For running this exercise, we should be inside the node dating folder and just should contain two folders: `data/` & `script/`.

In the **`data`** folder, you will find the following files:
* `data/bears_cytb.nex` # *an alignment in NEXUS format of 1,000 bp of cytochrome b sequences for 8 bear species*
* `data/bears_taxa.tsv` # *a tab-separated table listing every bear species (both
fossil and extant) and their occurrence age ranges. For extant taxa, the minimum age is 0.0 (i.e. the present).*

In the **`script`** folder, you will find the following files:
* `scripts/sub_GTRG.Rev` # *Sequence substitution model GTR + Γ*
* `scripts/clock_relaxed_lognormal.Rev` # *Relaxed clock model Uncorrelated*
* `scripts/tree_BD_nodedate.Rev` # *Birth Death tree prior*
* `scripts/MCMC_dating_ex3.Rev` # *Master Rev file (loads the data, the other scripts specifying models and monitors the output of the MCMC sampler)*
* `scripts/MCMC_dating_ex3_prior.Rev` # *Master Rev file under the priors (not data)*

## Justification for Models and Priors

*Calibrations*

The file `bears_taxa.tsv` contains information about the stratigraphic ranges for 20 bear species, including 12 extinct species. We’re not going to use all of the information from this file in this exercise, because the node dating approach to calibration limits the amount of data we can take advantage of, but we’ll use some of this information to constrain the age of two nodes We’re going to add two node calibrations: one on **the root** and one on the **internal node** for the clade Ursinae.
* The oldest first appearance of a crown group bear in our dataset is Ursus americanus at 1.84 Ma.
* We will take advantage of a previous estimate (∼49 Ma) for the age of caniforms, which is the clade containing bears and other “dog-like” mammals, from [dos Reis et al. 2012](https://royalsocietypublishing.org/doi/abs/10.1098/rspb.2012.0683). We will assume that the age of crown bears cannot be older than this.

*The clock model*

Remember the clock (or branch-rate) model describes how rates of substitution vary (or not) across the tree. We are going to use the uncorrelated exponential relaxed clock model. In this model rates for each branch will be drawn independently from an exponential distribution.

*The tree prior*

Since all of the taxa included in the analysis in this exercise are living species we’ll use a birth-death model that doesn’t incorporate the fossil recovery process, however, we will add calibration information from the fossil record to generate timetrees on a non- arbitrary timescale.

*The master Rev script*

We are also interested in the age of the most-recent-common ancestor (MRCA) of all living bears. To monitor the age of this node in our MCMC sample, we must use the clade function to identify the node.The analysis is set to run two chains for **10,000 generations** (10 to 30 minutes, really short analysis).

## Starting the analysis

* 1- On the terminal, and start RevBayes by typing ```rb```
RevBayes will open, then you will upload the master Rev script were all the needed information to perform the whole analysis is in there. Type the following and the analysis will automatically start:

```
source("scripts/MCMC_dating_ex3.Rev")
```

* 2- By the time you are waiting for the results of this analyses, start a second analysis under the prior. It is always useful to examine the output of your MCMC analysis in the absence of information from the sequence data (i.e. without calculating the likelihood that comes from the substitution model). Type on the terminal ```rb```

```
source("scripts/MCMC_dating_ex3_prior.Rev")
```

When both analyses are finished, let’s check at the newly generated **`ouput`** folder, which will contain:
* `output/bears_nodedate.log` # *log files to upload on Tracer of the prior*
* `bears_nodedate_run_1.log` # *log data generated by chain 1*
* `bears_nodedate_run_1.trees` # *trees generated by chain 1*
* `bears_nodedate_run_2.log` # *log data generated by chain 2*
* `bears_nodedate_run_2.trees` # *trees generated by chain 2*
* `bears_nodedate.mcc.tre` # *Maximum Clade Credibility tree*
* `bears_nodedate.trees` # *sampled trees*
* `bears_nodedate_prior.log` # *log files to upload on Tracer of the prior*
* `bears_nodedate_prior_run_1.log` # *log data generated by chain 1 (prior, no data)*
* `bears_nodedate_prior_run_1.trees` # *trees generated by chain 1 (prior, no data)*
* `bears_nodedate_prior_run_2.log` # *log data generated by chain 2 (prior, no data)*
* `bears_nodedate_prior_run_2.trees` # *trees generated by chain 2 (prior, no data)*
* `bears_nodedate_prior.mcc.tre` # *Maximum Clade Credibility tree under prior*
* `bears_nodedate_prior.trees` # *sampled trees under the prior*



* 3- Analysing results with [**Tracer**](http://tree.bio.ed.ac.uk/software/tracer/)
Evaluate and Summarize Your Results
In this section, we will evaluate the *mixing* and *convergence* of our MCMC simulation using the program **Tracer**.

 ![question](img/tracer_load_file.png)

Immediately upon loading your log files, you will see the list of **Trace Files** on the left-hand side (you can load multiple files). The bottom left section, called **Traces**, provides a list of every parameter in the log file, along with the mean and the effective sample size (ESS) for the posterior sample of that parameter. The ESS statistic provides a measure of the number of independent draws in our sample for a given parameter. This quantity will typically be much smaller than the number of generations of the chain. In Tracer, poor to fair values for the ESS will be coloured red and yellow. You will likely see a lot of red and yellow numbers because the MCMC runs in this exercise are too short to effectively sample the posterior distributions of most parameters. A much longer analysis is provided in the output directory.

Look through the various parameters and statistics in the list of **Traces**.
- **Are there any parameters that have really low ESS? Why do you think that might be?**

 ![question](img/tracer_load_file.png)
