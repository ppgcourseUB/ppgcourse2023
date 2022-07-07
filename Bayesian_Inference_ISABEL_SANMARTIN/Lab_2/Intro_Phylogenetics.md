# Lab 2: Phylogenetic Inference in RevBayes

We are going to program a Bayesian inference phylogenetic analysis in RevBayes using the `Rev` language. The three main components of a phylogenetic model (that we've learned about so far) are the *tree topology*, the *branch lengths*, and the *substitution model*. We will program them in this order.

## Data description and access

We are going to reconstruct the phylogeny of genus *Bufo* using a molecular dataset with two loci. Move into the `Lab_2` folder, where you can find all scripts and files needed for the practice.  Inside this folder, create a second folder and name it `data`. Copy the *cytB.nex* and *16S.nex* files containing the sequences for these two mitochondrial markers from the pggcourse website into the `lab_2/data` folder. 

```
mkdir data
cp -p cytb.nex data/
cp -p 16s.nex data/
```

## Launching RevBayes
Launch RevBayes by typing `rb-mpi` in the cluster command line (or `./rb` if your are using the Terminal). This should launch RevBayes and give you a command prompt (the `>` character); this means RevBayes is waiting for input.

## Constructing the phylogenetic model *interactively*

### Reading the sequence data

We have special functions that read alignments in `.nexus` or `.fasta` formats (among others). We read the alignments of the _cytb_ and _16s_ genes using the `readDiscreteCharacterData` function:

```
# read in cytb and 16S data
data_cytb = readDiscreteCharacterData("data/cytb.nex")
data_16s  = readDiscreteCharacterData("data/16s.nex")
```

We *concatenate* the two datasets into a single alignment. We also extract some useful information from the alignment that we will need later.
```
# concatenate them into one alignment
data = concatenate(data_cytb, data_16s)

# get some useful information about the data
taxa = data.taxa()
num_taxa = data.ntaxa()
num_branches = 2 * num_taxa - 3
num_sites = data.nchar()
```


### Creating a vector of movement proposals for `MCMC`
REMEMBER that in Bayesian Inference, a simulation approach called *Markov Chain Monte Carlo* is used to approximate the posterior probability of the phylogenetic model, as we cannot estimate this quantity analytically. For this, a composite *continuous-time Markov Chain* process (`CTMC`) with as many variables as parameters in the model is constructed, which has as its stationary distribution the posterior probability we wish to estimate. Using *Monte Carlo* simulations coupled with the *Metropolis-Hastings* algorithm, we update starting random values for each parameter in such a way that we approximate the stationary distribution by maximizing the likelihood in each step of the `CTMC` (climbing the "mountain") but occasionally allowing for a decrease in the (ln) Lk value. 

To perform this *updating*, we need movement proposals for each parameter, whose type depend on the parameter itself, e.g., *sliding window* for uniform parameters, *scale values* for rate parameters. Each movement proposal is accompanied by a *weight*, a constant variable with an assigned fixed value, which determines how often the parameter is updated during the MCMC analysis: the higher the weight, the more frequent the movement.

First, we create an empty vector, called `move_index`, which we will populate with the different *movement proposals* we need for each parameter in the model. The alternative option is to define a *move* after we define a parameter, but using this vector allows us to keep track of how many moves we have to create.

```
# make a move index variable
move_index = 0
```

### Specifying the tree topology
We use a uniform prior on the tree topology. This implies that all topologies have the same prior probability. To update the topology as we run the MC chain, we use two type of movements, *Nearest-Neighbor Interchange* (NNI) and *Subtree Pruning and Regrafting* (SPR).

We assign a large *weight* to this movement because the topology is a difficult parameter to estimate.

```
# We assume a uniform prior on topology.
topology ~ dnUniformTopology(taxa)
moves[++move_index] = mvNNI(topology, weight=10.0)
moves[++move_index] = mvSPR(topology, weight=10.0)
```

Now, you can see why we need an empty movement vector. The command `++move_index` increases the value of `move_index` each time we use it. This makes it easy to add moves to the vector without remembering for ourselves how many moves we have included in our analysis.

### Specifying the branch lengths
Now, we create a vector of branch lengths `br_lens`, one per branch in the phylogram, using a `for` loop. Each `br_lens[i]` is a stochastic variable, whose value is drawn from an exponential prior distribution with rate parameter 10 (this assigns a higher prior to shorter branch lengths). We specify a *scale* move on each `br_lens[i]`; this move multiplies or divides the parameter by a random number. We also keep track of the *tree length* (`TL`), a deterministic variable which is the sum of the branch lengths.
```
for(i in 1:num_branches){
  br_lens[i] ~ dnExponential(10.0)
  moves[++move_index] = mvScale(br_lens[i], weight=1.0)
}
TL := sum(br_lens)
```

### Creating the phylogram
A "phylogram" is a tree with branch lengths measured in the expected number of substitutions. We create a deterministic variable `phylogeny` to wrap the topology and branch lengths together into a phylogram using the function `treeAssembly`. Notice that because this is a deterministic node, we do not need to define movement proposals, as we did for the stochastic variables `br_lens` and `topology`.

```
phylogeny := treeAssembly(topology, br_lens)
```
### Specifying the substitution model
The last component of the model is the substitution model (the `Q matrix`), which comprises two different types of parameters: the _exchangeability rates_ (`er`) and the _stationary frequencies_ (`pi`). In this analysis, we will use the _GTR_ substitution model, in which there are six rates and four frequencies, 10 stochastic parameters. Notice that in a _GTI_ (time-irreversible) model, there would four frequencies but 12 rates. 

We assign a prior *Dirichlet* distribution to the stationary frequencies, a multinomial distribution in which each frequency has its own prior but the sum of priors equals to 1. Below, all priors are made equally probable. We assign the same type of Dirichlet prior to the rate parameters. Each of these priors is updated simultaneously in the same movement proposal.

Finally, we construct a deterministic variable, the `Q matrix` by wrapping together the `er` and `pi` parameters using the function `fnGTR`. 

```
pi ~ dnDirichlet(v(1,1,1,1))
moves[++move_index] = mvDirichletSimplex(pi, weight=1.0)

er ~ dnDirichlet(v(1,1,1,1,1,1))
moves[++move_index] = mvDirichletSimplex(er, weight=1.0)

Q := fnGTR(er, pi)
```

### Specifying the likelihood function for the `CTMC` process

We construct a stochastic node, the `seq` variable, and assign it a prior that is updated using the `dnPhyloCTMC` distribution; this function takes as arguments the *phylogeny* and *Q matrix* variable values. We need to specify the type of data (*DNA*) because a CTMC can also be used with morphological data.

We them "clamp" the seq variable to the observed data, i.e., the alignment. This tells RevBayes to treat the random variable as an observation. OBS: If we want to sample from the prior, we need to comment out (#) the last line.

```
seq ~ dnPhyloCTMC(tree=phylogeny, Q=Q, type="DNA")
seq.clamp(data) # attach the observed data
```

### Running the analysis
To run the MCMC analysis, we wrap up the model object `my_model`, as dependent on stochastic variable *Q*.
```
#########################
# Make the model object #
#########################

my_model = model(Q)
```
Define monitors that will print to the screen the results and also print to a file (the most important). There are two types of files: `bufo.log` will print the parameter values; `bufo.trees` will print the trees. We decide how often these parameters are printed with the `printgen` command.

```
######################
# Make the monitors. #
######################

monitors[1] = mnModel(filename="output/bufo.log",printgen=100, separator = TAB)
monitors[2] = mnFile(filename="output/bufo.trees",printgen=100, separator = TAB, phylogeny)
monitors[3] = mnScreen(printgen=100, TL)
```
Create an MCMC analysis object with the number of generations (number of steps in the CTMC) to run.

```
################
# Run the MCMC #
################

analysis = mcmc(my_model, monitors, moves)
analysis.run(generations=100000)
```

### Summarizing the tree
Finally, we summarize the marginal distribution of trees by building the MAP (_maximum a posteriori_) tree, which is the tree topology that has the highest posterior probability, and computing the posterior probabilities for each node in the MAP tree.

```
# start by reading in the tree trace
treetrace = readTreeTrace("output/bufo.trees", treetype="non-clock")
map_tree = mapTree(treetrace,"output/bufo_MAP.tree") 

# exit the program
q()
```

## Running the analysis in script mode (inside RevBayes)
If your want to run the analysis above without typing every command at the RevBayes prompt, you can source the script, as we learnt in *Lab_1*. Launch Source the `my_phylogenetic_analysis.Rev` script, which contains all commands above (OBS!! Remember the quotation marks!).

```
source("my_phylogenetic_analysis.Rev")
```

## Running the analysis in batch mode from the cluster
Finally, we can use a bash file if we are running the script within the cluster. Below is an example of a bash file to run the script `my_phylogenetic_analysis.Rev`. Because we are running in a cluster, we need to include a command for the output to be saved. You can find this script `myScript-cluster.sh` inside `Lab_2` folder.

```
#!/bin/bash                                                                                                             

#SBATCH -p normal                                                                                                       
#SBATCH -n 8                                                                                                            
#SBATCH -c 1                                                                                                            
#SBATCH --mem=6GB                                                                                                       
#SBATCH --job-name orthofinder-job01                                                                           \        
#SBATCH -o %j.out                                                                                                       
#SBATCH -e %j.err                                                                                                       

module load revbayes

mpirun -np 8 rb-mpi my_phylogenetic_analysis.Rev
```

## Visualizing the phylogeny and clade support

Once we quit RevBayes, we can visualize the output: the phylogeny and the clade statistical support (the *nodal posterior probabilities*), using the software `FigTree` (http://tree.bio.ed.ac.uk/software/figtree/). Open the file `bufo.trees` using `FigTree` and check clade support by clicking on *node labels: posterior* on the left pane.

## Evaluating stationarity and mixing of parameters

We can also evaluate how well our MCMC analysis approximated the posterior distribution of each parameter in the model by using the software `Tracer`  (https://beast.community/tracer). Download `Tracer` into your computer (it does not work in the cluster) and open it by double clicking: it is a graphical programme. Use *Import Trace File* from the Tracer menu to read the output file `bufo.log`. Check the shape of the posterior distributions, the ESS values (should be > 200), and search for potential correlations between parameters.
![This is an image] (../images/Tracer.png)


# Exercises

Repeat the above Bayesian phylogenetic analyses with the three alignments for the _Viburnum_ dataset located at the pggcourse and named `matK.nex`, `ITS.nex`, and `trnS-G.nex`. Inspect the posterior probabilities in FigTree. How do the support values compare to those you got for the _Bufo_ phylogeny?

How well did the above MCMC perform? Are the ESS values for all of the continuous parameters greater than 200?

Run your _Viburnum_ analysis again, or compare your tree with its posterior probabilities to those of your neighbor. Are they similar or different? If they are different, can you hypothesize why they might be different?

(Optional) Try setting up an analysis of the _Viburnum_ dataset where each sequence gets its own *Q matrix*, but all share the same topology and branch lengths. You can accomplish this by creating multiple *Q* matrices and `seq` objects, each corresponding to a different gene region. Can you express this model in graphical form?
