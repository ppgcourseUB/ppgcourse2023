# Process raw data

## 1. Obtaining sequence data

We will download the following CYTB sequences from the NCBI web server and save them in using the following format `genus_species_cytb.fasta`, where we will only use the first 3 letters of the genus and the first three letters of the species. Please follow the links below to download the corresponding sequences:

* [*Ailuropoda_melanoleuca*, NCBI RefSeq: NC_009492.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_009492.1?from=15529&to=16668&report=fasta).
* [*Tremarctos_ornatus*, NCBI RefSeq: NC_009969.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_009969.1?from=15030&to=16169&report=fasta).
* [*Helarctos_malayanus*, NCBI RefSeq: NC_009968.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_009968.1?from=15054&to=16193&report=fasta).
* [*Ursus_americanus*, NCBI RefSeq: NC_003426.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_003426.1?from=15113&to=16252&report=fasta).
* [*Ursus_thibetanus*, NCBI RefSeq: NC_008753.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_008753.1?from=15139&to=16278&report=fasta).
* [*Ursus_arctos*, NCBI RefSeq: NC_003427.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_003427.1?from=15306&to=16445&report=fasta).
* [*Ursus_maritimus*, NCBI RefSeq: NC_003428.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_003428.1?from=15301&to=16440&report=fasta).
* [*Melursus_ursinus*, NCBI RefSeq: NC_009970.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_009970.1?from=15086&to=16225&report=fasta).

## 2. Inferring alignment and phylogeny

### Processing raw input sequence data

First, we need to concatenate all the [sequences we have downloaded](ind_cytb) (see links provided above) in a unique FASTA file. In order to do this, we can run the following code from directory raw_data:

```sh
# Run from `raw_data` directory
for i in ind_cytb/*fasta
do
name=$( echo $i | sed 's/..*\///' | sed 's/\_cytb..*//' )
sed 's/^>..*/\>'${name}'/' $i >> ind_cytb/unaln_raw.fasta
done 
```

Now, we will generate a one-line FASTA file that is easier to parse:

```sh
# Run from `raw_data`
mkdir unaln_seq
../src/one_line_fasta.pl ind_cytb/unaln_raw.fasta
name1=$( echo ind_cytb/unaln_raw.fasta | sed 's/\.fasta//' )
name2=$( echo $name1 | sed 's/..*\///' )
mv $name1"_one_line.fa" unaln_seq/$name2".fasta"
```

### Alignment and tree inference

#### 1. Bayesian approach to co-infer alignment and tree (`BAli-Phy`)

We can run [`BAli-Phy`](https://doi.org/10.1093/bioinformatics/btab129) to co-infer phylogeny + alignment under a Bayesian framework given that the alignment is not too big. All the diagnostics and summary reports have been generated following the [`BAli-Phy` tutorials and recommendations](https://www.bali-phy.org/README.html#output). You can get started by following the next code snippet:

```sh
# Run from `raw_dat`
mkdir -p aln_seq/baliphy
cd aln_seq/baliphy
mkdir -p baliphy
cd baliphy
bali-phy ../unaln_raw.fasta -S hky85+Rates.gamma[5] -n cytb & # this will create a dir called `cytb_r1-1`
bali-phy ../unaln_raw.fasta -S hky85+Rates.gamma[5] -n cytb & # this will create a dir called `cytb_r1-2`

# Text written based on the BAli-Phy manual: https://www.bali-phy.org/README.html#mixing_and_convergence
# We can explore whether the chains may have converged by computing the ASDSF and the MSDFSF
#  - The SDSF value is the SD across runs of the posterior probabilities (PP) for that split.
#    By averaging the SDSF values across splits, we obtain the ASDSF value (Huelsenbeck and Ronquist, 2001).
#    An acceptable value of the ASDSF is < 0.01.
#  - The maximum of the SDSF values is the MSDSDF, which represents the range of variation in PP across 
#    the runs for the split with the most variation.
trees-bootstrap cytb-1/C1.trees cytb-2/C1.trees > partitions_bs.txt
##> Our results show the following:
##> ASDSF[min=0.100] = 0.007     MSDSF = 0.020
##>
##> In that way, we may say that there are chances our chains have reached convergence. Nevertheless, we will
##> run another diagnostic.
#
# We will calculate potential scale reduction factors (PSRF) to check that different runs have similar 
# posterior distributions. Only numerical variables may have a PSRF.
# The PSDRF is a ratio of the width of the pooled distribution to the average width of each distribution,
# and ideally should be 1. The PSRF is customarily considered to be small enough it is is less than 
# 1.01.
# The command below will report the following:
# - PSRF-80%CI: PSRF based on the length of 80% credible intervals (Brooks and Gelman 1998).
# - PSRF-RCF: for each individual distribution, the 80% credible interval is found. Then, the probability
#   of that interval (which may be more than 80%) is divided by the probability of the same interval under
#   the pooled distribution. 
# 
# This diagnostic can detect when a parameter value has stabilized at different values in several independent
# runs, indicating a lack of convergence. This situation might occur if different runs of the Markov chain were
# trapped in different modes and failed to adequately mix between modes.
statreport cytb-1/C1.log cytb-2/C1.log > report_PSRF.txt
##> Our results are the following:
##>  Ne  >= 240    (posterior)
##>  min burnin <= 248    (rs07:mean_length)
##>  PSRF-80%CI <= 1.007    (hky85:pi[C])
##>  PSRF-RCF <= 1.01    (hky85:pi[A])
##>
##> The PSRF is less than 1.01, and thus we can say that this check has been passed!
##>
##> It seems both diagnostics show that chances are our chains have reached convergence. Therefore, we will stop
##> our runs now.
```

Once you have checked for convergence and it seems it is time to stop them, you can find their PID number and kill the jobs (i.e., run from the terminal `kill -9 PID`). To summarise the results, we can run the following commands:

```sh
# We will generate the greedy consensus tree as the tre is fully resolved.
# We will run the following options:
# `-s 10%` Skip the first 10% of tree
# `--greedy-consensus`
#
# Run from `baliphy`
trees-consensus -s 10% cytb-1/C1.trees cytb-2/C1.trees --greedy-consensus=cytb_greedy_bp.tree
# We can also check topology convergence
trees-bootstrap cytb-1/C1.trees cytb-2/C1.trees > convergence_treetopo.txt

# Now, we can compute the alignment using posterior decoding
cut-range  cytb-1/C1.P1.fastas cytb-2/C1.P1.fastas --skip=100 | alignment-chop-internal --tree cytb_greedy_bp.tree | alignment-max > P1-max.fasta

# Last, we will write a summary output
# in HTML format
# Run also from `baliphy` dir
bp-analyze cytb-1 cytb-2
```

We can now open the file [`cytb_greedy_bp.tree`](aln_seq/baliphy/cytb_greedy_bp.tree) with `FigTree` to visualise the inferred tree topology. We will root the tree at *ail_mel* (i.e., *Ailuropoda_melanoleuca*). We will now select option `Ordering: increasing` after ticking the box `Trees> Order nodes`. In that way, we will be able to observe the tree topology without branch lengths, which shall match the next topology in Newick format:

```text
((((((urs_ame,urs_thi),hel_mal),(urs_arc,urs_mar)),mel_urs),tre_orn),ail_mel);
```

Now, we will export the rooted tree. First, we shall click `File> Expot Trees>`, then select `Newick` in the `Tree file format` box, tick the box `Save as currently displayed`, and then we shall save the rooted tree with branch lengths as `cytb_rooted_bl.tree` inside directory `baliphy`.

#### 2. Heuristics for alignment inference (`mafft`) and maximum likelihood for phylogeny inference (`RAxML-NG`)

First, given that we have a very small dataset, we can run either `mafft` or `muscle5` with the automatic commands. If you had a larger alignment, you may want to try other algorithms to speed up alignment inference. You can run both aligners and then compare the output:

```sh
# Run from `raw_data`
mkdir -p aln_seq/{mafft,muscle5}
cd aln_seq/mafft 
mafft --auto ../../unaln_seq/unaln_raw.fasta > cytb_aln_mafft.fasta
cd ../muscle5
muscle5.1 -align ../../unaln_seq/unaln_raw.fasta -output cytb_aln_muscle5.fasta
```

Given that both aligners have inferred the same alignment, then we can proceed with `RAxML-NG` to infer the best-scoring maximum-likelihood (ML) tre:

```sh
# Run from `raw_data`
mkdir tree
cd tree
raxml-ng --msa ../aln_seq/muscle5/cytb_aln_muscle5.fasta --model HKY+G5 --prefix cytb --threads 2 --seed 12345 --tree pars{25},rand{25}
```

We can now open the file [`cytb.raxml.bestTree`](tree/cytb.raxml.bestTree) with `FigTree` to visualise the inferred tree topology. We will root the tree at *ail_mel* (i.e., *Ailuropoda_melanoleuca*). As some branch lengths seem to have length 0, we will select option `Ordering: increasing` after ticking the box `Trees> Order nodes` and then option `Transform: cladogram` after ticking the box `Trees>Transform branches`. In that way, we will be able to observe the tree topology without branch lengths, which shall match the next topology in Newick format:

```text
((((((urs_thi,urs_ame),hel_mal),(urs_arc,urs_mar)),mel_urs),tre_orn),ail_mel);
```

Now, we will untick the box `Transform branches` and export the rooted tree. First, we shall click `File> Expot Trees>`, then select `Newick` in the `Tree file format` box, tick the box `Save as currently displayed`, and then we shall save the rooted tree with branch lengths as `cytb_rooted_bl.tree` inside directory `tree`.

----

You can see that we have inferred the same molecular alignment and tree topology regardless the approach chosen. Sometimes, however, that may not be the case and you will need to carry out further diagnostics. For this tutorial, however, we can now proceed to further format these input files for timetree inference. Please [click this link to continue with data formatting](../00_data_formatting/README.md)!
