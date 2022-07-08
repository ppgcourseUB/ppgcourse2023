# Variant calling from high-throughput (low-coverage) sequencing data

**Instructor**: Filipe G. Vieira

In this session you will learn how to do:
* calculate genotype likelihoods
* allele frequency estimation
* variant calling
* genotype calling

## Software description and hints

For most of the examples, we will use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by Thorfinn Korneliussen and Anders Albrechtsen (and many contributors) at the University of Copenhagen.
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

According to its website *ANGSD is a software for analyzing next generation sequencing data. The software can handle a number of different input types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes. This is especially useful for low and medium depth data. The software is written in C++ and has been used on large sample sizes. This program is not for manipulating BAM/CRAM files, but solely a tool to perform various kinds of analysis. We recommend the excellent program SAMtools for outputting and modifying bamfiles.*

## Data description and access

For the exercises we are going to use mapped data (BAM files) from chr11 of the 1000 Genomes Consortium.

**NOTICE that the data to be used in this practice is already available on the server. We have downloaded them to a data folder (/data/datasets/Calling_variants_from_low-coverage_ngs_data/) before the session, to avoid excessive download traffic during class. YOU DON'T HAVE TO COPY ANYTHING IN  YOUR HOME, environmental variables are set to read the data in that folder.**

## Intro to ANGSD

First, we will learn **how to build a command line in ANGSD**, with the specific example of calculating genotype likelihoods.

To see a full list of options in ANGSD type:
```
angsd
```
and you should see something like
```
Overview of methods:
        -GL             Estimate genotype likelihoods
        -doCounts       Calculate various counts statistics
        -doAsso         Perform association study
        -doMaf          Estimate allele frequencies
        -doError        Estimate the type specific error rates
        -doAncError     Estimate the errorrate based on perfect fastas
        -HWE_pval               Est inbreedning per site or use as filter
        -doGeno         Call genotypes
        -doFasta        Generate a fasta for a BAM file
        -doAbbababa     Perform an ABBA-BABA test
        -sites          Analyse specific sites (can force major/minor)
        -doSaf          Estimate the SFS and/or neutrality tests genotype calling
        -doHetPlas      Estimate hetplasmy by calculating a pooled haploid frequency

        Below are options that can be usefull
        -bam            Options relating to bam reading
        -doMajorMinor   Infer the major/minor using different approaches
        -ref/-anc       Read reference or ancestral genome
        -doSNPstat      Calculate various SNPstat
        -cigstat        Printout CIGAR stat across readlength
        many others

For information of specific options type:
        ./angsd METHODNAME eg
                ./angsd -GL
                ./angsd -doMaf
                ./angsd -doAsso etc
                ./angsd sites for information about indexing -sites files
Examples:
        Estimate MAF for bam files in 'list'
                './angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1'
```

If you specify an option, it will print more options:
```
$ANGSD/angsd -bam
parseArgs_bambi.cpp: bam reader:
        -bam/-b         (null)  (list of BAM/CRAM files)
        -i              (null)  (Single BAM/CRAM file)
        -r              (null)  Supply a single region in commandline (see examples below)
        -rf             (null)  Supply multiple regions in a file (see examples below)
        -remove_bads    1       Discard 'bad' reads, (flag >=256)
        -uniqueOnly     0       Discards reads that doesnt map uniquely
        -show           0       Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
        -minMapQ        0       Discard reads with mapping quality below
        -minQ           13      Discard bases with base quality below
        -trim           0       Number of based to discard at both ends of the reads
        -trim           0       Number of based to discard at 5 ends of the reads
        -trim           0       Number of based to discard at 3 ends of the reads
        -only_proper_pairs 1    Only use reads where the mate could be mapped
        -C              0       adjust mapQ for excessive mismatches (as SAMtools), supply -ref
        -baq            0       adjust qscores around indels (as SAMtools), supply -ref
        -checkBamHeaders 1      Exit if difference in BAM headers
        -doCheck        1       Keep going even if datafile is not suffixed with .bam/.cram
        -downSample     0.000000        Downsample to the fraction of original data
        -nReads         50      Number of reads to pop from each BAM/CRAMs
        -minChunkSize   250     Minimum size of chunk sent to analyses

Examples for region specification:
                chr:            Use entire chromosome: chr
                chr:start-      Use region from start to end of chr
                chr:-stop       Use region from beginning of chromosome: chr to stop
                chr:start-stop  Use region from start to stop from chromosome: chr
                chr:site        Use single site on chromosome: chr
```

### General Guidelines
- As a general guidance, `-GL 1`, `-doMaf 1/2` and `-doMajorMinor 1` should be the preferred choice when data uncertainty is high.
- Detecting variable sites based on their probability of being SNPs is generally a better choice than defining a threshold on the allele frequency.
- Every dataset is different, so dedicated cutoffs and filtering should be perform to assess robustness of your SNPs
- Not a lot of tools support genotype likelihoods, but there are some (e.g. [ngsTools](https://github.com/mfumagalli/ngsTools))

## Exercise 1 - Data QC

Even though GL are specially suitable for analyses of low coverage data, it is always advisable to perform QC checks and (if necessary) perform mild filtering. ANGSD can perform some basic filtering of the data:
* quality and depth, see [here](http://www.popgen.dk/angsd/index.php/Filters)
* SNP quality, see [here](http://popgen.dk/angsd/index.php/SnpFilters)
* sites, see [here](http://popgen.dk/angsd/index.php/Sites)

The most common filters are:

Parameter | Meaning
--- | ---
-uniqueOnly 1 | only uniquely mapping reads
-remove_bads 1 | not tagged as bad
-only_proper_pairs 1 | only proper pairs
-C 50 | reduces the effect of reads with excessive mismatches
-baq 1 | base alignment quality ([explained here](http://samtools.sourceforge.net/mpileup.shtml)) used to rule out false SNPs close to INDELS
-minMapQ 20 | minimum mapping quality
-minQ 20 | minimum base quality
-minInd 5 | use only sites with data from at least N individuals
-setMinDepth 30 | minimum total depth
-setMaxDepth 150 | maximum total depth

ANGSD can extract some QC info. You can use this bacth file to send the job using `sbath`:

```
#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=angsd                                                                                         
#SBATCH --error angsd-%j.err                                                                                     
#SBATCH --output angsd-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8 

# directories and commands
DATA=/data/datasets/Calling_variants_from_low-coverage_ngs_data/
REF=$DATA/hs37d5.fa.gz
ls $DATA/*.bam | sort -t "." -k 5,5 > samples.bam_list

# module load                                                                                                           
module load angsd   

# running the program
angsd -b samples.bam_list -ref $REF -r 11:21000000-22000000 -out 1000G_QC -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minQ 0 -doCounts 1 -doQsDist 1 -doDepth 1 -maxDepth 500 -nThreads 8
```

and then we can plot them:
```
Rscript Scripts/plotQC.R 1000G_QC
```

Remember that to visualize the plot you must download the file to your computer and open it locally.

>**QUESTION**
Should we use all samples? And which values would you choose as sensible thresholds on quality score and depth (minimum and maximum)?


## Exercise 3 - Genotype likelihoods

Now we are ready to calculate the ***genotype likelihoods*** for each site at each individual.

To do so you need to specify which genotype likelihood model to use.
```
$ANGSD/angsd -GL
-GL=0:
        1: SAMtools
        2: GATK
        3: SOAPsnp
        4: SYK
        5: phys
        6: Super simple sample an allele type GL. (1.0,0.5,0.0)
        -trim           0               (zero means no trimming)
        -tmpdir         angsd_tmpdir/   (used by SOAPsnp)
        -errors         (null)          (used by SYK)
        -minInd         0               (0 indicates no filtering)

Filedumping:
        -doGlf  0
        1: binary glf (10 log likes)    .glf.gz
        2: beagle likelihood file       .beagle.gz
        3: binary 3 times likelihood    .glf.gz
        4: text version (10 log likes)  .glf.gz
```
A description of these different implementation can be found [here](http://www.popgen.dk/angsd/index.php/Genotype_likelihoods).
The GATK model refers to the first GATK paper, SAMtools is somehow more sophisticated (non-independence of errors), SOAPsnp requires a reference sequence for recalibration of quality scores, SYK is error-type specific.
For most applications and data, GATK and SAMtools models should give similar results.

>**QUESTION**
How many genotype likelihoods do we expect per position? Why?

To calculate GL, we can modify our batch script by replacing the `angsd` command with:
```
angsd -b samples.bam_list -ref $REF -r 11:21000000-22000000 -out 1000G_GL -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -doCounts 1  -setMinDepth 30 -setMaxDepth 150 -GL 2 -doGlf 4 -nThreads 8
```
Parameter | Meaning
--- | ---
-GL 2 | genotype likelihood model as in GATK
-doGlf 2 | output in BEAGLE format

>**QUESTION**
What output files were created?
What information do they have?


## Exercise 3 - Estimation of allele frequencies

We now want to estimate allele frequencies at each site without relying on genotype calls.
In other words, at each site we want to to estimate (or count) how many copies of different alleles (two in case of biallelic variants) we observe in our sample (across all sequenced individuals).
However with low depth data direct counting of individually assigned genotypes can lead to biased allele frequencies.

ANGSD has an option to estimate **allele frequencies** taking into account data uncertainty from genotype likelihoods:
```
$ANGSD/angsd -doMaf
-doMaf  0 (Calculate persite frequencies '.mafs.gz')
        1: Frequency (fixed major and minor)
        2: Frequency (fixed major unknown minor)
        4: Frequency from genotype probabilities
        8: AlleleCounts based method (known major minor)
        NB. Filedumping is supressed if value is negative
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
        3: Using SFS as prior (still in development)
        4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an
Filters:
        -minMaf         -1.000000       (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
        -rmTriallelic   0.000000        (Remove sites with a pvalue lower)
Extras:
        -ref    (null)  (Filename for fasta reference)
        -anc    (null)  (Filename for fasta ancestral)
        -eps    0.001000 [Only used for -doMaf &8]
        -beagleProb     0 (Dump beagle style postprobs)
        -indFname       (null) (file containing individual inbreedcoeficients)
        -underFlowProtect       0 (file containing individual inbreedcoeficients)
NB These frequency estimators requires major/minor -doMajorMinor
```

Since the estimation of allele frequencies requires the specification of how to assign the major and minor alleles (if biallelic):
```
$ANGSD/angsd -doMajorMinor
        -doMajorMinor   0
        1: Infer major and minor from GL
        2: Infer major and minor from allele counts
        3: use major and minor from a file (requires -sites file.txt)
        4: Use reference allele as major (requires -ref)
        5: Use ancestral allele as major (requires -anc)
        -rmTrans: remove transitions 0
        -skipTriallelic 0
```

A possible command line to estimate allele frequencies might be (replace the angsd command in your batch fiile accordingly):
```
angsd -b samples.bam_list -ref $REF -r 11:21000000-22000000 -out 1000G_GL_MAF -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -doCounts 1  -setMinDepth 30 -setMaxDepth 150 -GL 2 -doMajorMinor 1 -doMaf 1 -nThreads 8
```
Parameter | Meaning
--- | ---
-doMajorMinor 1 | How to define major and minor alleles (if making a VCF, consider 4)
-doMaf 1 | How to calculate minor allele frequency

>**QUESTION**
What output files were created?
What information do they have?

```
chromo  position        major   minor   ref     knownEM nInd
11      21013539        A       C       A       0.000004        14
11      21013540        A       C       A       0.000004        13
11      21013541        A       C       A       0.000004        12
11      21013542        A       C       A       0.000003        12
11      21013543        A       C       A       0.000004        12
11      21013544        A       C       A       0.000006        12
11      21013545        A       C       A       0.000002        11
11      21013546        A       C       A       0.000003        12
11      21013547        A       C       A       0.000005        11
11      21013548        A       C       A       0.000005        12
11      21013549        A       C       A       0.000003        13
11      21013550        A       C       A       0.000001        13
11      21013551        A       C       A       0.000001        13
11      21014262        T       A       T       0.000001        11
11      21014263        T       A       T       0.000001        11
11      21014264        T       A       T       0.000001        11
11      21014265        T       A       T       0.000004        8

```
where `knownEM` specifies the algorithm used to estimate the (minor) allele frequency which is given under that column. The columns are: `chromosome`, `position`, `major allele`, `minor allele`, `reference allele`, `allele frequency`, `p-value for SNP calling` (if -SNP-pval was called), `number of individuals with data`.

You can notice that many sites have low allele frequency, probably reflecting the fact that that site is monomorphic. So, we may be interested in looking at allele frequencies only for sites that are actually variable in our sample or, in other words, perform **SNP calling**.


## Exercise 4 - SNP calling

There are two main ways to call SNPs using ANGSD:
```
        -minMaf         0.000000        (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
```
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

>**QUESTION**
Try calling SNPs and record how many sites are predicted to be variable for several scenario (e.g. vary cut-off and method).
Which sites are discordant between scenarios? What are their estimated allele frequencies?
Which frequencies are more difficult to estimate and therefore affect SNP calling?

SNP p-value | Number sites
--- | ---
0.05 | 8406
0.01 | 6546
0.001 | 5584
1e-4 | 5319
1e-6 | 5047


## Exercise 5 - Genotype calling

For low coverage samples, it is advisable to work always with genotype likelihoods but, sometimes, that is not always possible (most programs do not support them). For that reason, ANGSD can also call genotypes. To improve genotype calling, ANGSD uses the GL to calculate genotypes probabilities for each site on each individual, calling the genotype with the highest probability.

The option to call genotypes is `-doGeno`:
```
$ANGSD/angsd -doGeno
-doGeno 0
        1: write major and minor
        2: write the called genotype encoded as -1,0,1,2, -1=not called
        4: write the called genotype directly: eg AA,AC etc
        8: write the posterior probability of all possible genotypes
        16: write the posterior probability of called genotype
        32: write the posterior probabilities of the 3 gentypes as binary
        -> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
        -postCutoff=0.333333 (Only genotype to missing if below this threshold)
        -geno_minDepth=-1       (-1 indicates no cutof)
        -geno_maxDepth=-1       (-1 indicates no cutof)
        -geno_minMM=-1.000000   (minimum fraction af major-minor bases)
        -minInd=0       (only keep sites if you call genotypes from this number of individuals)

        NB When writing the posterior the -postCutoff is not used
        NB geno_minDepth requires -doCounts
        NB geno_maxDepth requires -doCounts
```

Therefore, if we set `-doGeno 2`, genotypes will be coded as `-1,0,1,2`(missing data, and number of alternate alleles). If we want to print the major and minor alleles as well then we set `-doGeno 3`.

To calculate the posterior probability of genotypes we need to define a model.
```
$ANGSD/angsd -doPost
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
        3: Using SFS as prior (still in development)
        4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an
```

If you have a population, you can use `-doPost 1` (frequencies under HWE as prior) but here, since we have individuals from different populations, we will use `-doPost 2` (uniform prior). Replace the angsd command in your batch script with this one):
```
angsd -b samples.bam_list -ref $REF -r 11:21000000-22000000 -out 1000G_CG_MAF -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -doCounts 1  -setMinDepth 30 -setMaxDepth 150 -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doPost 2 -doGeno 3 -doBcf 1 -nThreads 8
```
Parameter | Meaning
--- | ---
-doPost 2 | how to calculate genotype posterior probabilities
-doGeno 3 | format of output genotype file
-doBcf 1 | output BCF file

>**QUESTION**
What output files were created?
What information do they have?
How many sites have at least one missing genotype?

You can control how to set missing genotype when their confidence is low with `-postCutoff`.
For instance, we can set as missing genotypes when their (highest) genotype posterior probability is below 0.95, in this case the line to replace with in the script would be:

```
angsd -b samples.bam_list -ref $REF -r 11:21000000-22000000 -out 1000G_CG_MAF -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -doCounts 1  -setMinDepth 30 -setMaxDepth 150 -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doGeno 3 -doPost 2 -postCutoff 0.95 -nThreads 8
```

Setting this threshold depends on the mean sequencing depth of your data, as well as your application.
For some analyses you need to work only with high quality genotypes (e.g. measure of proportion of shared SNPs for gene flow estimate), while for others you can be more relaxed (e.g. estimate of overall nucleotide diversity).


## Exta Exercises - PCA and ngsAdmix
One of the most common exploratory analyses done, is probably a PCA plot and, right after, probably an admixture plot. ANGSD (and some extra tools) can do some of these analyses. Use the same template batch script and replace the angsd command with:
```
angsd -b samples.bam_list -ref $REF -r 11:21000000-22000000 -out 1000G_GL_PCA -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -doCounts 1  -setMinDepth 30 -setMaxDepth 150 -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doIBS 1 -doCov 1 -makeMatrix 1 -minMaf 0.05 -nThreads 8
```
Parameter | Meaning
--- | ---
-doGlf 2 | get GL in beagle format
-doIBS 1 | do IBS analysis
-doCov 1 | output the covariance matrix (for PCA)
-makeMatrix 1 | output the pairwise IBS matrix (distance between pairs of individuals; for MDS)
-minMaf 0.05 | remove sites with frequency lower than 0.05

**PCA**
```
R --vanilla --slave -e 'm <- as.matrix(read.table("1000G_GL_PCA.ibsMat")); mds <- cmdscale(as.dist(m)); plot(mds,lwd=2,ylab="Dist",xlab="Dist",main="multidimensional scaling",col=rep(1:3,each=5)); m <- as.matrix(read.table("1000G_GL_PCA.covMat")); e <- eigen(m); plot(e$vectors[,1:2],lwd=2,ylab="PC 2",xlab="PC 1",main="Principal components",col=rep(1:3,each=5),pch=16)'
mv Rplots.pdf 1000G_GL_PCA.pdf
```
Remember that to visualize this pdf file you must download the file to your computer and open it locally.

**ngsAdmix**

Here you have a template batch script to run `MGSadmix` in the cluster:

```
#!/bin/bash                                                                                                             

# define names                                                                                                          
#SBATCH --job-name=ngsadmix                                                                                        
#SBATCH --error ngsadmix-%j.err                                                                                     
#SBATCH --output ngsadmix-%j.out                                                                                    

# memory and CPUs request                                                                                               
#SBATCH --mem=6G                                                                                                        
#SBATCH --cpus-per-task=8   

# module load                                                                                                           
module load angsd         

# running the program
NGSadmix -likes 1000G_GL_PCA.beagle.gz -K 3 -P 4 -o 1000G_GL -minMaf 0.05 -P 8
```

```
R --vanilla --slave -e 'admix <- t(as.matrix(read.table("1000G_GL.qopt"))); barplot(admix,col=1:3,space=0.01,border=NA,xlab="Individuals",ylab="admixture")'
mv Rplots.pdf 1000G_GL_PCA.admix.pdf
```
Again, to visualize this pdf file you must download the file to your computer and open it locally.

>**QUESTION**
What output files were created?
What information do they have?
What would happen if we:
- Subsample the BAM files (`samtools view -s 12345.3`) and re-run analyses, both with genotype likelihoods and called genotypes
- What SNPs would be different? What are their frequencies?
- How do the PCA and admixture plots differ (use `admixture` instead `ngsAdmix`)?
