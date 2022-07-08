

# Linkage disequilibrium based methods to detect positive selection

**Instructor:** Rocío Caro Consuegra



The goal of this practical session is to learn how to compute two linkage disequilibrium based tests designed to detect positive selection: the integrated haplotype score (Voight et al., 2006) and the cross-population extended haplotype homozygosity (Sabeti et al., 2007). 

We are going to focus on an already known signal that affects the *SLC24A5* (solute carrier family 24 member 5) gene, located in chromosome 15, and which plays a role in human skin pigmentation (Lamason et al., 2005). The ancestral allele of one of the reported coding polymorphism (rs1426654) is almost fixed in African and East Asian populations (93% and 99%, respectively), whereas in Europeans, the derived allele is fixed (100%). 

### Data description and access

We are performing our tests on publicly available sequencing data generated for the 1000 Genomes Project (1KGP)  (<http://www.internationalgenome.org/>). You will find a trimmed version of the chromosome 15 for an European (CEU; Utah residents with Northern and Western European ancestry) and an African population (YRI; Yoruba in Ibadan, Nigeria) in the data/ directory of this practical session. We have already performed quality control (QC) with PLINK/1.9b, with filters such as: missingness per genotype, missingness per individual, relatedness and minor allele frequency, etc. We have also phased the dataset using SHAPEIT and polarised the dataset using custom Python scripts.

### Selscan

There are multiple softwares that you can used to compute LD-based selection scans, but we are focusing on Selscan (Szpiech et al., 2014; https://github.com/szpiech/selscan).

#### iHS: integrated haplotype score

For the calculation of iHS at a site, selscan first calculates the integrated haplotype homozygosity (iHH) for the ancestral (0) and for the derived (1) haplotypes separately. The unstandardized iHS is the natural logarithm of iHH1/iHH0. This means that in loci with positive (negative) iHS values, the derived (ancestral) allele seems adaptive for a particular selective pressure. Note that Voight et al., 2006 defined it the other way around (iHH0/iHH1).

```bash
selscan --ihs --hap <hapfile> --map <mapfile> --out <outfile>
# or
selscan --ihs --vcf <vcffile> --map <mapfile> --out <outfile>
```

For standardization, iHS scores are normalized by frequency bins across the genome by subtracting the means and dividing it by the standard deviations within that bin. 

```bash
norm --ihs --files <file1.*.out> ... <fileN.*.out>
```

##### Input files

The required input files are:

- **hap**: Variant information coded as 0s and 1s where each row represents an haploid copy and each column the sorted loci separated by whitespace. This format can be obtained from a VCF file. For the latest versions of selscan, you can directly use the **VCF** file. In any case, remember that data has to be phased and polarised (coded as ancestral/derived instead of reference/alternative as usual).
- **map**: Contains genetic and physical map information. The first column contains the chromosome, the second the SNP id, the third the genetic position and the fourth the physical position. In case there is no information about the genetic position, use the physical position again in its place (if left empty, selscan won't work).
- **out**: Gives a base name for the output file.

##### Computing iHS

Once we have all the input files ready, we can send the job as follows:

```bash
#!/bin/bash
#SBATCH -J iHS_tutorial
#SBATCH --cpus-per-task=8
#SBATCH --mem=6GB
#SBATCH -o iHS_tutorial.out
#SBATCH -e iHS_tutorial.err

module load selscan/v1.2.0a

# Computing iHS for Utah residents (CEPH) with Northern and Western European ancestry
HAP_FILE=CEU.recode.vcf.hap
MAP_FILE=CEU_filled.map
IHS_FILE=CEU
selscan --ihs --hap $HAP_FILE --map $MAP_FILE --out $IHS_FILE --threads 8

# Computing iHS for Yoruba in Ibadan, Nigeria
HAP_FILE=YRI.recode.vcf.hap
MAP_FILE=YRI_filled.map
IHS_FILE=YRI
selscan --ihs --hap $HAP_FILE --map $MAP_FILE --out $IHS_FILE --threads 8

```

The option --threads is used to parallelize jobs and speed up the process. Check in your own servers the capacity of the system.

And for the normalization:

```bash
# Normalizing iHS for Utah residents (CEPH) with Northern and Western European ancestry
NORM_IHS_FILE=CEU.ihs.out
norm --ihs --files $NORM_IHS_FILE

# Normalizing iHS for Yoruba in Ibadan, Nigeria
NORM_IHS_FILE=YRI.ihs.out
norm --ihs --files $NORM_IHS_FILE
```

##### Output files

Selscan generates two output files:

- An **.out** file organized in the following columns:

  - locus ID
  - physical position
  - derived allele frequency 
  - Integrated haplotype homozygosity of the derived allele (iHH1)
  - Integrated haplotype homozygosity of the ancestral allele (iHH0)
  - unstandardized iHS

- A **.log** file containing the runtime parameters. It will also include information about the excluded loci. Loci might be excluded because of:

  - a minor allele frequency below 0.05
  - reaching a gap over 200kbp
  - reaching the chromosome edge before the EHH has decayed below 0.05

  These are the parameters by default, but they can be modified  using --maf, --max-gap, or --cutoff respectively.

When normalizing, we obtain two extra files:

- **[100bins].norm**, which is organized as the **.out** file but with two extra columns:
  - standardized iHS
  - 1 if |standardized iHS| >= 2, and 0 otherwise

- In case we have defined an *.err and an *.out file for the process when sending it to the server, the **.err** file will also contain information about the normalization. In our case, this means a four columns table:
  - bin, indicating the frequency bin
  - num, indicating the number of SNPs that fall within that bin
  - mean
  - variance

  By default, normalization will be computed by 100 frequency bins in the range [0, 1]. This parameter can also be modified using --bins. In case we would rather use windows of a constant bp size, see --bp-win. Be careful because this will mean a varying number of SNPs within each window.

#### XP-EHH: cross-population extended haplotype homozygosity

For the calculation of XP-EHH between two populations A and B, selscan computes iHH for each population independently. The unstandardized XP-EHH is the natural logarithm of iHHA/iHHB.  Loci with positive values of XP-EHH will indicate a signal of selection in population A, and loci with negative values in population B.

```bash
selscan --xpehh --hap <pop1_hapfile> --ref <pop2_hapfile> --map <mapfile> --out
<outfile>
# or
selscan --xpehh --vcf <pop1_vcffile> --vcf-ref <pop2_vcffile> --map <mapfile> --out
<outfile>
```

For standardization, XP-EHH scores are normalized by subtracting the genome-wide mean and dividing it by the standard deviation.

```bash
norm --xpehh --files <file1.*.out> ... <fileN.*.out>
```

##### Input files

For XP-EHH, the required input files are the same than in iHS plus a **--ref** file. This is a .hap formatted file (or .vcf) from the reference population (population B). 

Note that **--hap** and **--ref** must contain the same loci. In this case, the files do not need to be coded as ancestral/derived (as it happened for the computation of iHS), they just need to be consistent between them. Finally **--map** might be the .map file of any of the two populations, as long as they contain the information about the exact same loci.

##### Computing XP-EHH

Once we have all the input files ready, we can compute XP-EHH as:

```bash
#!/bin/bash
#SBATCH -J XPEHH_tutorial
#SBATCH --cpus-per-task=8
#SBATCH --mem=6GB
#SBATCH -o XPEHH_tutorial.out
#SBATCH -e XPEHH_tutorial.err

module load selscan/v1.2.0a

# XPEHH computation for YRI vs CEU
HAP1=YRI.recode.vcf.hap
HAP2=CEU.recode.vcf.hap
MAP=CEU_filled.map
XPEHH_FILE=YRIvsCEU
selscan --xpehh --hap $HAP1 --ref $HAP2 --map $MAP --out $XPEHH_FILE --threads 8
```

And for the normalization:

```bash
# Normalization
NORM_XPEHH_FILE=YRIvsCEU.xpehh.out
norm --xpehh --files $NORM_XPEHH_FILE
```

##### Output files

Again, selscan generates two output files:

- An **.out** file organized in the following columns:
  - locus ID
  - physical position
  - Population A allele frequency 
  - Integrated haplotype homozygosity of the population A allele (iHHA)
  - Population B allele frequency 
  - Integrated haplotype homozygosity of the population B allele (iHH0)
  - unstandardized XP-EHH
- A **.log** file with the same parameters as in the iHS .log file.

When normalizing, we obtain two extra files:

- **.norm**, which is organized as the **.out** file but with two extra columns:
  - standardized XP-EHH
  - 1 if |standardized XP-EHH| >= 2, and 0 otherwise

In case we have defined an *.err and an *.out file for the process when sending it to the server, the **.err** file will also contain information about the normalization. For XP-EHH, this means a three column table with a single row containing the number of SNPs being analized, their XP-EHH mean and their XP-EHH variance.

### Bibliography

1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. *Nature*, *526*(7571), 68.

Lamason, R. L., Mohideen, M. A. P., Mest, J. R., Wong, A. C., Norton, H. L., Aros, M. C., ... & Sinha, S. (2005). SLC24A5, a putative cation exchanger, affects pigmentation in zebrafish and humans. *Science*, *310*(5755), 1782-1786.

Sabeti, P. C., Varilly, P., Fry, B., Lohmueller, J., Hostetter, E., Cotsapas, C., ... & Schaffner, S. F. (2007). Genome-wide detection and characterization of positive selection in human populations. *Nature*, *449*(7164), 913.

Szpiech, Z. A., & Hernandez, R. D. (2014). selscan: an efficient multithreaded program to perform EHH-based scans for positive selection. Molecular biology and evolution, 31(10), 2824-2827.

Voight, B. F., Kudaravalli, S., Wen, X., & Pritchard, J. K. (2006). A map of recent positive selection in the human genome. *PLoS biology*, *4*(3), e72.

