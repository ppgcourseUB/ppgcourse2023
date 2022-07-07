suppressPackageStartupMessages(library(mjolnir))

# Define input fastq files (only names of R1 files are needed)
R1_filenames <-c("MARKET_BARCODES_R1.fastq.gz")

# Define number of cores to be used in parallel.
cores <- 4

# Input names of the individual libraries to be used. It should be a 4-character name, matching the information of the ngsfilter files
lib_prefixes <- c("MRKT")

# Input name for the final combined library (should be a 4-character name)
lib <- "MRKT"

####################
# MJOLNIR pipeline #
####################

# RAN will distribute R1 & R2 fastq files into equal-sized pieces for parallel computing
mjolnir1_RAN(R1_filenames,cores,lib_prefixes,R1_motif="_R1.",R2_motif="_R2.") 

# FREYJA will do the paired-end alignment, demultiplexing & length filtering
mjolnir2_FREYJA(lib_prefixes,cores,Lmin=140, Lmax=190,lib)

# HELA will remove chimaeric sequences in a sample-by-sample basis, will change identifiers of remaining unique sequences & will generate a table of their abundances in each sample & a fasta file with unique sequences and their total abundance for ODIN
mjolnir3_HELA(lib,cores,remove_singletons = FALSE)

# ODIN will do the clustering & will generate a table with the abundances of each MOTU in each sample
mjolnir4_ODIN(lib,cores,d=1,generate_ESV=FALSE)

# THOR will asign the taxonomy to the representative sequence of each MOTU
mjolnir5_THOR(lib, cores, tax_dir="/home/ubuntu/taxo", ref_db="DUFA_MiFish_20220106.fasta ", taxo_db="taxo_NCBI_20210720", run_ecotag=T)

# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment from ODIN & THOR in a single table
mjolnir6_FRIGGA(lib)

# LOKI kill remove the pseudogenes and will keep track of the taxonomic information of the removed MOTUs
mjolnir7_LOKI(lib)

# RAGNAROC will change the names of the samples to recover the original names and will remove unnecessary columns
mjolnir8_RAGNAROC(lib, "MRKT_metadata.tsv", "MRKT_final_dataset.tsv", sort_MOTUs="taxonomy", remove_bacteria=T, remove_contamination=F, min_reads=2)

