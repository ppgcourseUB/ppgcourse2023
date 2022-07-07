# Metabarcoding: the MJOLNIR pipeline

In this session you will learn how to run the MJOLNIR pipeline to perform a metabarcoding analysis of a multiplexed amplicon library.

##Data required:

Available from: https://drive.google.com/drive/folders/1zkX01xJov68FFObgbaNYehWD_3dAbDe3?usp=sharing 

## Background and description of the data set:

The files MARKET_BARCODES_R1.fastq.gz and MARKET_BARCODES_R2.fastq.gz contain paired-end reads from an Illumina MiSeq run of a multiplexed library built from 80 samples of fish individuals taken from a fish market in Singapore. The extracted DNA from these individuals was amplified by PCR using the MiFish 12S metabarcoding primers: 
- Forward primer: GTCGGTAAAACTCGTGCCAGC,
- Reverse primer: CATAGTGGGGTATCTAATCCCAGTTTG

Individual  primers  were  tagged  in the 5' end with 7-bp twin-tags (both tags equal in each pair of forward and reverse primers). Sequences of the tags used for each fish sample are available in the file Dataset market fish.csv. 

We will run the MJOLNIR pipeline to retrieve the most abundant sequence in each sample with the goal of identifying the fish species.


##Tasks:

a) Write the ngsfilter_MRKT.tsv file needed to demultiplex these 80 samples, in a format that is compatible with MJOLNIR. 

b) Write the MRKT_metadata.tsv file needed by MJOLNIR_RAGNAROC, in a format that is compatible with MJOLNIR.

c) Create a new folder called "MRKT" in the server. Upload the fastq files and the two needed table files, and run the MJOLNIR pipeline for these samples using 4 computing cores. Use the adequate values for 12S MiFish for running FREYJA and ODIN (you can get this information from https://github.com/uit-metabarcoding/MJOLNIR ), and use the default values for all other steps. Use the reference and taxonomy databases for MiFish data available from: https://drive.google.com/drive/folders/1Y_NWhNctHs5VqBcwVeMl5QTN5-2EBuF3?usp=sharing  

d) How many MOTUs are there in the final dataset? How many total reads? 

e) Open the final table in R. Calculate the total reads for each sample and remove those samples with less than 50 reads (we will consider them as failed PCR amplifications). How many samples remain with more than 50 reads (successful identifications)? 

f) For the first ten fishes in the data set (Market_fish_1 to Market_fish_10), look at the taxonomic assignment for the most abundant MOTU in each sample.   Then check the best identity provided by MJOLNIR and if it is lower than 100%, then use the Genbank to look for a better match for the sequence.
