#This script convert VCF files to MSA files that can be use as BPP inout together with the IMAP.
#lets install SNPRelate
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SNPRelate")
library(SNPRelate)
library(adegenet)
library(vcfR)
library(dplyr)
library(dartR)
#see https://rdrr.io/cran/dartR/man/gl2bpp.html for a more detailed explanation
gl.install.vanilla.dartR()


#read in vcf as vcfR
setwd("~/")
vcfR <- read.vcfR("populations.snps.vcf")

#convert vcfR into a 'genind' object
data<-vcfR2genind(vcfR)

#convert to genlight
gen.gl<-vcfR2genlight(vcfR)

# filtering missing data
gen.gl <- gl.filter.callrate(gen.gl,threshold = 1)
gen.gl <- gl.filter.monomorphs(gen.gl)
#subsampling at random "n" number of loci from genlight data
gen.gl <- gl.subsample.loci(gen.gl,n=20)
gl2bpp(x = test,
       method = 2,
       outfile = "input_for_bpp.txt",
       imap = "Imap.txt")

