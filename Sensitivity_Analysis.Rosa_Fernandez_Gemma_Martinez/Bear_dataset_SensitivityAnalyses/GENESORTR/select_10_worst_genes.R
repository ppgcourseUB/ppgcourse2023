library(phangorn)
library(ape)
library(tibble)
library(dplyr)
library(tidyr)
library(phytools)
library(adephylo)

setwd('.')

#INPUT: fill this with file names
alignment <- 'sorted_alignment_50genes.fa'
partition <- 'sorted_alignment_50genes.txt'
species_tree <- 'bear_species_tree.tre'
gene_trees <- 'sorted_trees_50genes.tre'

#INPUT: is the alignment 'DNA' or 'AA'
type <- 'AA'

##INPUT: Desired number of genes to retain
#if n_genes == 'all' then the dataset is sorted but not subsampled.
n_genes <- '10'

data <- read.phyDat(alignment, format = 'fasta', type = type)
if(type == 'AA') {
  data <- as.AAbin(data) 
} else {
  data <- as.DNAbin(data)
}

partitions <- read.table(partition, sep = ' ')
names <- as.character(unlist(enframe(partitions[,(which(partitions[1,] == '=') - 1)], name = NULL), use.names = F))
partitions <- enframe(partitions[,ncol(partitions)], name = NULL)
partitions <- partitions %>% separate(value, into = c('Start', 'End'), sep = '-') %>% 
  mutate_if(is.character, as.numeric)

gene_trees <- read.tree(gene_trees)
species_tree <- read.tree(species_tree)

genes <- 1:length(gene_trees)

if(n_genes == 'all') {
  n_genes <- nrow(variables)
  cut <- F
} else {
  if(is.character(n_genes)){
    n_genes <- as.numeric(n_genes)
  }
  cut <- T
}

#subsample (or not if n_genes was left == 0)
#Change value of n_genes
sorted_names <- names[length(genes)-n_genes:length(genes)]
sorted_partitions <- partitions[length(genes)-n_genes:length(genes),]
sorted_data <- data[,1:sorted_partitions$End[n_genes]]
sorted_trees <- gene_trees[length(genes)-n_genes:length(genes)]

write.phyDat(phyDat(sorted_data, type = type), 
             file = paste0(getwd(), '/sorted_alignment_', 
                           n_genes, '_worst_genes.fa'), format = 'fasta')
partitions_tosave <- paste0(sorted_names, ' = ', sorted_partitions$Start, '-', sorted_partitions$End)
write(partitions_tosave, file = paste0(getwd(), '/sorted_alignment_', n_genes, '_worst_genes.txt'))
write.tree(sorted_trees, file = paste0(getwd(), '/sorted_trees_', n_genes, '_worst_genes.tre'))
