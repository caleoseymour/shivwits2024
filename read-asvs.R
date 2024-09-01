## Cale Seymour
## June 2023

## Read in NOT-YET-DECONTAMINATED OTUs from file.

library('data.table')
library('magrittr')
library('phyloseq')
library('phytools')
library('ape')

otus = fread('hides-processing/table-250f-250r-dada2.filter.tsv', sep='\t', header=TRUE)
colnames(otus)[1] = 'SV'
otumat = otus[,-1] %>% as.matrix()
rownames(otumat) = otus[,SV]

## Read in the tree as well.
physeq_otu_table = phyloseq(otu_table(otumat, taxa_are_rows=TRUE))
tree = read.tree('hides-processing/rep-seqs-250f-250r-dada2.filter.dash.gappyout.nwk') %>% midpoint.root() %>% ladderize()
physeq_phy_tree = phyloseq(phy_tree(tree))