## Cale Seymour
## June 2023

## This one is short. Read in the metadata.

library('data.table')
library('magrittr')
library('phyloseq')

## Read in metadata file.
metadata = read.table('experimental-design.tsv', sep='\t',
    stringsAsFactors=FALSE, quote='', comment.char='',
    header=TRUE) %>% as.data.table()
    
metadata[,MaterialFac := factor(Material, levels = c('water', 'hide'))]
metadata[,TreatmentFac := factor(Treatment, levels = c('untreated', 'lime'))]
levels(metadata$MaterialFac) = c('Water', 'Hide')
levels(metadata$TreatmentFac) = c('Untreated', 'Lime')
    
colnames(metadata) = colnames(metadata) %>% make.names()
    
setkeyv(metadata, 'X.SampleID')

## Convert to dataframe.
samdata = metadata %>% as.data.frame()
rownames(samdata) = metadata[,X.SampleID]

## Create phyloseq object.
physeq_sample_data = sample_data(samdata)