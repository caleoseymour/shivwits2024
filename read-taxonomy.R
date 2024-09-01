## Cale Seymour
## June 2023

## Read in the taxonomy & make it look pretty.

library('data.table')
library('magrittr')

## Taxonomy

## Read in taxon table.
taxondt = fread('hides-processing/taxonomy-250f-250r-dada2.tsv', sep='\t', header=TRUE)

## Split taxonomy on ';' to separate taxa, and then identify the rank using the '__' prefix tag.
taxondt = taxondt[,.(taxon = strsplit(Taxon, '; ') %>% do.call('c',.)), by=.(SV=get('Feature ID'))]
taxondt[,rank := gsub('__.*','', taxon)]

## Grab all unique ranks and SVs to generate the rectangular matrix.
ranks = taxondt[,rank] %>% unique()
svs = taxondt[,SV] %>% unique()
form = expand.grid(factor(ranks, levels=ranks), svs) %>% as.data.table() %>% setkey(Var2, Var1)

## Apply order to the ranks.
taxondt[,rank := factor(rank, levels=ranks)]
taxondt %>% setkey(SV, rank)
taxondt = taxondt[form,]

## Where the rank does not exist, make it blank.
taxondt[is.na(taxon), taxon := '']

## Assign if unclassified.
taxondt[,unct := taxon %>% gsub('.__','Unclassified ',.), by=SV]
taxondt[,cl := !(taxon == ''
         | grepl('.__$', taxon)
         | grepl('uncultured', taxon)
         | grepl('unidentified', taxon)
         | grepl('Incertae_Sedis', taxon)
         | grepl('metagenome', taxon))
         , by=SV]

## Cumulative sum to identify the first unclassified rank per ASV.
taxondt[,uncg := cumsum(cl) %>% as.character()]

## Remove duplicated cumsum levels and key the table by unique ASVs
uncguide = taxondt[!duplicated(uncg), .(unct, uncg)] %>% setkey(uncg)
taxondt %>% setkey(uncg)

## Where unclassified, enforce that the displayed taxon is that of the first unclassified rank.
taxondt[cl == FALSE, unct := uncguide[uncg, unct]]
taxondt[cl == FALSE, taxon := paste0(rank, '__', unct)]

## Remove prefix tags, since we no longer need them at this point.
taxondt[,taxon := gsub('.__','',taxon)]

## Apply standardized ranks.
taxonomy = dcast(taxondt, SV ~ rank, value.var = 'taxon')[,.(SV, Kingdom = d, Phylum = p, Class = c, Order = o, Family = f, Genus = g, Species = s)] %>% setkey(SV)
taxondf = taxonomy[,-1] %>% as.data.frame() %>% as.matrix()

## If taxondf is NA, set to 'unclassified'. (This shouldn't happen?)
taxondf[is.na(taxondf)] = 'unclassified'
taxondf = gsub('.__','',taxondf)
class(taxondf) = union(class(taxondf), 'character')
rownames(taxondf) = taxonomy[,SV]

## Make a phyloseq object.
physeq_tax_table = tax_table(taxondf)
