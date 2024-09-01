## Cale Seymour
## June 2023

## Calculate alpha diversity metrics from hides processing dataset.

library('data.table')
library('magrittr')
library('vegan')
library('phyloseq')
library('ggplot2')
library('picante')
library('ggpubr')

library('cowplot')

if (!exists('decon_physeq')) source('decontaminate.R')

dir.create('hides-processing/alpha-diversity')

physeq.faithpd = function(physeq, split = TRUE)
## Quick faith's PD calculation using picante and phyloseq.
## Pass it a phyloseq object with a tree.
{
    if (!split) {
        OTU = taxa_sums(physeq)
    }
    else if (split) {
        OTU = as(otu_table(physeq), "matrix")
        if (taxa_are_rows(physeq)) {
            OTU = t(OTU)
        }
    }
    
    picante::pd(samp = OTU, tree = phy_tree(physeq), include.root=FALSE)
}

estimate_richness_faithpd = function(physeq, split = TRUE)
## Access faith's pd function to return a dataframe with all diversity metrics.
{
    richnesses = estimate_richness(physeq, split)
     
    faithpd = matrix(physeq.faithpd(physeq, split)[,1], ncol=1, dimnames=list(sample_names(physeq), 'Faiths.PD'))
    
    out = cbind(richnesses, faithpd)
    return(out)
}

rarefied_physeq = rarefy_even_depth(decon_physeq, sample.size = sample_sums(decon_physeq) %>% min(), rngseed = 1002924, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
diversity = estimate_richness_faithpd(rarefied_physeq)
rn = sample_names(rarefied_physeq) %>% setNames(sample_names(rarefied_physeq) %>% make.names())
rownames(diversity) = rn[rownames(diversity)] ## Fix the names :(
fwrite(diversity %>% as.data.table(keep.rownames=TRUE), 'hides-processing/alpha-diversity/alpha-diversity.tsv.xls', sep='\t')

boxplotdt = merge(rarefied_physeq@sam_data, diversity, by = 'row.names') %>% as.data.table(., keep.rownames=TRUE)
measures = c('Observed', 'Shannon', 'Faiths.PD', 'Simpson')

plots = lapply(measures, function(m)
{
    ggboxplot(boxplotdt, x = 'Day', y = m, fill = 'Treatment') +
          facet_grid(.~ Material, space='free', scales = 'free') +
          ylab(m) +
          xlab('Material + Treatment') +
          ggtitle(m)
})

sink('hides-processing/alpha-diversity/alpha-diversity-anovas.txt')

    for (measure in measures)
    {
        ftext = paste0(measure, ' ~ Material + Treatment')
        cat('ANOVA FORMULA: ', ftext, '\n')
        form = ftext %>% as.formula()
        try({aov(data = boxplotdt, formula = form) %>% summary() %>% print()})
        cat('---\n\n')
    }

sink()

pdf('hides-processing/alpha-diversity/alpha-diversity.pdf', width = 10, height = 7.5)
    plot(plot_grid(plotlist=plots))
dev.off()
