## Cale Seymour
## July 2023

## This script makes taxon barplots from the amplicon data for the hide
## processing project.

library('data.table')
library('magrittr')
library('ggplot2')
library('ggh4x')
library('magrittr')
library('cowplot')
library('scales')
library('phyloseq')

if (!exists('decon_physeq')) source('decontaminate.R');

dir.create('hides-processing/barplots')

colllapsePercent = 4

bptheme = theme(
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.background = element_blank(),
            plot.background = element_rect(fill='white'),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(size=8, angle=90, vjust=0.5),
            axis.text.y = element_text(size=8),
            axis.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=8, face='bold',angle=90),
            legend.text = element_text(size=8, margin = margin(t = 0, b = 0, unit = "pt")),
            legend.title = element_text(size=8, face='bold'),
            legend.spacing.x = unit(0, 'cm'),
            legend.spacing.y = unit(0, 'cm'),
            panel.border = element_rect(color='black'),
            panel.spacing = unit(0, 'cm'),
            strip.background = element_rect(color='white',fill='#A0A0A0'),
            strip.background.x = element_rect(color='white',fill='#A0A0A0'),
            strip.background.y = element_rect(color='white',fill='#A0A0A0'),
            strip.placement = 'outside',
            legend.key.size = unit(0.05,"cm"),
            legend.key.height = unit(0.05, 'cm'), #change legend key height
            legend.key.width = unit(0.35, 'cm'),
            legend.position = c(0,1),
            legend.justification = c(0,1),
            legend.background=element_blank(),
            legend.box = "vertical",
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
            strip.text = element_text(color = 'white', angle=0, hjust=0.01,face="bold",size=11),
            strip.text.y = element_text(color = 'white', angle=0, hjust=0.01, face="bold", size=8),
            strip.text.x = element_text(color = 'white', angle=0, hjust=0.01, face="bold", size=8))
            
collapseRanks = rank_names(decon_physeq)[c(-1,-7)]
collapsePercentText = paste0('<', colllapsePercent, '%')

tdt = decon_physeq@tax_table %>% data.frame(check.names=FALSE) %>% as.data.table(., keep.rownames=TRUE)
if (attributes(decon_physeq@otu_table)$taxa_are_rows) otus = decon_physeq@otu_table %>% data.frame(check.names=FALSE) %>% as.data.table(., keep.rownames=TRUE)
if (!attributes(decon_physeq@otu_table)$taxa_are_rows) otus = decon_physeq@otu_table %>% t() %>% data.frame(check.names=FALSE) %>% as.data.table(., keep.rownames=TRUE)
samdt = decon_physeq@sam_data %>% data.frame(check.names=FALSE) %>% as.data.table(., keep.rownames=TRUE)
colnames(otus)[1] = 'SV'
colnames(tdt)[1] = 'SV'

## Create a list for barplot objects.
barplots = list()

for(collapseRank in collapseRanks)
{
    ## Read in sample-counts.
    lcdt = otus %>% melt(id.var='SV')
    colnames(lcdt) = c('SV', 'sample', 'count')
    lcdt = lcdt[count > 0,]
    
    ## Add the rank of interest to the counts matrix.
    setkey(lcdt, SV)
    setkey(tdt, SV)
    lcdt[tdt, Rank := get(collapseRank)]

    ## Agglomerate at the rank of interest.
    lcgdt = lcdt[, .(gcount = sum(count)), by=c('sample', 'Rank')]
    lcgdt[,proportion := gcount/sum(gcount), by=sample]
    
    ## Filter.
    flcgdt = rbind(
    lcgdt[proportion >= colllapsePercent/100,],
    lcgdt[proportion < colllapsePercent/100, .(Rank = collapsePercentText, gcount = sum(gcount), proportion = sum(proportion)), by=sample]
    )

    ## Merge all data.
    bpdt = merge(x=flcgdt, y=samdt, by.x="sample", by.y="rn")

    ## Read palette in from file.
    pal = scan('include/palette.txt', what='character',sep='\n')

    ## Aggregate to get the priority of the phylum within the dataset.
    bpagg = aggregate(bpdt$gcount, by=list(Rank=bpdt$Rank), FUN=sum)
    bpaggsort = bpagg[rev(order(bpagg$Rank != collapsePercentText, !grepl('Unclassified', bpagg$Rank), bpagg$x)),]

    ## Assign colors from the palette.
    bpaggsort$color = pal[1:nrow(bpaggsort)]

    ## Store these to make our lives easier
    bpFills = setNames(bpaggsort$color, bpaggsort$Rank)
    bpFills[names(bpFills) == collapsePercentText] = '#C0C0C0'
    bpTaxa = bpaggsort$Rank
    bpLabels = gsub('^D_[0-9]__','',bpaggsort$Rank)
    bpLabels = gsub('\\(.*','',bpLabels)

    bpdt[,Rank := factor(Rank, levels = rev(bpTaxa))]
    ## Build a stacked barplot object.
    bp = ggplot(data=bpdt) +
        geom_bar(aes(x=Replicate, y=proportion, fill=Rank), stat="identity") +
        theme(legend.position = 'bottom',
            axis.text.x = element_text(size=8, angle=90, hjust = 1),
            axis.text.y = element_text(size=8),
            legend.text = element_text(size=8),
            legend.title = element_text(size=8)
        ) +
        ylab('Cumulative Abundance') +
        guides(fill=guide_legend(ncol=ceiling(length(bpLabels)/100))) +
        scale_fill_manual(name=collapseRank, breaks=bpTaxa, values=bpFills, labels=bpLabels) +
        scale_y_continuous(labels = percent_format(accuracy = 1)) +
        facet_nested(. ~ MaterialFac + Day + TreatmentFac, scales='free_x', space='free_x') +
        
        theme_bw() + bptheme
    
    pdf(file = NULL)
    leg = get_legend(bp)
    dev.off()
    
    ## Create the plot
    pp = plot_grid(bp + theme(legend.position = 'none'), leg, rel_widths = c(4,2.5))
    
    ## add the plot to the list of barplots.
    barplots = c(barplots, list(pp))
    
    ## Plot the barplot object.
    outfname = paste0('hides-processing/barplots/taxon-barplot-', collapseRank, '.pdf')
    #png('taxon-barplot.png', width=640*(300/72), height=960*(300/72), res=300)
    pdf(outfname, height = 6, width = 6.5) 
        plot(pp)
    dev.off()
}

names(barplots) = collapseRanks