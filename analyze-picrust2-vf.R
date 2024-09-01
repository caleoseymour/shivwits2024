## Cale Seymour
## July 2023

## Parse picrust2 data and run aldex2 to id significant virulence factors
library('data.table')
library('magrittr')
library('ggpubr')
library('ALDEx2')
library('ggplot2')

output_dir = 'hides-processing/q2-picrust2_output/aldex2'

run_aldex2 = function(otus, grouping, aldex.analysis = NULL, output_dir='.', force_reanalyze = FALSE)
{
    dir.create(output_dir)
    g = setNames(make.names(grouping), names(grouping))
    ug = unique(g) %>% sort()
    
    stopifnot(length(ug) == 2)
    
    tt = g %>% unique() %>% paste0(collapse='_vs_')
    
    ## Look to find the file for the analysis, and read it if it exists already
    if (is.null(aldex.analysis) | force_reanalyze)
    {
        outfname_analysis = paste0(output_dir, '/', tt, '-analysis.txt')
        if (file.exists(outfname_analysis) & !force_reanalyze)
        {
            aldex.dt = fread(outfname_analysis, sep ='\t', header=TRUE)
            aldex.analysis = aldex.dt[,-1] %>% as.data.frame()
            rownames(aldex.analysis) = aldex.dt[,rn]
        } else {
            
            aldex.analysis = aldex(otus[,names(g)], g, mc.samples=1000, test="t", effect=TRUE, include.sample.summary=FALSE, denom="all", verbose=TRUE, paired.test=FALSE)
            aldex.dt = aldex.analysis %>% as.data.table(keep.rownames=TRUE)
            fwrite(aldex.dt, outfname_analysis, sep='\t')
        }
    }
    
    aldex.dt = as.data.table(aldex.analysis, keep.rownames=TRUE)
    
    outfname_diagnostic = paste0(output_dir, '/', tt, '-diagnostic.pdf')
    
    pdf(outfname_diagnostic, height = 7.5, width = 10)
        
        aldex.plot(aldex.analysis, type="MA", test="welch", xlab="Log-ratio abundance", ylab="Difference")
        aldex.plot(aldex.analysis, type="MW", test="welch", xlab="Dispersion", ylab="Difference")

        par(mfrow=c(1,2))
        plot(aldex.analysis$effect, aldex.analysis$we.ep, log="y", cex=0.7, col=rgb(0,0,1,0.2),
          pch=19, xlab="Effect size", ylab="P value", main="Effect size plot")
        points(aldex.analysis$effect, aldex.analysis$we.eBH, cex=0.7, col=rgb(1,0,0,0.2),
          pch=19)
        abline(h=0.05, lty=2, col="grey")
        legend(15,1, legend=c("P value", "BH-adjusted"), pch=19, col=c("blue", "red"))

        plot(aldex.analysis$diff.btw, aldex.analysis$we.ep, log="y", cex=0.7, col=rgb(0,0,1,0.2),
          pch=19, xlab="Difference", ylab="P value", main="Volcano plot")
        points(aldex.analysis$diff.btw, aldex.analysis$we.eBH, cex=0.7, col=rgb(1,0,0,0.2),
          pch=19)
        abline(h=0.05, lty=2, col="grey")
    dev.off()
    
    aldex.dt %>% setkey(rn)
    aldex.dt[kos_vf_classification, pathogenesis := classification]
    aldex.dt[pathogenesis_kos[,value], pathogenesis := 'pathogen-associated']
    aldex.dt[pathogenesis_pathways[,ko], pathogenesis := 'toxin or invasion']
    #aldex.dt[kos_vfdb_classification, vfID := vfid]
    
    #found.by.all = which(aldex.analysis$we.eBH < 0.05 & aldex.analysis$wi.eBH < 0.05)
    
    ## use vals found by either the wilcoxon or welch's test.
    found.by.one = which(aldex.analysis$we.eBH < 0.05 | aldex.analysis$wi.eBH < 0.05)
    effect.size = which(abs(aldex.analysis$effect) >= 1)
    
    aldex.dt.signif = aldex.dt[intersect(found.by.one, effect.size),][order(effect),]
    #aldex.dt.signif = aldex.dt[found.by.one,][order(effect),]
    aldex.dt.signif[,otu := factor(rn, levels=rn)]
    
    aldex.dt.signif[,enriched.sign := sign(effect)]
    enriched.groups = lapply(paste0('rab.win.', ug), match, colnames(aldex.analysis)) %>% setNames(ug)
    enriched.group.signs = c(enriched.groups[[1]] - enriched.groups[[2]], enriched.groups[[2]] - enriched.groups[[1]]) %>% sign() %>% setNames(ug)
    
    for(i in 1:length(enriched.group.signs))
    {
        aldex.dt.signif[enriched.sign == enriched.group.signs[i], enriched.group := names(enriched.group.signs)[i]]
    }
    # aldex.dt.signif.pathogen = aldex.dt.signif[vfPA == 'pathogen-associated',]
    # aldex.dt.signif.pathogen = aldex.dt.signif[!is.na(vfID),]

    
  
    
    pdt = aldex.dt.signif[pathogenesis %in% c('pathogen-associated', 'toxin or invasion'),]
    #print(pdt)
    #print(ug)
    p = ggplot(pdt) + 
        geom_col(aes(x = effect, y = otu, fill = enriched.group), color='black', width = 0.8) +
        #geom_text(aes(x = effect, y = otu, label = rn, color = vfPA), x = sign(pdt$effect) * 0.05, hjust = ifelse(pdt$effect > 0, 1, 0), size = 8 * (5/14)) +
        geom_text(aes(x = effect, y = otu, label = rn), x = -1 * sign(pdt$effect) * 0.1, hjust = ifelse(pdt$effect > 0, 1, 0), size = 8 * (5/14)) +
        labs(x = 'Median effect size', y = NULL) +
        theme_pubclean()
   
    if(nrow(pdt) > 0)
    {
        outfname_diffabund = paste0(output_dir, '/', tt, '-diffabund.pdf')
        pdf(outfname_diffabund, height = 10*(0.125 + (nrow(pdt) / 95)), width = 7.5)
            plot(p)
        dev.off()
    }
}

##  Read ref from tables.
kos_vf_classification = fread('ko-virulence-factors.tsv.xls', sep='\t', header=TRUE) %>% setkey(kegg)
pathogenesis_pathways = fread('pathogenesis-pathways.txt', sep='\t', header=TRUE)
pathogenesis_kos = fread('pathogenesis-module-kos.txt', sep='\t', header=TRUE)

## Read in metadata file.
metadata = read.table('experimental-design.tsv', sep='\t',
    stringsAsFactors=FALSE, quote='', comment.char='',
    header=TRUE) %>% as.data.table() %>% setkey(X.SampleID)

## Read in picrust KO outputs.
kos = fread('hides-processing/q2-picrust2_output/ko_metagenome.tsv', sep='\t', header=TRUE)
colnames(kos)[1] = 'KO'
ko_melt =  kos %>% melt(id.var='KO')
ko_melt[,percent := value/sum(value),by=variable]
setkey(ko_melt, KO)

## Split into comparisons & perform ALDEx2
komat = kos[,-1] %>% as.matrix() %>% round()
rownames(komat) = kos[,KO]


groupings = metadata[,setNames(interaction(Material, Treatment, Day) %>% as.character(), X.SampleID)]

## List of comparisons to make.
water.lime.day = groupings[groupings %in% c('water.lime.Day 0', 'water.lime.Day 6')]
try({run_aldex2(otus = komat, grouping = water.lime.day, output_dir = output_dir)})

water.untreated.day = groupings[groupings %in% c('water.untreated.Day 0', 'water.untreated.Day 6')]
try({run_aldex2(otus = komat, grouping = water.untreated.day, output_dir = output_dir)})

water.day6 = groupings[groupings %in% c('water.untreated.Day 6', 'water.lime.Day 6')]
try({run_aldex2(otus = komat, grouping = water.day6, output_dir = output_dir)})

hide.untreated.day = groupings[groupings %in% c('hide.untreated.Day 0', 'hide.untreated.Day 6')]
try({run_aldex2(otus = komat, grouping = hide.untreated.day, output_dir = output_dir)})

hide.lime.day = groupings[groupings %in% c('hide.untreated.Day 0', 'hide.lime.Day 6')]
try({run_aldex2(otus = komat, grouping = hide.lime.day, output_dir = output_dir)})

hide.day6 = groupings[groupings %in% c('hide.untreated.Day 6', 'hide.lime.Day 6')]
try({run_aldex2(otus = komat, grouping = hide.day6, output_dir = output_dir)})

water.untreated.day = groupings[groupings %in% c('water.lime.Day 0', 'water.untreated.Day 6')]
try({run_aldex2(otus = komat, grouping = water.untreated.day, output_dir = output_dir)})