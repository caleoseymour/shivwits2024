## Cale Seymour
## July 2023

## This script makes NMDS ordinations from the amplicon data for the hide
## processing project.

library('data.table')
library('magrittr')
library('vegan')
library('phyloseq')
library('ggplot2')
library('ggrepel')
library('picante')
library('broom')

if (!exists('decon_physeq')) source('decontaminate.R')

percent_decon_physeq = prune_taxa(which(taxa_sums(decon_physeq) > 0) %>% names(), decon_physeq) %>% transform_sample_counts(., function(x) x/sum(x))
samdt = decon_physeq@sam_data %>% data.frame(check.names=FALSE) %>% as.data.table(., keep.rownames=TRUE)
setkey(samdt, rn)

dir.create('hides-processing/ordinations/')

ordination_theme = theme(
        panel.background = element_blank(),
        plot.background = element_rect(fill='white'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=8, hjust=0.5),
        axis.text.y = element_text(size=8, hjust = 1),
        axis.title = element_blank(),
        axis.title.x = element_text(size=8, face='bold'),
        axis.title.y = element_text(size=8, face='bold',angle=90),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8, face='bold'),
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        panel.border = ggplot2::element_rect(color='black'),
        panel.spacing = unit(0, 'cm'),
        strip.background = ggplot2::element_rect(color='black',fill='white'),
        strip.background.x = ggplot2::element_rect(color='black',fill='white'),
        strip.background.y = ggplot2::element_rect(color='black',fill='white'),
        strip.placement = 'outside',
        legend.key.size = unit(0.1,"cm"),
        legend.key.height = unit(0.25, 'cm'), #change legend key height
        legend.key.width = unit(0.25, 'cm'),
        legend.position = c(1,0),
        legend.justification = c(1, 0),
        legend.background=element_blank(),
        legend.box = "horizontal",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
        )

## Run Adonis

## Generate distance matrices.
unifrac = UniFrac(percent_decon_physeq, weighted=TRUE, normalized=TRUE, parallel=FALSE)

if (attributes(decon_physeq@otu_table)$taxa_are_rows) dist.bray = vegdist(percent_decon_physeq@otu_table %>% t(), method='bray', diag=TRUE, upper=TRUE)
if (!attributes(decon_physeq@otu_table)$taxa_are_rows) dist.bray = vegdist(percent_decon_physeq@otu_table, method='bray', diag=TRUE, upper=TRUE)

if (attributes(decon_physeq@otu_table)$taxa_are_rows) nmds.bray = metaMDS(percent_decon_physeq@otu_table %>% t(), distance='bray', trace=0)
if (!attributes(decon_physeq@otu_table)$taxa_are_rows) nmds.bray = metaMDS(percent_decon_physeq@otu_table, distance='bray')

bray.longdist = tidy(dist.bray) %>% as.data.table()
bray.longdist = rbind(bray.longdist, bray.longdist[,.(item1 = item2, item2 = item1, distance)]) %>% unique()

unifrac.longdist = tidy(unifrac) %>% as.data.table()
unifrac.longdist = rbind(unifrac.longdist, unifrac.longdist[,.(item1 = item2, item2 = item1, distance)]) %>% unique()

split_batches = split(samdt, by = c('Day','Material','Treatment')) %>% lapply(function(x) return(x[,rn]))
comparisons = combn(names(split_batches), 2) %>% t() %>% as.data.table()

## print bray curtis adonis to file.
sink('hides-processing/ordinations/pairwise-adonis-braycurtis.txt')
    for (i in 1:nrow(comparisons))
    {
        c1 = comparisons[i,V1]
        cg1 = split_batches[c1] %>% unlist()
        
        c2 = comparisons[i,V2]
        cg2 = split_batches[c2] %>% unlist()
        
        cat('ADONIS: ', cg1, ' vs. ', cg2, '\n')
        
        cbdt = bray.longdist[ item1 %in% cg1 & item2 %in% cg2
                              | item2 %in% cg1 & item1 %in% cg2
                              | item2 %in% cg1 & item1 %in% cg1
                              | item2 %in% cg2 & item1 %in% cg2
                              ,] %>% dcast( item2 ~ item1, value.var = 'distance')
                              
        cbdm = cbdt[,-1] %>% as.matrix()
        rownames(cbdm) = cbdt[,item2]
        cbdist = cbdm %>% as.dist(diag=FALSE, upper=FALSE)
        
        cbsd = data.frame(sample = c(cg1, cg2), group = c(rep(c1, length(cg1)), rep(c2, length(cg2))))
        rownames(cbsd) = cbsd$sample
            
        adonis2(cbdist ~ group, data = cbsd, by="terms") %>% print()
        cat('---\n\n')
    }
sink()

## print unifrac adonis to file.
sink('hides-processing/ordinations/pairwise-adonis-unifrac.txt')
    for (i in 1:nrow(comparisons))
    {
        c1 = comparisons[i,V1]
        cg1 = split_batches[c1] %>% unlist()
        
        c2 = comparisons[i,V2]
        cg2 = split_batches[c2] %>% unlist()
        
        cat('ADONIS: ', cg1, ' vs. ', cg2, '\n')
        
        cbdt = unifrac.longdist[ item1 %in% cg1 & item2 %in% cg2
                               | item2 %in% cg1 & item1 %in% cg2
                               | item2 %in% cg1 & item1 %in% cg1
                               | item2 %in% cg2 & item1 %in% cg2
                               ,] %>% dcast( item2 ~ item1, value.var = 'distance')
                              
        cbdm = cbdt[,-1] %>% as.matrix()
        rownames(cbdm) = cbdt[,item2]
        cbdist = cbdm %>% as.dist(diag=FALSE, upper=TRUE)
        
        cbsd = data.frame(sample = c(cg1, cg2), group = c(rep(c1, length(cg1)), rep(c2, length(cg2))))
        rownames(cbsd) = cbsd$sample
            
        adonis2(cbdist ~ group, data = cbsd) %>% print()
        cat('---\n\n')
    }
sink()



## print bray curtis adonis to file.
sink('hides-processing/ordinations/hides-adonis-braycurtis.txt')
    hide_samples = samdt[Material == 'hide', rn]
        
    cat('ADONIS: HIDES ONLY')
    
    cbdt = bray.longdist[item1 %in% hide_samples & item2 %in% hide_samples,] %>% dcast( item2 ~ item1, value.var = 'distance')
                          
    cbdm = cbdt[,-1] %>% as.matrix()
    rownames(cbdm) = cbdt[,item2]
    cbdist = cbdm %>% as.dist(diag=FALSE, upper=FALSE)
    
    adonis2(cbdist ~ Treatment + Day, data = samdt[rn %in% hide_samples,], by="margin") %>% print()
sink()

sink('hides-processing/ordinations/hides-adonis-unifrac.txt')
    hide_samples = samdt[Material == 'hide', rn]
        
    cat('ADONIS: HIDES ONLY')
    
    print(hide_samples)

    cbdt = unifrac.longdist[item1 %in% hide_samples & item2 %in% hide_samples,] %>% dcast( item2 ~ item1, value.var = 'distance')
                          
    cbdm = cbdt[,-1] %>% as.matrix()
    rownames(cbdm) = cbdt[,item2]
    cbdist = cbdm %>% as.dist(diag=FALSE, upper=FALSE)
    
    adonis2(cbdist ~ Treatment + Day, data = samdt[rn %in% hide_samples,], by="margin") %>% print()

sink()

sink('hides-processing/ordinations/water-adonis-unifrac.txt')
    water_samples = samdt[Material == 'water', rn]
        
    cat('ADONIS: WATER ONLY')
    
    print(water_samples)

    
    cbdt = unifrac.longdist[item1 %in% water_samples & item2 %in% water_samples,] %>% dcast( item2 ~ item1, value.var = 'distance')
                          
    cbdm = cbdt[,-1] %>% as.matrix()
    rownames(cbdm) = cbdt[,item2]
    cbdist = cbdm %>% as.dist(diag=FALSE, upper=FALSE)
    
    adonis2(cbdist ~ Treatment + Day, data = samdt[rn %in% water_samples,], by="margin") %>% print()

sink()

sink('hides-processing/ordinations/water-adonis-braycurtis.txt')
    water_samples = samdt[Material == 'water', rn]
        
    cat('ADONIS: WATER ONLY')
    
    print(water_samples)
    
    cbdt = bray.longdist[item1 %in% water_samples & item2 %in% water_samples,] %>% dcast( item2 ~ item1, value.var = 'distance')
                          
    cbdm = cbdt[,-1] %>% as.matrix()
    rownames(cbdm) = cbdt[,item2]
    cbdist = cbdm %>% as.dist(diag=FALSE, upper=FALSE)
    
    adonis2(cbdist ~ Treatment + Day, data = samdt[rn %in% water_samples,], by="margin") %>% print()

sink()

sink('hides-processing/ordinations/water-day0v6untreated-adonis-braycurtis.txt')
    water_timecourse_samples = c(samdt[Material == 'water' & Day == 'Day 0' & Treatment == 'lime', rn],
        samdt[Material == 'water' & Day == 'Day 6' & Treatment == 'untreated', rn])
        
    cat('ADONIS: DAY 0 VS DAY6 UNTREATED WATER')
    
    print(water_timecourse_samples)
    
    cbdt = bray.longdist[item1 %in% water_timecourse_samples & item2 %in% water_timecourse_samples,] %>% dcast( item2 ~ item1, value.var = 'distance')
                          
    cbdm = cbdt[,-1] %>% as.matrix()
    rownames(cbdm) = cbdt[,item2]
    cbdist = cbdm %>% as.dist(diag=FALSE, upper=FALSE)
    
    adonis2(cbdist ~ Day, data = samdt[rn %in% water_timecourse_samples,]) %>% print()

sink()


sink('hides-processing/ordinations/water-day0v6untreated-adonis-unifrac.txt')
    water_timecourse_samples = c(samdt[Material == 'water' & Day == 'Day 0' & Treatment == 'lime', rn],
        samdt[Material == 'water' & Day == 'Day 6' & Treatment == 'untreated', rn])
        
    cat('ADONIS: DAY 0 VS DAY6 UNTREATED WATER')
    
    print(water_timecourse_samples)
    
    cbdt = unifrac.longdist[item1 %in% water_timecourse_samples & item2 %in% water_timecourse_samples,] %>% dcast( item2 ~ item1, value.var = 'distance')
                          
    cbdm = cbdt[,-1] %>% as.matrix()
    rownames(cbdm) = cbdt[,item2]
    cbdist = cbdm %>% as.dist(diag=FALSE, upper=FALSE)
    
    adonis2(cbdist ~ Day, data = samdt[rn %in% water_timecourse_samples,]) %>% print()

sink()


sink('hides-processing/ordinations/water-day0v6lime-adonis-braycurtis.txt')
    water_timecourse_samples = c(samdt[Material == 'water' & Day == 'Day 0' & Treatment == 'lime', rn],
        samdt[Material == 'water' & Day == 'Day 6' & Treatment == 'lime', rn])
        
    cat('ADONIS: DAY 0 VS DAY6 LIME WATER')
    
    print(water_timecourse_samples)
    
    cbdt = bray.longdist[item1 %in% water_timecourse_samples & item2 %in% water_timecourse_samples,] %>% dcast( item2 ~ item1, value.var = 'distance')
                          
    cbdm = cbdt[,-1] %>% as.matrix()
    rownames(cbdm) = cbdt[,item2]
    cbdist = cbdm %>% as.dist(diag=FALSE, upper=FALSE)
    
    adonis2(cbdist ~ Day, data = samdt[rn %in% water_timecourse_samples,]) %>% print()

sink()


sink('hides-processing/ordinations/water-day0v6lime-adonis-unifrac.txt')
    water_timecourse_samples = c(samdt[Material == 'water' & Day == 'Day 0' & Treatment == 'lime', rn],
        samdt[Material == 'water' & Day == 'Day 6' & Treatment == 'lime', rn])
        
    cat('ADONIS: DAY 0 VS DAY6 LIME WATER')
    
    print(water_timecourse_samples)
    
    cbdt = unifrac.longdist[item1 %in% water_timecourse_samples & item2 %in% water_timecourse_samples,] %>% dcast( item2 ~ item1, value.var = 'distance')
                          
    cbdm = cbdt[,-1] %>% as.matrix()
    rownames(cbdm) = cbdt[,item2]
    cbdist = cbdm %>% as.dist(diag=FALSE, upper=FALSE)
    
    adonis2(cbdist ~ Day, data = samdt[rn %in% water_timecourse_samples,]) %>% print()

sink()

## Run NMDS.
ord.bray = nmds.bray$points %>% as.data.table(keep.rownames=TRUE)
colnames(ord.bray)[1] = 'sample'
stress.bray = nmds.bray$stress

nmds.unifrac = metaMDS(unifrac, trace=0)
ord.unifrac = nmds.unifrac$points %>% as.data.table(keep.rownames=TRUE)
colnames(ord.unifrac)[1] = 'sample'
stress.unifrac = nmds.unifrac$stress
ord.bray %>% setkey(sample)

## Print NMDS to file.
pd.bray = samdt[ord.bray,]
pd.unifrac = samdt[ord.unifrac,]

pd.unifrac.mean = pd.unifrac[,.(MDS1 = mean(MDS1), MDS2 = mean(MDS2)), by = .(Label, Treatment)]

punifrac = ggplot(pd.unifrac) +
    geom_polygon(aes(x = MDS1, y = MDS2, fill = Treatment, color=Treatment, group = interaction(Treatment, Material, Day)), alpha=0.3, linewidth = 0.75) +
    geom_point(aes(x = MDS1, y = MDS2, color = Treatment, shape = Material), size = 2) +
    geom_point(data=pd.unifrac.mean, aes(x = MDS1, y = MDS2), size = 2, color='black') +
    geom_text_repel(data=pd.unifrac.mean, aes(x = MDS1, y = MDS2, label = Label), size = 8 * (5/14)) +
    #geom_text_repel(aes(x = MDS1, y = MDS2, label = Label), size = 8 * (5/14)) +
    geom_text(x = -0.08, y = -0.10, label = paste0('Stress=', round(stress.unifrac, 3)), hjust=0, size = 8 * (5/14)) +

    theme_bw() + ordination_theme
            
pdf('hides-processing/ordinations/ordination-NMDS-unifrac.pdf', height = 8, width=8) 
    plot(punifrac)
dev.off()

pbray = ggplot(pd.bray) +
    geom_polygon(aes(x = MDS1, y = MDS2, fill = Treatment, color=Treatment, group = interaction(Treatment, Material, Day)), alpha=0.3, linewidth = 0.75) +
    geom_point(aes(x = MDS1, y = MDS2, color = Treatment, shape = Material), size = 2) +
    geom_text_repel(aes(x = MDS1, y = MDS2, label = Label), size = 8 * (5/14)) +
    geom_text(x = 1.3, y = -0.66, label = paste0('Stress = ', round(stress.bray, 3)), hjust=0, size = 8 * (5/15)) +

    theme_bw() + ordination_theme
            
pdf('hides-processing/ordinations/ordination-NMDS-bray.pdf', height = 4, width=6.5) 
    plot(pbray)
dev.off()