## Cale Seymour
## June 2023

## Run lefse analyses on hide processing data.

library('data.table')
library('magrittr')
library('microbiomeMarker')
assignInNamespace("fix_na_tax", value = function(x) x, ns = "microbiomeMarker")
library('phyloseq')
library('ggplot2')
library('dplyr')

dir.create('hides-processing/lefse/')

make_lefse_plot = function(physeq, group, enrich_group = NULL, levels = NULL)
{
    lef = run_lefse(physeq, group = group, multigrp_strat = TRUE)

    marker = lef %>% microbiomeMarker:::marker_table()
    marker$logp = -log10(marker$padj)

    if (is.null(enrich_group)) enrich_group = (sample_data(physeq)[,group] %>% unlist() %>% unique())[1];
    
    wti = marker$enrich_group == enrich_group
    marker$ef_lda[wti] = 0 - marker$ef_lda[wti]

    marker = data.frame(marker) %>% arrange(., .data$ef_lda)
    feat = marker$feature
    marker$feature = factor(feat, levels = feat)

    marker$feature_label = marker$feature %>% as.character() %>% sapply(microbiomeMarker:::get_features_labels, 1, 60)
    marker$level = marker$feature_label %>% gsub('__.*','',.)

    return(marker)
}


plot_lefse_plot = function(marker, group, enrich_group = NULL, levels = NULL)
{
    ## Remove levels as appropriate.
    if (!is.null(levels))
    {
        marker = marker[marker$level %in% levels,] %>% as.data.table()
    }
    marker = marker[!grepl('unclassified', marker$feature_label),]

    label_x = "LDA score (log10)"
    p = ggplot(marker, aes(ef_lda, feature)) + 
        geom_col(aes(fill = enrich_group), color='black', width = 0.8) +
        geom_text(aes(label = feature_label), x = 0 - sign(marker$ef_lda) * 0.05, hjust = ifelse(marker$ef_lda > 0, 1, 0), size = 8 * (5/14)) +
        labs(x = label_x, 
             y = NULL, fill = "Enriched group") +
        #scale_fill_manual(values = c('#008000', '#FF0000'), breaks = c('water', 'lime')) +
        #scale_fill_manual(breaks = c('ko', 'wildtype')) +
        scale_x_continuous(expand = c(0, 0), breaks=c(-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6), limits = c(-6, 6)) +
        theme_bw() +
        theme(
            legend.position=c(0.01,0.99),
            panel.background = element_blank(),
            plot.background = element_rect(fill='white'),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(size=8, vjust=1, hjust = 0.5),
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            axis.title.x = element_text(size=8, face='bold'),
            axis.title.y = element_text(size=8, face='bold'),
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
            legend.key.size = unit(0.25,"cm"),
            legend.key.height = unit(0.25, 'cm'), #change legend key height
            legend.key.width = unit(0.25, 'cm'),
            legend.justification = c(0,1),
            #legend.background=element_blank(),
            legend.box = "vertical",
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
        )
}

physeq_Day6_water = prune_samples(sample_names(decon_physeq)[decon_physeq@sam_data$Day == 'Day 6' & decon_physeq@sam_data$Material == 'water'], decon_physeq)
pdf('hides-processing/lefse/lefse_water-Day6_water-vs-lime.pdf', height = 10, width = 7.5)
    lefse_Day6_water = make_lefse_plot(physeq_Day6_water, group = 'Treatment', enrich_group = 'lime')
    plot_lefse_plot(lefse_Day6_water, group = 'Treatment', enrich_group = 'lime', levels = 'g') %>% plot()
dev.off()

fwrite(lefse_Day6_water, 'hides-processing/lefse/lefse_water-Day6_water-vs-lime.tsv.xls', sep='\t')


# physeq_Day0_water = prune_samples(sample_names(decon_physeq)[decon_physeq@sam_data$Day == 'Day 0' & decon_physeq@sam_data$Material == 'water'], decon_physeq)
# pdf('lefse_water-Day0_water-vs-lime.pdf', height = 10, width = 7.5)
    # lefse_Day0_water = make_lefse_plot(physeq_Day0_water, group = 'Treatment', enrich_group = 'lime')
    # plot_lefse_plot(lefse_Day0_water, group = 'Treatment', enrich_group = 'lime', levels = 'g')
# dev.off()

# physeq_Day0_hides = prune_samples(sample_names(decon_physeq)[decon_physeq@sam_data$Day == 'Day 0' & decon_physeq@sam_data$Material == 'water'], decon_physeq)
# pdf('lefse_hides-Day0_water-vs-lime.pdf', height = 10, width = 7.5)
    # lefse_Day0_hides = make_lefse_plot(physeq_Day0_hides, group = 'Treatment', enrich_group = 'lime'))
# dev.off()

physeq_Day6_hides = prune_samples(sample_names(decon_physeq)[decon_physeq@sam_data$Day == 'Day 6' & decon_physeq@sam_data$Material == 'hide'], decon_physeq)
pdf('hides-processing/lefse/lefse_hides-Day6_water-vs-lime.pdf', height = 10, width = 7.5)
    lefse_Day6_hides = make_lefse_plot(physeq_Day6_hides, group = 'Treatment', enrich_group = 'lime', levels = 'g')
    plot_lefse_plot(lefse_Day6_hides, group = 'Treatment', enrich_group = 'lime', levels = 'g') %>% plot()
dev.off()
fwrite(lefse_Day6_hides, 'hides-processing/lefse/lefse_hides-Day6_water-vs-lime.tsv.xls', sep='\t')

## Lefse within material & treatment on day
# physeq_lime_hides = prune_samples(sample_names(decon_physeq)[decon_physeq@sam_data$Treatment == 'lime' & decon_physeq@sam_data$Material == 'hide'], decon_physeq)
# pdf('lefse_hides-untreated_Day0-vs-Day6.pdf', height = 10, width = 7.5)
    # lefse_lime_hides = make_lefse_plot(physeq_lime_hides, group = 'Day', enrich_group = 'Day 6')
# dev.off()


# physeq_untreated_hides = prune_samples(sample_names(decon_physeq)[decon_physeq@sam_data$Treatment == 'untreated' & decon_physeq@sam_data$Material == 'hide'], decon_physeq)
# pdf('hides-processing/lefse/lefse_hides-untreated_Day0-vs-Day6.pdf', height = 10, width = 7.5)
   # lefse_untreated_hides = make_lefse_plot(physeq_untreated_hides, group = 'Day', enrich_group = 'Day 6', levels = 'g')
# dev.off()

# physeq_untreated_water = prune_samples(sample_names(decon_physeq)[decon_physeq@sam_data$Treatment == 'untreated' & decon_physeq@sam_data$Material == 'water'], decon_physeq)
# pdf('lefse_untreated-water_Day0-vs-Day6.pdf', height = 10, width = 7.5)
    # lefse_untreated_water = make_lefse_plot(physeq_untreated_water, group = 'Day', enrich_group = 'Day 6', levels = 'g')
# dev.off()

# physeq_lime_water = prune_samples(sample_names(decon_physeq)[decon_physeq@sam_data$Treatment == 'lime' & decon_physeq@sam_data$Material == 'water'], decon_physeq)
# pdf('hides-processing/lefse/lefse_lime-water_Day0-vs-Day6.pdf', height = 10, width = 7.5)
    # lefse_lime_water = make_lefse_plot(physeq_lime_water, group = 'Day', enrich_group = 'Day 6', levels = 'g')
# dev.off()