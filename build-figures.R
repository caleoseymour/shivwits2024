## Cale Seymour
## July 2023

## This script acts as an accessor for all of the others when I want to rebuild
## the figures.

library('ggplot2')
library('cowplot')
library('svglite')

dir.create('hides-processing/figures')

## Build figure 1.
if (!exists('barplots')) source('barplots.R')
if (!exists('punifrac')) source('ordinations.R')

figure1 = plot_grid(punifrac, barplots[['Genus']], ncol = 1, rel_heights = c(3.5, 4), labels = c('A', 'B'))
pdf('hides-processing/figures/barplot-unifrac.pdf', height = 7.5, width = 6.5)
    plot(figure1)
dev.off()

svglite('hides-processing/figures/barplot-unifrac.svg', height = 7.5, width = 6.5)
    plot(figure1)
dev.off()


## Build other figures.
source('alpha-diversity.R')
source('lefse.R')