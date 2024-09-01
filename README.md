# shivwits2024
This repository contains all of the code used to process microbiome data for the 2024 Shivwits lime hide-tanning experiments. Precise methods are available in the final publication.

## Description of files
Filename | Description
------------ | -------------
pipeline.sh | Run the whole analysis (you'll need to match the input filenames).
build-figures.R | Make the figures for the paper.
analyze-picrust2-vf.R | Analyze virulence-associated KOs using picrust2 data.
decontaminate.R | Run contaminant-ASV removal scripts (note that this repos DOES NOT include Sourcetracker 1.0.1, which should be placed in /include/sourcetracker-1.0.1)
read-asvs.R | Read asvs from file.
read-metadata.R | Read metadata from file.
read-taxonomy.R | Read taxonomy from file.
barplots.R | Generate barplots.
ordinations.R | Generate ordinations.
alpha-diversity.R | Run alpha div and plot to output folder.
lefse.R | Run lefse.
