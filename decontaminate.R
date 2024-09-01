## Cale Seymour
## June 2023

## Assess contamination and run the sourcetracker-based decontamination protocol
library('decontam')
library('data.table')
library('magrittr')

source('include/runSourceTracker.R')

if (!exists('physeq_otu_table')) source('read-asvs.R');
if (!exists('physeq_sample_data')) source('read-metadata.R');
if (!exists('physeq_tax_table')) source('read-taxonomy.R');

## build inpust for source tracker.
outdir = 'hides-processing/decontaminated-asvs'
filebase = 'decon.asvs'
dir.create(outdir)

experimental_vars = c('Material', 'Treatment', 'Day')

## Remove mitochondria and chloroplasts.
physeq = phyloseq(physeq_otu_table, physeq_sample_data, physeq_tax_table, physeq_phy_tree)
mitochondria_chloroplast = subset_taxa(physeq, Family == 'Mitochondria' | Order == 'Chloroplast') %>% taxa_names()

prune_physeq = prune_taxa(setdiff(taxa_names(physeq), mitochondria_chloroplast), physeq)

contamdf_freq = isContaminant(physeq, method="frequency", conc="DNA.Conc.ng.ul")
contaminant_freq_asvs = contamdf_freq[contamdf_freq$contaminant == TRUE,] %>% rownames()
percent_contamination = sample_sums(prune_taxa(contaminant_freq_asvs, prune_physeq)) / sample_sums(prune_physeq) %>% as.data.frame()
colnames(percent_contamination) = 'Percent.Contamination.Decontam'
fwrite(percent_contamination %>% as.data.table(keep.rownames=TRUE), paste0(outdir, '/percent-contamination-decontam.txt'), sep = '\t')

## SourceTracker
st_ps = prune_physeq
st_otus = as.data.frame(st_ps@otu_table) %>% t()
rarefaction = sample_sums(st_ps) %>% min()
taxa = as.data.frame(st_ps@tax_table)
metadata = as.data.frame(st_ps@sam_data)
metadata$Env = apply(metadata[,experimental_vars], 1, paste, collapse = '.') %>% make.names()

## Run sourcetracker.
if (!file.exists('hides-processing/decontaminated-asvs/full_results/decon.asvs_Unknown_contributions-absolute.txt')) st_results = runSourceTracker(st_otus, metadata, outdir, filebase, rarefaction);

## Read in the decontaminated ASV file.
decontaminated_asvs = fread('hides-processing/decontaminated-asvs/full_results/decon.asvs_Unknown_contributions-absolute.txt', sep='\t', header=TRUE)
decontaminated_asvs = decontaminated_asvs[,colnames(decontaminated_asvs) %>% setdiff(., mitochondria_chloroplast), with=FALSE] ## Drop mitochondria_chloroplast columns, if they still exist.

decontaminated_asv_matrix = decontaminated_asvs[,-1] %>% as.matrix()
rownames(decontaminated_asv_matrix) = decontaminated_asvs[,SampleID]

decon_physeq_otu_table = otu_table(decontaminated_asv_matrix, taxa_are_rows = FALSE)

## Return decontaminated phyloseq object.
decon_physeq = phyloseq(decon_physeq_otu_table, physeq_sample_data, physeq_tax_table, physeq_phy_tree)

## Run decontam to evaluate the difference in contamintion after filtering.
decontam_contaminants_decon = isContaminant(decon_physeq, method="frequency", conc="DNA.Conc.ng.ul")
decontam_contaminants_decon_asvs =  decontam_contaminants_decon[decontam_contaminants_decon$contaminant == TRUE,] %>% rownames()
contamdf_freq = isContaminant(decon_physeq, method="frequency", conc="DNA.Conc.ng.ul")
contaminant_freq_asvs = contamdf_freq[contamdf_freq$contaminant == TRUE,] %>% rownames()
percent_contamination = sample_sums(prune_taxa(contaminant_freq_asvs, decon_physeq)) / sample_sums(decon_physeq) %>% as.data.frame()
colnames(percent_contamination) = 'Percent.Contamination.Decontam'
fwrite(percent_contamination %>% as.data.table(keep.rownames=TRUE), paste0(outdir, '/percent-contamination-postfilter.txt'), sep = '\t')

## Create a standard format biom table, for reimport into qiime2 (for picrust)
std_format_fname = 'hides-processing/decontaminated-asvs/full_results/biom-format-tsv.txt'
std_format_tsv = decon_physeq_otu_table %>% t() %>% as.data.table(keep.rownames=TRUE)
colnames(std_format_tsv)[1] = '#OTU ID'

## Write formatted table to file.
cat('#Formatted for biom\n', file = std_format_fname)
fwrite(std_format_tsv, std_format_fname, sep='\t', append=TRUE, col.names=TRUE)
