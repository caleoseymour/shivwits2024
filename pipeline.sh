#!/bin/bash

## Cale Seymour
## June 2023
## This pipeline runs qiime2.
. include/q2tils.sh

unzip 2023_06_01_051523BH2illcus515F_Analysis_Pipeline.Zip 

READ_DIR=2023_06_01_051523BH2illcus515F_Raw_Data_UDI/demux
SAMPLE_METADATA=2023_06_01_051523BH2illcus515F_Raw_Data_UDI/sample-metadata.tsv
OUTPUT_DIR='hides-processing'
#export Q2TILS_LOGFILE="bzB4QKgZFTdUq.logfile.txt"
#export Q2TILS_UID="bzB4QKgZFTdUq"


[ ! -d $READ_DIR ] && unzip 2023_06_01_051523BH2illcus515F_Raw_Data_UDI.zip # Raw data from the sequencer.
export OMP_NUM_THREADS=1


## Run fastqc.
mkdir $OUTPUT_DIR/fastqc
for readfile in $(find $READ_DIR -name "*.fastq*")
do
    rf=$(basename $readfile)
    fastqc -o $OUTPUT_DIR/fastqc -f fastq $readfile
done

## Build and check manifest.
q2tils_make_manifest $READ_DIR $SAMPLE_METADATA
stop_on_error $?

q2tils_check_manifest $READ_DIR $READ_DIR/manifest.txt $SAMPLE_METADATA
stop_on_error $?

## Import files.
q2tils_import_reads $READ_DIR $OUTPUT_DIR
stop_on_error $?

## Try to guess the read truncation values so I don't have to go back every time
## and manually reassign it.
## (This assigns the environment vars $Q2TILS_FORWARD_TRUNC and $Q2TILS_REVERSE_TRUNC)
q2tils_analyze_read_truncation $OUTPUT_DIR/paired-end-demux.qzv $OUTPUT_DIR
stop_on_error $?

## Run dada2.
q2tils_dada2 $OUTPUT_DIR/paired-end-demux.qza $OUTPUT_DIR
stop_on_error $?

## Classify.
q2tils_classify_nb $OUTPUT_DIR/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza $OUTPUT_DIR
stop_on_error $?

## Export for the next steps.
q2tils_export $OUTPUT_DIR/taxonomy-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza taxonomy.tsv
mv taxonomy-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.tsv $OUTPUT_DIR/taxonomy-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.tsv

q2tils_export $OUTPUT_DIR/table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza feature-table.biom
mv table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.biom $OUTPUT_DIR/table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.biom

q2tils_export $OUTPUT_DIR/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza dna-sequences.fasta
mv rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.fasta $OUTPUT_DIR/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.fasta


q2tils_qiime2_execute biom convert --to-tsv -i $OUTPUT_DIR/table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.biom -o $OUTPUT_DIR/table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.tsv

## Mothur is better for alignments, so we're using that.
q2tils_mothur_align_procedure $OUTPUT_DIR/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.fasta $OUTPUT_DIR/taxonomy-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.tsv

## Filter table to remove ASVs unassigned at the domain level.
cp $OUTPUT_DIR/table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.tsv $OUTPUT_DIR/table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.filter.tsv
for accno in $(cat $Q2TILS_UID.temp.unassigned.accnos)
do
    sed -i "/^$accno\t.*/d" $OUTPUT_DIR/table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.filter.tsv
done

## Replace dots in the alignment with dashes, and then run the alignment through fasttree (for unifrac the tree doesn't need to be perfect, but it's easy to make a good alignment).
awk '/^[^>]/{gsub("\\.","-");print $0;next;}1' $Q2TILS_UID.temp.pick.filter.fasta > $OUTPUT_DIR/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.filter.dash.fasta
trimal -in $OUTPUT_DIR/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.filter.dash.fasta -out $OUTPUT_DIR/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.filter.dash.gappyout.fasta -gappyout
FastTree -fastest -noml -nome -bionj -nt $OUTPUT_DIR/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.filter.dash.gappyout.fasta > $OUTPUT_DIR/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.filter.dash.gappyout.nwk44

## Stopping point -- build the figures for the paper.
Rscript build-figures.R

## Read in the decontaminated ASV file.
asv_file='hides-processing/decontaminated-asvs/full_results/biom-format-tsv.txt'
biom_file='hides-processing/decontaminated-asvs/full_results/biom-format-biom.biom'
q2tils_qiime2_execute biom convert -i $asv_file -o $biom_file --to-hdf5

## Reimport decontaminated OTU table.
q2tils_qiime2_execute qiime tools import \
  --input-path $biom_file \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path $biom_file.qza

## Run picrust2
q2tils_qiime2_execute qiime picrust2 full-pipeline \
   --i-table $biom_file.qza \
   --i-seq hides-processing/rep-seqs-250f-250r-dada2.qza \
   --output-dir hides-processing/q2-picrust2_output \
   --p-threads 1 \
   --p-max-nsti 2 \
   --verbose

## Export picrust2 results.
q2tils_export hides-processing/q2-picrust2_output/pathway_abundance.qza feature-table.biom
mv pathway_abundance.biom hides-processing/q2-picrust2_output/pathway_abundance.biom
q2tils_qiime2_execute biom convert --to-tsv -i hides-processing/q2-picrust2_output/pathway_abundance.biom -o hides-processing/q2-picrust2_output/pathway_abundance.tsv

q2tils_export hides-processing/q2-picrust2_output/ko_metagenome.qza feature-table.biom
mv ko_metagenome.biom hides-processing/q2-picrust2_output/ko_metagenome.biom
q2tils_qiime2_execute biom convert --to-tsv -i hides-processing/q2-picrust2_output/ko_metagenome.biom -o hides-processing/q2-picrust2_output/ko_metagenome.tsv

q2tils_export hides-processing/q2-picrust2_output/ec_metagenome.qza feature-table.biom
mv ec_metagenome.biom hides-processing/q2-picrust2_output/ec_metagenome.biom
q2tils_qiime2_execute biom convert --to-tsv -i hides-processing/q2-picrust2_output/ec_metagenome.biom -o hides-processing/q2-picrust2_output/ec_metagenome.tsv

## Plot picrust2 results, pathogenesis focus.
Rscript analyze-picrust2-vf.R
