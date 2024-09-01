#!/bin/bash

## Cale Seymour
## 2023

## A collection of utility functions for working with EMP-protocol Illumina
## MiSeq data. Uses Qiime2.

export Q2TILS_DEBUG=1   ## Debug (does nothing right now.)
export Q2TILS_REBUILD=0 ## Force recreation of files that already exist.
export Q2TILS_FORWARD_TRUNC='' ## Forward read truncation default.
export Q2TILS_REVERSE_TRUNC='' ## Reverse read truncation default.
export Q2TILS_PROCESSORS=1 ## Number of processors to use for high-level processes.


export Q2TILS_UID=$(tr -dc A-Za-z0-9 </dev/urandom | head -c 13 ; echo '') ## Unique identifier for this call.
export Q2TILS_LOCATION="$( dirname -- "${BASH_SOURCE[0]}" )" ## Location from which the script was run.
export Q2TILS_LOGFILE="$Q2TILS_UID.logfile.txt"
## Basic bash functions.
stop_on_error ()
## Function that kills the script if you provide it with an input of anything
## except 0 as the first argument. Can be used to interpret error codes,
## since "0" is considered to be a success.
{
    if [ $1 -ne 0 ];
    then
        echo "Last commmand exited with nonzero error code: $1"
        exit
    fi
    return 0
}

printer()
## Logging function. Match arguments and print.
{
    for i in "$*"
    do 
        printf "$i " >> $Q2TILS_LOGFILE
        printf "$i "
    done
    printf "\n" >> $Q2TILS_LOGFILE
    printf "\n"
    return 0
}

logger()
## Logging function. Match arguments and print with timestamp.
{
    printf "[@$(date +'%m/%d/%Y %H:%M:%S')]: "
    printer "$*"
    return 0
}

debugger()
## Debug function. Match arguments and print with timestamp if the debug flag == 1.
{
    [ $Q2TILS_DEBUG -eq 1 ] && logger $*
    return 0
}

make_dir_if_not_exists()
## Make a directory if it doesn't already exist.
{
    [[ ! -d $1 ]] && mkdir $1
    return $?
}

q2tils_qiime2_execute()
## This function runs a qiime2 command in the selected conda env.
## All arguments are passed to the conda environment within quotes.
{
    ## Build the string to execute.
    exec_str=""
    for i in $*
    do 
        exec_str="$exec_str $i"
    done
    
    printer "... Running cmd: '$exec_str' (selected conda env: $Q2TILS_QIIME_ENV)"
    conda run --no-capture-output -n $Q2TILS_QIIME_ENV bash -c "$exec_str" &>> $Q2TILS_LOGFILE
}

q2tils_mothur_execute()
## Execute a command using mothur.
{
    exec_str="set.logfile(name=$Q2TILS_LOGFILE, append=T);"
    for i in $*
    do 
        exec_str="$exec_str $i"
    done
    
    printer "... Running cmd: '$exec_str' using mothur $Q2TILS_MOTHUR_VERSION"
    mothur "#$exec_str" >/dev/null
}

q2tils_mothur_align_procedure()
## Align reads to the SILVA reference DB. Requires most other steps to have been run already.
{
    infile=$1
    taxonomy=$2
    cp $infile $Q2TILS_UID.temp.fasta
    output="${infile%.*}.filter.fasta"
    
    logger 'Doing some work in mothur to generate a good alignment...'
    printer "... Aligning sequences in '$infile'..."
    
    [ ! -f  $Q2TILS_UID.temp.align ] && q2tils_mothur_execute "align.seqs(candidate=$Q2TILS_UID.temp.fasta, template=$Q2TILS_LOCATION/silva.nr_v138_1.515f_806r.align, processors=$Q2TILS_PROCESSORS, flip=F)"
    [ ! -f  $Q2TILS_UID.temp.align ] && return 1
    
    printer "... Removing domain-level unassigned sequences based on taxonomy in '$taxonomy'..."
    if [ ! -z $taxonomy ] && [ -f $taxonomy ];
    then 
        grep -P '\tUnassigned\t' $taxonomy | cut -f 1 > $Q2TILS_UID.temp.unassigned.accnos
        grep -P '\td__Eukaryota;' $taxonomy | cut -f 1 >> $Q2TILS_UID.temp.unassigned.accnos
        q2tils_mothur_execute "remove.seqs(accnos=$Q2TILS_UID.temp.unassigned.accnos, fasta=$Q2TILS_UID.temp.align)"
    fi
    
    printer "... Filtering vertical gaps from the final alignment..."
    [ ! -f  $Q2TILS_UID.temp.filter.pick.fasta ] && q2tils_mothur_execute "filter.seqs(fasta=$Q2TILS_UID.temp.pick.align, vertical=T, processors=$Q2TILS_PROCESSORS)"
    [ ! -f  $Q2TILS_UID.temp.filter.pick.fasta ] && return 1
    
    mv $Q2TILS_UID.temp.filter.pick.fasta $output
    printer "... Filtered alignment '$output'."
}

q2tils_export()
## Export the file $2 from the qza object $1
{
    term=$2
    inf=$1
    ext=$(echo $term | sed -e 's/.*\.//g')
    outf=$(basename -s .qza $inf).$ext
    zipath=$(unzip -l $inf | grep "$term" | sed 's/.* //g')
    unzip -p $inf $zipath > $outf
}

q2tils_make_manifest()
## This function attempts to guess the type of data that we need to import.
## $1 is the directory that contains the reads.
## $2 is the sample metadata file. First column has to be the sample name.
## Basically, we iterate over the first column of the sample metadata, and then
## we attempt to match each sample name to a read file. If we match two read
## files, then we switch to paired mode.
{
    dirname=$1
    metadata=$2
    
    ## Default directory is current dir.
    [ -z $1 ] && dirname=$PWD
    
    logger "Making a manifest in the dir '$dirname'..."
    
    ## Count the number of readfiles in the directory.
    nreads=$(find $dirname -name "*.fastq*" -exec printf %c {} + | wc -c)
    
    if [ $nreads -eq 0 ];
    then
        printer "... No reads were found in the directory '$dirname'."
        return 1
    fi
    
    ## If no metadata file exists, then we're going to make a fake one by
    ## stripping the locus tags from the read files.
    if [ -z $2 ];
    then
        echo "#SampleID" > $dirname/.metadata.txt
        metadata=$dirname/.metadata.txt
        
        printer "... A metadata file was not provided. We're going to try to guess the sample names by looking for forward reads."
        find $dirname -name "*.fastq*" -print | \
            grep "_S[0-9]_L[0-9]*_R1_[0-9]*\.fastq" | \
            sed -e "s/.*\///g" | \
            sed -e 's/_S[0-9]_L[0-9]*_R1_[0-9]*\.fastq.*$//g' \
        >> $dirname/.metadata.txt
    fi
    printer "... Using samples from the metadata file, '$metadata'."

    ## Create some manifest files.
    echo sample-id,absolute-filepath,direction > $dirname/.manifest-paired.txt
    echo sample-id,absolute-filepath,direction > $dirname/.manifest-forward.txt
    
    ## Pre-set the mode that we're going to use.
    mode="paired"
    
    ## Read the metadata file and iterate over samples.
    for sampleid in $(grep -v '^#' $metadata | awk 'BEGIN { FS="\t" } { print $1 }')
    do
        sanitizedid=$(echo $sampleid | sed -e 's/_/-/g')
        if [ $sanitizedid != $sampleid ];
        then
            printer "... Underscores in $sampleid were replaced with dashes: the new sample name is '$sanitizedid'."
        fi
        
        nreads=$(find $dirname -name "$sampleid*.fastq.gz" -exec printf %c {} + | wc -c)
        if [ $nreads -eq 0 ];
        then
            printer "... ERROR: No reads were found for the sample '$sampleid'. Exiting."
            return 1
        elif [ $nreads -eq 2 ];
        then
            ## Todo: add support for multi-run reads.
            forwardReads=$(find $dirname -name "$sampleid*_R1_*.fastq*" | xargs -0 realpath)
            echo "$sanitizedid,$forwardReads,forward" >> $dirname/.manifest-paired.txt
            echo "$sanitizedid,$forwardReads,forward" >> $dirname/.manifest-forward.txt
            
            reverseReads=$(find $dirname -name "$sampleid*_R2_*.fastq*" | xargs -0 realpath)
            echo "$sampleid,$reverseReads,reverse" >> $dirname/.manifest-paired.txt

        elif [ $nreads -eq 1 ];
        then
            [ $mode = "paired" ] && printer "... WARNING: Only one read file was found for the sample $sampleid. Switching to forward-only mode."
            mode="forward"
            forwardReads=$(find $dirname -name "$sampleid*.fastq*" -print)
        else
            printer "... ERROR: More than two matching read files were found for the sample '$sampleid'. Exiting."
            return 1
        fi
    done
    
    ## Clean up the manifests
    mv $dirname/.manifest-$mode.txt $dirname/manifest.txt
    rm $dirname/.manifest-*.txt
    
    printer "... Made a $mode manifest in the directory '$dirname'."
    return 0
}

q2tils_check_manifest()
## This function checks to see if a manifest file is going to be OK.
## $1 is the directory with the reads
## $2 is the manifest file
## $3 is the sample metadata file.
## In essence, all we want to do is check to see if all three are fully 
## concordant and return 0 if so.
{
    dirname=$1
    manifest=$2
    metadata=$3
    
    logger "Checking the manifest, '$manifest' for problems..."

    
    ## Collect validate the header and the csv format.
    [[ $(head -n 1 $manifest) -ne "sample-id,absolute-filepath,direction" ]] && printer "... Found a problem with the manifest: Malformed header. The first line should be 'sample-id,absolute-filepath,direction'."
    
    ## Make sure each line has three fields, and make sure the read files exist.
    tail +2 $manifest | while IFS=, read sampleid absolutefilepath direction; do
        if [[ -z $sampleid ]];
        then
            printer "... Found a problem with the manifest: Missing sample id."
        fi
        
        if [[ ! -f  $absolutefilepath ]];
        then
            printer "... Found a problem with the manifest: Read file '$absolutefilepath' does not exist."
            return 1
        fi
        if [[ $direction -ne "forward" ]] && [[ $direction -ne "reverse" ]];
        then
            printer "... Found a problem with the manifest: Direction must be one of 'forward' or 'reverse'."
            return 1
        fi
    done
    
    
    ## Next, we compare the number of read files to the number of samples.
    ## The idea is that we should either have 1 or 2 read files per sample;
    ## no more and no less.
    nreads=$(find $dirname -name "*.fastq*" -exec printf %c {} + | wc -c)
    nsamples=$(grep -vc '^#' $metadata)
    
    if [[ $nreads == $((nsamples*2)) ]];
    then
    ## Paired mode -- two readfiles per sample.
        for sampleid in $(grep -v '^#' $metadata | awk 'BEGIN { FS="\t" } { print $1 }')
        do
            sanitizedid=$(echo $sampleid | sed -e 's/_/-/g')
            if [[ $(grep -cF "$sanitizedid" $manifest) != 2 ]];
            then
                printer "... Found a problem with the manifest: There are not exactly two read files found for the sample '$sampleid'."
                return 1
            fi
        done
    elif [[ $nreads == $nsamples ]];
    then
    ## Single mode -- one readfile per sample expected.
        for sampleid in $(grep -v '^#' $metadata | awk 'BEGIN { FS="\t" } { print $1 }')
        do
            sanitizedid=$(echo $sampleid | sed -e 's/_/-/g')
            if [[ $(grep -c "$sanitizedid" $manifest) != 1 ]];
            then
                printer "... Found a problem with the manifest: There is not exactly one read file found for the sample '$sampleid'."
                return 1
            fi
        done
    else
        printer "... Found a problem with the manifest: There are $nreads read files and $nsamples samples. There should be exactly one or exactly two read files per sample."
        return 1
    fi
    
    for readfile in $(find $dirname -name "*.fastq*")
    do
        rf=$(basename $readfile)
        if [[ $(grep -cF "$rf" $manifest) != 1 ]];
        then
            printer "... Found a problem with the manifest: There is not exactly one sample id found for the readfile '$rf'."
            return 1
        fi
    done
    
    printer "... The manifest is OK."
    return 0
}

q2tils_import_reads()
## This is a catch-all function for importing fastq reads.
## $1 is the directory with the reads (and the manifest).
## $2 is the directory with the output files.
{
    ## This script does not currently support single-end or multiplexed reads.
    ## It also does not support dual-index reads where the index files are not
    ## yet merged. Need to add support for that in the future.
    dirname=$1
    output=$2
    
    if [[ -z $output ]];
    then
        output='.'
    else
        make_dir_if_not_exists $2
    fi

    logger "Importing demux-ed sequences to the directory '$output'..."
    if [[ ! -f $dirname/manifest.txt ]];
    then
        printer "...ERROR: Didn't find a manifest in '$dirname'. Run q2tils_make_manifest on the directory to build one."
        return 1
    fi
    
    if [[ ! -f $output/paired-end-demux.qza ]] || [[ $Q2TILS_REBUILD == 1 ]];
    then
        q2tils_qiime2_execute qiime tools import \
          --type 'SampleData[PairedEndSequencesWithQuality]' \
          --input-path $dirname/manifest.txt \
          --output-path $output/paired-end-demux.qza \
          --input-format PairedEndFastqManifestPhred33
        stop_on_error $?
    fi
    
    printer "... Imported reads to '$output/paired-end-demux.qza'."
}

q2tils_analyze_read_truncation()
## The demultiplexed reads are truncated based on a PHRED score threshold.
## This function identifies the appropriate read truncation length by looking
## at the seven-number summaries contained in each of the reverse or forward
## read files.
## $1 is a paired-end-demux file.
## $2 is the output directory.
{
    input=$1
    output=$2/$(basename $input | sed 's/\.qza$/.qzv/')
    
    logger "Setting truncation thresholds for demux: '$input'."
    
    make_dir_if_not_exists $2
    
    if [[ ! -f $output ]] || [[ $Q2TILS_REBUILD == 1 ]];
    then
        q2tils_qiime2_execute qiime demux summarize \
          --i-data $input \
          --o-visualization $output
        stop_on_error $?
    fi

    ## Parse the seven number summaries to find the position at which the median
    ## sequence phred score is < 20.
    
    q2tils_export $output forward-seven-number-summaries.tsv
    mv paired-end-demux.qzv.tsv $2/forward-seven-number-summaries.tsv
    

    col=0
    for i in $(grep "^50\%" $2/forward-seven-number-summaries.tsv | sed 's/.*\%\t//')
    do
        inti=$(echo $i | sed 's/\.[0-9]*$//')
        value=$((inti + 0))
        [[ $value < 20 ]] && break
        col=$((col + 1))
    done
    
    export Q2TILS_FORWARD_TRUNC=$(head -n 1 $2/forward-seven-number-summaries.tsv | awk -v x=$col '{print $x}')
    printer "... Forward trunc: $Q2TILS_FORWARD_TRUNC"

    q2tils_export $output reverse-seven-number-summaries.tsv
    mv paired-end-demux.qzv.tsv $2/reverse-seven-number-summaries.tsv
    
    col=0
    for i in $(grep "^50\%" $2/reverse-seven-number-summaries.tsv | sed 's/.*\%\t//')
    do
        inti=$(echo $i | sed 's/\.[0-9]*$//')
        value=$((inti + 0))
        [[ $value < 20 ]] && break
        col=$((col + 1))
    done
    
    export Q2TILS_REVERSE_TRUNC=$(head -n 1 $2/reverse-seven-number-summaries.tsv | awk -v x=$col '{print $x}')
    printer "... Reverse trunc: $Q2TILS_REVERSE_TRUNC"
    
    return 0
}

q2tils_dada2()
## Run dada2 using the stored truncaion vals.
## $1 is the input of paired-end demultiplexed reads.
## $2 is the output directory.
{
    logger "Executing DADA2..."
    input=$1
    output=$2
    
    if [[ -z $output ]];
    then
        output='.'
    else
        make_dir_if_not_exists $2
    fi
    
    ## Run integrity check and run command.
    if [[ ! -f $output/table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza ]] || \
       [[ ! -f $output/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza ]] || \
       [[ ! -f $output/dada2-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-stats.qza ]] || \
       [[ $Q2TILS_REBUILD == 1 ]];
    then
    
        logger "Running DADA2 on read file '$input'."
        printer "... Forward trunc: $Q2TILS_FORWARD_TRUNC"
        printer "... Reverse trunc: $Q2TILS_REVERSE_TRUNC"

        if [[ -z $Q2TILS_FORWARD_TRUNC ]] || [[ -z $Q2TILS_REVERSE_TRUNC ]];
        then
            printer "... Forward or reverse truncation values are invalid. Make sure to run 'q2tils_analyze_read_truncation' prior to this function, or set them manually."
            return 1
        fi
        
        q2tils_qiime2_execute qiime dada2 denoise-paired \
            --i-demultiplexed-seqs $input \
            --p-trunc-len-f $Q2TILS_FORWARD_TRUNC \
            --p-trunc-len-r $Q2TILS_REVERSE_TRUNC \
            --o-table $output/table-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza \
            --o-representative-sequences $output/rep-seqs-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza \
            --o-denoising-stats $output/dada2-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-stats.qza \
            --verbose
        stop_on_error $?
    fi
    
    printer "... Output dada2 files into $output."
    return 0
}

q2tils_classify_nb()
## Classify using a naiive-bayes classifier.
## $1 is the location of the input dada2 reads.
## $2 is the file to which the classification should be output.
## $3 is the classifier location. Defaults to the 515F/806R classifier that is
## built when this script is sourced.
{
    "Executing sequence classification..."
    input=$1
    output=$2
    classifier=$3
    
    if [[ -z $output ]];
    then
        output='.'
    else
        make_dir_if_not_exists $2
    fi
    
    ## Default to the 515f-806r classifier that gets auto-built.
    if [[ -z $classifier ]];
    then
        classifier=$Q2TILS_LOCATION/nb-classifier.$Q2TILS_QIIME_ENV.silva-138-99.515f-806r.qza
    fi
    
    if [ ! -f $output/taxonomy-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza ];
    then
        printer "... Classifying '$input' using the classifier '$classifier'..."
        q2tils_qiime2_execute qiime feature-classifier classify-sklearn \
            --i-classifier $classifier \
            --i-reads $input \
            --o-classification $output/taxonomy-"$Q2TILS_FORWARD_TRUNC"f-"$Q2TILS_REVERSE_TRUNC"r-dada2.qza \
            --p-reads-per-batch 1000
        stop_on_error $?
    fi
    printer "... Output taxonomy into $output."
    return 0
}

printer "... logfile for this run is '$Q2TILS_LOGFILE'."

## Grab the most recent qiime2 version and save it to an environment var
logger "Checking pipeline dependencies..."
if [ -z $(which conda) ];
then
    printer "... no 'conda' distribution was found. Please install one. It doesn't matter which one you use as long as you can call it via the command 'conda'."
    return 1
fi

printer "... Using conda: '$(which conda)'"

export Q2TILS_QIIME_ENV=$(conda env list | grep "^qiime2-[0-9]*\.[0-9]" | tail -n 1 | sed 's/ .*//g')
stop_on_error $?

if [ -z $Q2TILS_QIIME_ENV ];
then
    printer "... I didn't find a qiime2 conda environment. Make sure you follow the instructions for installation via miniconda."
    exit
fi

printer "... Using qiime2 conda env: '$Q2TILS_QIIME_ENV'."

## Check for the sequence classifier base arifacts.
if [[ ! -f $Q2TILS_LOCATION/nb-classifier.$Q2TILS_QIIME_ENV.silva-138-99.515f-806r.qza ]] ;
then
    QIIME_ENV_YEAR=$(echo $Q2TILS_QIIME_ENV | sed 's/^qiime2-//')
    logger "You don't have a 515F/806R sequence classifier!"
    printer "... Fetching the seqs and taxonomy for a 515F/806R classifier from the qiime2 folks. Make sure to cite Bokulich et al. 2020: https://doi.org/10.1101/2020.10.05.326504"
     
    wget -qO $Q2TILS_LOCATION/silva-138-99-seqs-515-806.qza https://data.qiime2.org/"$QIIME_ENV_YEAR"/common/silva-138-99-seqs-515-806.qza
    stop_on_error $?
    
    wget -qO $Q2TILS_LOCATION/silva-138-99-tax-515-806.qza https://data.qiime2.org/"$QIIME_ENV_YEAR"/common/silva-138-99-tax-515-806.qza
    stop_on_error $?
    
    printer "... Building the qiime2 sequence classifier. The location of the classifier will be: '$Q2TILS_LOCATION/nb-classifier.$Q2TILS_QIIME_ENV.silva-138-99.515f-806r.qza'"
    q2tils_qiime2_execute qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads $Q2TILS_LOCATION/silva-138-99-seqs-515-806.qza \
      --i-reference-taxonomy $Q2TILS_LOCATION/silva-138-99-tax-515-806.qza \
      --o-classifier $Q2TILS_LOCATION/nb-classifier.$Q2TILS_QIIME_ENV.silva-138-99.515f-806r.qza \
      --p-classify--chunk-size 8000
fi

if [ -z $(which mothur) ];
then
    printer "... I couldn't find a mothur executable. The pipeline will still work, but mothur required to align the sequences to the SILVA reference."
else
    export Q2TILS_MOTHUR_VERSION=$(mothur -v | grep 'version=' | sed 's/.*version=//')    


    if [ ! -f $Q2TILS_LOCATION/lane1349.silva.filter ];
    then
        printer "... Downloading the SILVA Lane mask filter for mothur..."
        wget -qO $Q2TILS_LOCATION/lane1349.silva.filter https://mothur.s3.us-east-2.amazonaws.com/wiki/lane1349.silva.filter
    fi
    if [ ! -f $Q2TILS_LOCATION/silva.nr_v138_1.515f_806r.align ];
    then
        printer "... You don't have a mothur-compatible 515f/806R SILVA reference alignment. I'm going to build one."
        wget -qO $Q2TILS_LOCATION/silva.nr_v138_1.tgz https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_1.tgz
        tar -xf $Q2TILS_LOCATION/silva.nr_v138_1.tgz silva.nr_v138_1.align -O > $Q2TILS_LOCATION/silva.nr_v138_1.align
        q2tils_mothur_execute "pcr.seqs(fasta=$Q2TILS_LOCATION/silva.nr_v138_1.align, start=11895, end=25318, keepdots=F, processors=$Q2TILS_PROCESSORS)"
        rm $Q2TILS_LOCATION/silva.nr_v138_1.align
        rm $Q2TILS_LOCATION/silva.nr_v138_1.scrap.pcr.align
        rm $Q2TILS_LOCATION/silva.nr_v138_1.tgz
        q2tils_mothur_execute "screen.seqs(fasta=$Q2TILS_LOCATION/silva.nr_v138_1.pcr.align, start=1, end=25318, processors=$Q2TILS_PROCESSORS)"
        mv $Q2TILS_LOCATION/silva.nr_v138_1.pcr.align $Q2TILS_LOCATION/silva.nr_v138_1.515f_806r.align
        stop_on_error $?
        
        printer "... 515f/806R SILVA reference alignment is built."
    fi
fi

## Check for the sequence classifier.
# if [[ ! -f $Q2TILS_LOCATION/silva-138-99-nb-weighted-classifier.qza ]] ;
# then
     # printer "You don't have a 515F/806R sequence classifier!"
     # printer "... Fetching a 515F/806R classifier from the qiime2 folks. Make sure to cite Bokulich et al. 2020: https://doi.org/10.1101/2020.10.05.326504"

     # wget -qO $Q2TILS_LOCATION/silva-138-99-seqs-515-806.qza https://data.qiime2.org/"$(echo $Q2TILS_QIIME_ENV | sed 's/^qiime2-//')"/common/silva-138-99-nb-weighted-classifier.qza
     # stop_on_error $?
# fi
