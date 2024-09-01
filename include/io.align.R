## I wrote a set of fucntions to read and write alignments because most of the
## libraries out there are massive and I don't need them to properly read and
## write fasta/f*a files.

read.align = function(fname, type = 'list', truncate_at = Inf,
    comment.char = ';', nthreads=NULL)
## function read.align
## Desc: Reads in the fasta file [fname] and parses it into the type [type].
##       works as fast as I could get it to. There are improvements that could
##       be made if the user had prior knowledge about the structure of the 
##       alignment, but this will read in ANY fasta file and do so very fast.
## Args: fname = A filename to read in.
##       type  = the type that the sequences should be cast to (list, df, etc.)
##               Defaults to 'list'. Can also handle 'vector' and 'data.frame'.
##       truncate_at = Remove everything after n chars? -1 will truncate at the
##       first space.
##       comment.char = Ignore lines starting with this character. Default: ";"
##       nthreads = Number of threads to use with data.table::fread
{
    ## Read in data with fread (fastest reading, saves a LOT of time.)
    defthreads = data.table::getDTthreads()
    if (!is.null(nthreads))
    {
        data.table::setDTthreads(nthreads)
    }
    
    txt = data.table::fread(fname, sep='', stringsAsFactors=FALSE,
        data.table=FALSE, header=FALSE, check.names=FALSE)[,1]
    
    data.table::setDTthreads(defthreads)

    ## Remove comment lines so they are not parsed.
    is.comment = stringr::str_sub(txt,1,1) == ';'
    txt = txt[!is.comment]
    
    ## Find label positions. Is the first character a ">"?
    is.label = stringr::str_sub(txt,1,1) == '>'
    is.label[which(is.na(is.label))] = FALSE
    txt[which(is.na(is.label))] = ''
    seqnames = stringr::str_sub(txt[is.label], 2)
    
    
    ## Throw out the names if truncate at is 0.
    if (truncate_at == 0)
    {
        seqnames = NULL
    }
    else if (is.finite(truncate_at))
    {
        if (truncate_at == -1)
        {
            ## Strip the first space and everything after it.
            seqnames = stringr::str_remove(seqnames, ' .*')
        }
        
        ## Truncate sequence names to n letters
        seqnames = stringi::stri_sub(seqnames, from = 1, to = truncate_at)
    }
    
    ## Find lines that are sequences.
    seqlines = txt[!is.label]
    
    ## Find which index each sequence part belongs to
    index = cumsum(is.label)[!is.label]
    
    ## Separate sequences
    seqs = tapply(seqlines, index, stringr::str_c, collapse='')
    
    ## Remove spaces from inside sequences
    seqs = stringr::str_remove_all(seqs, ' ')
    
    ## Give sequences names
    names(seqs) = seqnames
    
    ## Cast sequences as the type we are interested in.
    ## as.DNAString or as.RNAString et cetera from biostrings might work here.
    if (type == 'vector')
    {
        return(seqs)
    } else {
        ## Cast return as the user defined type.
        return(get(paste0('as.',type))(seqs))
    }
}

## Dice a string into <width>-length substrings
str_dice = function(s, width)
{
   L = stringr::str_length(s)
   stringr::str_sub(s, start=seq(1L, L-width+1, width),
                       end=seq(width, L, width))
}

write.align = function(x, fname, colw = Inf, append = FALSE, nthreads = NULL)
## function write.align
## Desc: Writes the (flexible) object stored in x to the filename fname.
##       might work fast, but I haven't really done any writing of large fasta
##       files yet. Not optimized for small fasta files.
{
    ## Cast seqs as a list
    z = as.list(x)
    
    ## Fix the names.
    if (is.null(names(x)))
    {
        names(z) = paste0('>',rownames(x))
    } else {
        names(z) = paste0('>',names(x))
    }
    
    ## Split on column width.
    if (is.finite(colw))
    {
        z = lapply(z, str_dice, width = colw)
    }

    ## Add sequences to the top of the vector
    z = lapply(1:length(z), function(i) return(c(names(z)[i], z[[i]])))
    
    ## Collapse into a single matrix
    outz = matrix(unlist(z), ncol=1)
    
    ## fwrite to file.
    defthreads = data.table::getDTthreads()
    if (!is.null(nthreads))
    {
        data.table::setDTthreads(nthreads)
    }
    data.table::fwrite(outz, file = fname, append = append, quote = FALSE,
                       row.names = FALSE, col.names = FALSE, eol='\n')
    data.table::setDTthreads(defthreads)
}


read.fastq.dt = function(fname, nthreads=NULL)
{
    defthreads = data.table::getDTthreads()
    if (!is.null(nthreads))
    {
        data.table::setDTthreads(nthreads)
    }
    
    txt = data.table::fread(fname, sep='', stringsAsFactors=FALSE,
        header=FALSE, check.names=FALSE)
    l = nrow(txt)
    
    dats = txt[,.(name=V1[1:l %% 4 == 1], seq=V1[1:l %% 4 == 2],
                  addt = V1[1:l %% 4 == 3], qual=V1[1:l%%4 == 0])]
    data.table::setDTthreads(defthreads)
    
    return(dats)
}

write.fastq.dt = function(x, fname)
{
    z = x[,.(name, seq, addt, qual)]
    data.table::fwrite(z, fname, sep='\n', col.names=FALSE, compress='gzip')
}
