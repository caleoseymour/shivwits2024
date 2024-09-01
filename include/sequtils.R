## Collection of functions for dealing with sequence data.

#------------------------------------------------------------------------------#
## Note: I use a function lapply a lot. lapply is an iterator function that, in
## R, is faster than using a for loop. The format is lapply(x, function) where x
## is an teratable object and function is the function applied to the whole
## object. These two peices of code are equivalent:

## Using lapply:
#     o = lapply(x, func);

## Using a for loop:
#     o = list();
#     for(i in 1:length(x))
#     {
#          o[[i]] = func(x);    
#     }

## The latter is much faster. There are other iterators that are used as well.
## data table objects can be iterated using the following format:

#     d = data.table::data.table(x = <vector>, y = <vector>)
#     d[,o := fun(x)]

## This creates a column in d named o, then populates it with fun(x)

#------------------------------------------------------------------------------#

## Replace U's with T's.
rna2dna = function(seqs)
{
    s = stringr::str_replace_all(seqs, pattern = 'U', replacement = 'T')
    s = stringr::str_replace_all(s, pattern = 'u', 't')
    
    return(setNames(s, names(seqs)))
}

## Replace T's with U's.
dna2rna = function(seqs)
{
    s = stringr::str_replace_all(seqs, pattern = 'T', replacement = 'U')
    s = stringr::str_replace_all(s, pattern = 't', replacement = 'u')
    
    return(setNames(s, names(seqs)))
}

## Sanitize alignment, changing unknowns to gaps.
sanitize.align = function(seqs)
{
    s = toupper(seqs)
    s = stringr::str_replace_all(seqs, '.', '-')
    s = stringr::str_replace_all(seqs, '?', '-')
    
    return(s)
}

## Reverse-complement a sequence.
reverse_complement = function(ntstring)
{
    guide = c('A' = 'T',
              'T' = 'A',
              'G' = 'C',
              'C' = 'G')
    
    ## Split into letters
    spl = strsplit(rna2dna(ntstring), '')
    
    ## Replace letters using guide, reversing.
    rc = lapply(spl, function(x) paste0(rev(guide[x]), collapse=''))
    
    ## Combine again.
    return(do.call('c', rc))
}

## Count the length of a set of sequences
countsites.align = function(seqs) return(unlist(lapply(seqs, nchar)))

## Cast sequences as a data table (two columns, column 1 is "name", column 2
## is the sequence.)
seqs2dt = function(seqs)
{
    return(data.table::data.table(name = names(seqs),
        seq = unlist(unname(seqs))))
}

## Ensure that all sequences are the same length. Return lengths if they are.
samelen.align = function(seqs)
{
    seqlens = countsites.align(seqs)
    
    ## Make sure all sequences are the same length stop if not.
    if(!all(seqlens == seqlens[1]))
    {
        ## If any other sequences have the same length as the first, then
        ## we consider those the "good" ones. Otherwise, the first must be the
        ## offending sequence.
        if (any(seqlens)[-1] == seqlens[1])
        {
            offender = names(seqs)[seqlens != seqlens[1]]
        } else {
            offender = names(seqs)[1]
        }
        
        ## Build string listing the offenders.
        if (length(offender) == 1)
        {
            offendstring = offender
        } else {
            offendstring = paste0(offender, ' and ', length(offender)-1,
                ' others.')
        }
        
        ## Error out
        err = paste0('Alignment length differs. Offending sequence(s):',
            offendstring)
        stop(err)
    }
    
    ## Otherwise, return sequence length
    return(seqlens[1])
}

## Remove all gaps from a sequence object
dealign.align = function(seqs, gapchar = '-', type=typeof(seqs))
{
    ## Use this to cast the final product.
    newtypefun = get(paste0('as.',type))
    
    ## cast as a data table
    seqdt = seqs2dt(seqs)
    
    ## replace sites
    for(s in gapchar)
    {
        seqdt[,seq := stringr::str_remove_all(seq, s)]
    }
    
    ## Assign names to object
    sout = setNames(newtypefun(seqdt[, seq]), seqdt[, name])
    
    return(sout)
}

## This function takes an alignment and returns the position of the aligned
## residue in the non-aligned sequence.
profile.align = function(seqs, gapchar = '-')
{
    seql = samelen.align(seqs)[1]
    
    ## cast as a data table
    seqdt = seqs2dt(seqs)
    
    ## Split into matrix
    ssplitdt = seqdt[, .(residue = unlist(strsplit(seq,'')),
        in.alignment = 1:seql), by = name]
    
    ## Identify 
    ssplitdt[, in.sequence := cumsum(!(residue %in% gapchar)), by = name]
    ssplitdt[, keep := !duplicated(in.sequence) & in.sequence > 0, by = name]
    
    sout = ssplitdt[keep == TRUE,.(name, residue, in.alignment, in.sequence)]
    return(sout)
}

## Remove gap-only sites using data.table.
reduce.align = function(seqs, gapchar = '-', quiet=FALSE, type = typeof(seqs))
{   
    seql = samelen.align(seqs)[1]
    
    ## cast as a data table
    seqdt = seqs2dt(seqs)
    
    ## get individual sites as a two-column table.
    ssplitdt = seqdt[, .(residue = unlist(strsplit(seq, '')),
        site = 1:seql), by = name]

    ## Identify sites that have data that aren't gaps by removing all
    ## gap sites, then returning a unique vector of site numbers.
    ngapsites = ssplitdt[!(residue %in% gapchar), sort(unique(site))]

    ## Report the number of sites removed from the final alignment
    if (!quiet)
    {
        message(paste0(length(ngapsites), '/', seql, ' sites retained.'))
    }

    ## Extract non-gap sites and paste them together as sequences.
    spastedt = ssplitdt[site %in% ngapsites,
        .(seq = paste0(residue, collapse='')), by = name]
        
    ## Cast as user defined type and return
    sout = setNames(get(paste0('as.', type))(spastedt[,seq]), spastedt[,name])
    return(sout)
}

## Get the consensus from a list of sequences
consensus = function(seqs, maxdegen = 3)
{
    bases = c('A','G','T','C')
    ntmatrix = do.call('rbind', strsplit(seqs, ''))
    site_proportions = apply(ntmatrix, 2,
        function(x) table(factor(x, levels = bases))) / length(seqs)
    passing_consensus = site_proportions >= 0.25
    
    degens = sum(colSums(passing_consensus) > 1)
    if (degens > maxdegen)
    {
        return(NA)
    }
    
    consensus_bases = apply(passing_consensus, 2, function(x) bases[x])
    #print(consensus_bases)
    if(is.list(consensus_bases))
    {
        consensus_matrix = do.call('expand.grid', consensus_bases)
        consensus_probes = apply(consensus_matrix, 1, paste0, collapse='')
    } else {
        consensus_probes = paste0(consensus_bases, collapse = '')
    }
    
    return(consensus_probes)
}   