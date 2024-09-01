## Functions for sequence chemistry.

Tm = function(ntstring)
## Melting temperature of the sequences
{
    C = stringi::stri_count_fixed(ntstring, 'C')
    G = stringi::stri_count_fixed(ntstring, 'G')
    
    return(69.3 + (41 * (G + C) - 650) / nchar(ntstring))
}

Ta = function(fprimer, rprimer, amplicon)
## Annealing temperature of an amplicon as defined by forward primer, reverse
## primer, and intervening sequence
{
    return(0.3 * min(Tm(fprimer), Tm(rprimer)) + 0.7 * Tm(amplicon) - 14.9)
}

GC = function(ntstring)
## GC content
{
    C = stringi::stri_count_fixed(ntstring, 'C')
    G = stringi::stri_count_fixed(ntstring, 'G')
    
    return((G + C) / nchar(ntstring))
}

homopolymer = function(ntstring, homop = 3)
## Identify strings of repeated letters that are homop in length or greater.
{
    if(homop < 2) stop('Homopolymer must be > 1')
    
    ## Split sequences into letters.
    spl = strsplit(ntstring, '')
    
    count_repeats = function(x)
    {
        hmat = rbind(x, do.call('rbind', lapply(1:(homop-1),
            function(z) c(x[-1*(1:z)], rep('-',z)))))
        any(apply(hmat, 2, function(y) all(y == y[1])))
    }
    
    ## Detect substring homopolymers.
    cbm = unlist(lapply(spl, count_repeats))
    
    return(cbm)
}

