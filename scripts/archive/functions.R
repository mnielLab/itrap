library(readr)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(ggplot2)
library(magrittr)
library(Hmisc)
library(stringi)
library(purrr)
library(seqinr)

# ----- FUNCTIONS  -----------------

read_fasta <- function(f, patternToRemoveFromName = ".*-") {
  read.fasta(f,
             as.string = TRUE,
             strip.desc = TRUE,
             set.attributes = FALSE,
             forceDNAtolower = FALSE) %>%
    set_names( names(.) %>% str_remove(patternToRemoveFromName) )
}

rc <- function(x) { x %>% Hmisc::translate("ACGT", "TGCA") %>% stri_reverse }

kmer <- function(x,k,i) substr(x, i, i+k-1)


seq2kmer <- function(x, k) {
  n <- nchar(x)
  j <- seq(k, n)-k+1
  setNames(j, sapply(j, function(i) kmer(x, k, i) ))
}

uniqkmers <- function(x) {
  n <- names(x) %>% str_remove("\\*")
  isdup <- duplicated(n) | duplicated(n, fromLast = TRUE)
  x[!isdup]
}


kmermap <- function(kmer_x, kmer_p, length_p, return.all = FALSE) {
  
  names(kmer_p)[1] <- paste0("*", names(kmer_p)[1])
  names(kmer_p)[length(kmer_p)] <- paste0(names(kmer_p)[length(kmer_p)], "*")
  
  kmer_x <- uniqkmers(kmer_x)
  kmer_p <- uniqkmers(kmer_p)
  
  star_names <- names(kmer_p)
  names(kmer_p) <- str_remove(names(kmer_p), "\\*")
  
  starts <- sapply(names(kmer_p), function(i) unname( kmer_x[i] - kmer_p[i] + 1 ) )
  names(starts) <- star_names
  
  ends <- starts + length_p - 1
  if(return.all)
    return(cbind(starts, ends))
  
  if( all( is.na(starts) ) )
    from <- NA
  else
    # from <- min(starts, na.rm = TRUE)
    from <- starts[!is.na(starts)][1]
  
  if( all( is.na(ends)    ) ) 
    to <- NA
  else
    # to <- max(ends, na.rm = TRUE)
    to <- rev(ends[!is.na(ends)])[1]
  return(c(from, to))
  
}

# wrapper for kmermap
findprimer <- function(x, p, k = 8) {
  kmer_x <- seq2kmer(x, k)
  kmer_p <- seq2kmer(p, k)
  
  kmermap(kmer_x, kmer_p, nchar(p))
  
}


bestMatch <- function(x, db, k = 5) {
  n <- nchar(db[[1]])
  nx <- nchar(x)
  
  res <- lapply(db, function(p) {
    kmermap(
      seq2kmer(x, k),
      seq2kmer(p, k),
      nchar(p),
      return.all = TRUE) %>% t()
  })
  
  resfilt <- lapply(res, function(j) {
    # k <- nchar(colnames(j)[1])
    nas <- sum(is.na(j)) / 2
    # A fraction of kmers corresponding to max 3 mutation (each affecting k kmers),
    # but proportional to actual kmers in list
    # (some are removed if they were not unique)
    allowed_nas <- ncol(j) * 3 * k / (n-k+1)
    if(nas <= allowed_nas)
      return(j)
  })
  
  resfilt <- resfilt[sapply(resfilt, length)>0]
  
  maxKmerSupport <- sapply(resfilt, function(j) {
    lapply(as.data.frame(j), function(x) {
      if(any(is.na(x)))
        NULL
      else
        seq(x[1], x[2]) %>% tabulate(nbins = nx)
      } ) %>%
      do.call(rbind, .) %>% colSums %>% max

  })
  
  resfilt <- resfilt[order(maxKmerSupport, decreasing = TRUE)][1]
  
  if (length(resfilt) == 1) {
    
    hit <- resfilt[[1]]
  
    starts <- hit[1,]
    ends   <- hit[2,]
  
    from <- starts[!is.na(starts)][1]
    to <- rev(ends[!is.na(ends)])[1]
    
    out <- list(c(from, to))
    names(out) <- names(resfilt)
    
    # return( out )
    return( names(out) )
    
  } else {
    
    # return( c(NA, NA) )
    return( NA )
    
  }
  
}


kmerMatch <- function(x, p, k = 5) {
  n <- nchar(p)
  nx <- nchar(x)
  
  res <- kmermap(
      seq2kmer(x, k),
      seq2kmer(p, k),
      n,
      return.all = TRUE) %>% t()
  
  maxKmerSupport <- lapply(as.data.frame(res), function(x) {
      if(any(is.na(x)))
        NULL
      else
        seq(x[1], x[2]) %>% tabulate(nbins = nx)
      }) %>%
      do.call(rbind, .) %>% colSums %>% max

  return(maxKmerSupport)
  
}



# the absence of an asteirsk in 'to' means we know this number did not come from the last kmer, ie. we have "fuzzy" info on the ending of this mapping
