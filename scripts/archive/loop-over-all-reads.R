#!/home/tuba-nobackup/shared/R/R-3.3.0/bin/R

print(R.Version()$version.string)

library(Biostrings)
library(magrittr)
library(stringr)
library(dplyr)
library(readr)
library(tibble)
library(stringi)


rc <- function(x) { x %>% Hmisc::translate("ACGT", "TGCA") %>% stri_reverse }


TEMPLATESWITCH <- "TTTCTTATATGGG"
PRIMER_R4 <- "CGAGTACCATGGGCGTAAC"
ANNEAL <- "GTGTGACCTTCCCCAAAAGGCGTAG"
FWD_PRIMER <- "GAAGTTCCAGCCAGCGTC"
REV_PRIMER <- paste0(TEMPLATESWITCH, PRIMER_R4)

epia <- readDNAStringSet(filepath = "data/exp3_MHC/barcode_library/a.fa")
epib <- readDNAStringSet(filepath = "data/exp3_MHC/barcode_library/b.fa")
samp <- readDNAStringSet(filepath = "data/exp3_MHC/barcode_library/sample.fa")
#FWD_PRIMER <- readDNAStringSet(filepath = "data/exp3_MHC/barcode_library/forward_primer.a.fa")

BARCODE <- list(
  SAMPLE = samp,
  A_PRIM = FWD_PRIMER,
  A_N6 = "NNNNNN",
  A_OLIGO = epia,
  ANNEAL = ANNEAL,
  B_OLIGO = reverseComplement(epib),
  B_N6 = "NNNNNN",
  B_PRIM = rc(REV_PRIMER)
)



findBarcodes <- function(db, x) {
  
  plain_x <- class(x)!="QualityScaledDNAStringSet"
  plain_db <- class(db)!="DNAStringSet"
  N_db <- plain_db && str_detect(db[1] %>% as.character, "^N*$")
  
  i <- 1
  # sapply(seq_along(x), function(i) {
  # important!: seqs[[i]] will drop the quality encoding, so must use seqs[i] ...!
  if(plain_x) 
    subj <- x[[i]]
  else
    subj <- x[i]
  
  if(!N_db) {      
    pa <- pairwiseAlignment(
      db, subj,
      type = "global-local",
      gapOpening = 4, gapExtension = 1)
  } else {
    pa <- NA
  }
  
  if(plain_db) {
    hit <- pa
    nextScore <-  NA
    hit_name <- NA
  } else {
    sc <- score(pa)
    isMax <- sc == max(sc)
    
    if(sum(isMax) != 1)
      return( list(c(read = names(subj), NA, NA, NA, NA, NA, NA, NA)) )
    
    hit <- pa[isMax]
    
    # next best
    sc <- sc[!isMax]
    isNext <- sc == max(sc)
    nextScore <- sc[isNext][1]
    hit_name  <- names(alignedPattern(hit))
  }
  
  if(N_db) {
    hit_start <- NA
    hit_end   <- NA
    hit_len <- nchar(db[1] %>% as.character())
    hit_score <- NA
  } else {
    hit_start <- start(subject(hit))
    hit_end   <- end(subject(hit))
    hit_len <- nchar(hit)
    hit_score <- score(hit)
  }   
  
  # nchar(hit)
  
  return(  list(
    read = names(subj),
    start = hit_start,
    end = hit_end,
    score = hit_score,
    next_best = nextScore,
    length = hit_len,
    hit = hit_name,
    actual_length = hit_end - hit_start + 1)
  ) 
  
  # }) %>% do.call(rbind, .) 
  
}


# update start and end for the N seqs, based on expected based on length
fixNpositions <- function(x) {
  for (i in is.na(x$start) %>% which) {
    
    if(i == 1 || i == nrow(x)) { next }
    
    x$start[i] <- x$start[i+1] - x$length[i]
    x$end[i]   <- x$end[i-1]   + x$length[i]
    observed_length <- x$end[i] - x$start[i] + 1
    
    if(observed_length != x$length[i]) {
      previous_seq <- x$score[i-1] / x$actual_length[i-1]
      next_seq     <- x$score[i+1] / x$actual_length[i+1]
      
      # overwrite start and end based on the coords found in j
      if(previous_seq >= next_seq) {
        j <- i-1 
        x$start[i] <- x$end[j] + 1
        x$end[i]   <- x$end[j] + x$length[i]
      } else {
        j <- i+1 
        x$end[i]   <- x$start[j] - 1
        x$start[i] <- x$start[j] - x$length[i]
      }
      observed_length[i] <- x$end[i] - x$start[i] + 1
    }      
  }
  return(x)
}


# cat tmp/barracoda_in/IONTORRENT.R1.gems.no_umi.no_adapters.revcomp.fq | paste - - - - | grep GTGTGACCTT | tr '\t' '\n' > grep-anneal-reads.fq
# cat tmp/barracoda_in/IONTORRENT.R1.gems.no_umi.no_adapters.revcomp.fq | paste - - - - | grep -v GTGTGACCTT | tr '\t' '\n' > grep-anneal-reads-inverse.fq
# cat grep-anneal-reads-inverse.fq | paste - - - - | grep AAAGGCGTAG | tr '\t' '\n' > grep-anneal-reads-inverse-grep-end-of-anneal.fq
# cat grep-anneal-reads.fq grep-anneal-reads-inverse-grep-end-of-anneal.fq > grep-anneal-reads-both-ends.fq

nlines <- system("cat data/exp3_MHC/processed/tmp/grep-anneal-reads-both-ends.fq | wc -l", intern = TRUE) %>% as.numeric
nreads <- nlines/4

n <- 100
from <- seq(1, nreads, by= n)
to   <- c(seq(n, nreads, by= n), nreads)

for (i in seq_along(from)) {
 message("Analysing reads ", from[i], " to ", to[i])
   seq_i <- readQualityScaledDNAStringSet("data/exp3_MHC/tmp/grep-anneal-reads-both-ends.fq", skip = from[i]-1, nrec = n)
  
  df <- lapply(seq_along(seq_i), function(s) {
    seq_length <- nchar(as.character(seq_i[s]))
    aligned <- lapply(BARCODE, findBarcodes, seq_i[s]) %>% bind_rows(.id = "seq")
    aligned %<>% fixNpositions
    aligned %<>% mutate(dist_to_previous = start - c(0, end[-length(end)]) - 1,
                        dist_to_next     = c(start[-1], seq_length +1) - end - 1,
                        seq_length = seq_length)
  }) %>% bind_rows()
  
  # df %<>% mutate(seq = factor(seq, levels = names(BARCODE)))
  
  df %<>% rowwise() %>%  mutate(dist_to_best = min(abs(dist_to_previous), abs(dist_to_next)))
  
  write_tsv(df, path = "data/exp3_MHC/tmp/grep-anneal-reads.analysed.tsv", append = TRUE)
  
}
  
