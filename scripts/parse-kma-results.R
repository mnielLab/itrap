#!/home/tuba-nobackup/shared/R/R-3.6.1/bin/Rscript 

#################
### LIBRAREIS ###
#################

library(Biostrings)
library(magrittr)
library(stringr)
library(dplyr)
library(readr)
library(tibble)
library(stringi)
library(parallel)

################################
###  COMMAND LINE ARGUMENTS  ###
################################

# RUN
#/home/tuba-nobackup/shared/R/R-3.6.1/bin/Rscript parse-kma-results.R ../data/helle/mapping.frag ../data/barcode_library/sample.fa ../data/barcode_library/a.fa ../data/barcode_library/b.fa ../data/barcode_library/barcode-information.fa ../results/kma_parser_parallel.tsv

cArgs <- commandArgs(TRUE)

# Output from KMA
kma_file     <- cArgs[1]

# Barcode information
sample_file  <- cArgs[2]
oligo_A_file <- cArgs[3]
oligo_B_file <- cArgs[4]
barcode_file <- cArgs[5]

# Output filename
out_file     <- cArgs[6]

# Check if the correct number of input arguments were given
if (length(cArgs) < 6 | length(cArgs) > 6) {
  stop("The number of input arguments is not correct. There should be 6 input arguments! 
       Input arguments are: kma_file sample_file oligo_A_file oligo_B_file barcode_file out_file")
}

#################
### FUNCTIONS ###
#################

# reverse complement function 
rc <- function(x) { x %>% Hmisc::translate("ACGT", "TGCA") %>% stri_reverse }

# findBarcodes function makes a pairwise alignment to find 
# the different elements (eg. sample, primers, UMIs, and oligos) in the sequence   

findBarcodes <- function(db, x) {
  
  plain_x <- class(x)!="QualityScaledDNAStringSet"
  plain_db <- class(db)!="DNAStringSet"
  N_db <- plain_db && str_detect(db[1] %>% as.character, "^N*$")
  
  i <- 1
  # important!: seqs[[i]] will drop the quality encoding, so must use seqs[i] ...!
  if(plain_x) {
    subj <- x[[i]]
    names(subj) <- names(x)
  } else {
    subj <- x[i]
  }
  
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
      return( list(read = names(subj), start = NA, end = NA, score = NA, next_best = NA, length = NA, hit = NA, actual_length = NA) )
    
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
  
}

# The fixNpositions will update start and end for the N seqs, 
# based on expected based on length

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


fx <- function(kma_hits) {
  results <- sapply(1:nrow(kma_hits), function(i) {
    
    # message("Analysing reads ", from[i], " to ", to[i])
    seq_i <- rc(kma_hits[i,"seq"])
    names(seq_i) <- kma_hits[i,"read"]
    
    s <- 1
    seq_length <- nchar(as.character(seq_i[s]))
    aligned <- lapply(set_names(names(BARCODE), names(BARCODE)), function(ELEM) {
      if(ELEM=="A_OLIGO") {
        A <- kma_hits[i, "A"]
        DB <- BARCODE$A_OLIGO[A]
      } else if(ELEM=="B_OLIGO") {
        B <- kma_hits[i, "B"]
        DB <- BARCODE$B_OLIGO[B]
      } else {
        DB <- BARCODE[[ELEM]]
      }
      findBarcodes(DB, seq_i)
    }) %>% bind_rows(.id = "seq")
    
    aligned %<>% fixNpositions
    aligned %<>% mutate(dist_to_previous = start - c(0, end[-length(end)]) - 1,
                        dist_to_next     = c(start[-1], seq_length +1) - end - 1,
                        seq_length = seq_length)
    aligned %<>% rowwise() %>%  mutate(dist_to_best = min(abs(dist_to_previous), abs(dist_to_next), na.rm = TRUE))
    aligned$dist_to_best[aligned$dist_to_best == Inf] <- NA
    
    aligned_seqs <- apply(aligned, 1, function(x) {
      if(is.na(x["dist_to_best"])) {
        return("")
      } else if(x["dist_to_best"] %>% as.numeric() %>% is_greater_than(5)) {
        return("")
      }
      start <- as.numeric(x["start"])
      end <- as.numeric(x["end"])
      str_sub(seq_i, max(0,start), max(0,end))
    }) %>% set_names(aligned$seq)
    aligned_seqs["sample_id"] <- ""
    if(aligned_seqs["SAMPLE"] != "")
      aligned_seqs["sample_id"] <- aligned[1, "hit"]
    
    return(c(kma_hits[i,c("read","gem", "ab", "sample")], aligned_seqs))
    
  }) %>% t()
  return(results)
}

##########################
###  PASS INPUT FILES  ###
##########################

### READ BARCODE DATA ###

# Read in all the barcode information
samples      <- readDNAStringSet(filepath = sample_file)
a_oligos     <- readDNAStringSet(filepath = oligo_A_file)
b_oligos     <- readDNAStringSet(filepath = oligo_B_file)
barcode_info <- readDNAStringSet(filepath = barcode_file)

# Read data from barcode_info
FWD_PRIMER <- barcode_info$A_primer %>% as.character()

TS         <- barcode_info$TemplateSwitch %>% as.character()
B_PRIMER   <- barcode_info$B_primer %>% as.character()
REV_PRIMER <- paste0(TS, B_PRIMER)

ANNEAL     <- barcode_info$Annealing %>% as.character()
A_N        <- barcode_info$A_N %>% as.character()
B_N        <- barcode_info$B_N %>% as.character()

# Make a list of all the barcode information
BARCODE <- list(
  SAMPLE = samples,
  A_PRIM = FWD_PRIMER,
  A_N6 = A_N,
  A_OLIGO = as.character(a_oligos),
  ANNEAL = ANNEAL,
  B_OLIGO = as.character(reverseComplement(b_oligos)),
  B_N6 = B_N,
  B_PRIM = rc(REV_PRIMER)
)

## Make sure barcode names match between what comes out of KMA and what is in the BARCODES object
names(BARCODE$A_OLIGO) <- str_extract(names(BARCODE$A_OLIGO), "A\\d+$")
names(BARCODE$B_OLIGO) <- str_extract(names(BARCODE$B_OLIGO), "B\\d+$")

### READ KMA DATA ###

### Read in the KMA data ###
kma <- read_tsv(kma_file, col_names = c("seq", "v1", "v2", "v3", "v4", "ab", "read", "brc", "label"))
kma %<>% tidyr::separate(ab, into = c("ab", "sample"), sep = "_", extra = "merge")
kma %<>% mutate(A = str_extract(ab, "A\\d+"), B = str_extract(ab, "B\\d+"))

# number of samples ids(akeys) from KMA. Set to 1 when you analyse without the sample ids
max_hits_thresholds <- 5  

# Clean and modify the KMA dataframe   
kma_hits_ALL <- as.data.frame(kma[kma$v1 <= max_hits_thresholds, ])
kma_hits_ALL$gem <- str_extract(kma_hits_ALL$read, "[ACTG]{16}")
kma_hits_ALL$read <- str_replace(kma_hits_ALL$read, " .*", "")
kma_hits_ALL$sample <- str_replace(kma_hits_ALL$sample, " .*", "")

print('Number of rows in KMA data')
print(nrow(kma_hits_ALL))

#####################
### MAIN FUNCTION ###
#####################

ncores <- 20 # Number of cores to be used
nlines <- 500 # define subchunks of N number of lines to minimize the RAM useage 

# Check if number of cores is too large!
# max_ncores is the number of cores avaliable on the current machine
# minus 3 (always leave 3 cores unused)
max_ncores <- detectCores()-3 
if (ncores > max_ncores){
  ncores <- max_ncores
}

# split data into the number of avaliable cores
ind <- rep(1:ncores, length.out = nrow(kma_hits_ALL)) %>% sort
kma_hits_parallel <- split(kma_hits_ALL, ind)

# Use mclapply to parallelized the process to speed up the calculations
parallel_results <- mclapply(X = kma_hits_parallel, FUN = function(chunk) {
  
  # Split up the data chunks into smaller chunk lists with only 
  # N number og lines defined by nlines 
  
  n <- nrow(chunk)
  chunk_ind <- sort(rep(c(1:ceiling(n/nlines)), nlines))[1:n]
  chunk_list <- split(chunk, chunk_ind)
  
  #chunk_results <- lapply(chunk_list[1:2], function(kma_hits) {
  
  # Each sequence in chunk list is alinged to the each element 
  # in the barcode design (eg. sample, primers, UMIs, and oligos) 
  
  chunk_results <- lapply(chunk_list, fx) %>% do.call(rbind, .)
  
  return(chunk_results)
}, mc.cores = ncores)

# combine final results
all_results <- do.call(rbind, parallel_results)
write.table(all_results, out_file, row.names = FALSE, quote = FALSE, sep = "\t")

print('DONE')





