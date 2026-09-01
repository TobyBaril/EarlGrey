#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("-i", "--in_seq"), default=NA, type = "character", help="Input sequence (required)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if(is.na(opt$in_seq)){
  stop("Path to in sequence must be supplied")
}

# make empty variable function
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(BSgenome))

# read in sequence
in_seq <- readDNAStringSet(opt$in_seq)

# sort by package used, round number and family number
sorted_ranges <-
  tibble(seqnames = names(in_seq), start = 1, end = width(in_seq)) %>%
    dplyr::mutate(pkg = substr(seqnames, 1, 3),
                  rnd = as.integer(sub(".*-", "", sub("_family.*", "", seqnames))),
                  family = as.integer(sub("#.*", "", sub(".*_family-", "", seqnames)))) %>%
    dplyr::arrange(pkg, rnd, family) %>%
    plyranges::as_granges()

# Get sequence
out_seq = getSeq(in_seq, sorted_ranges)
names(out_seq) <- seqnames(sorted_ranges)

# Write to file
writeXStringSet(out_seq, opt$in_seq)