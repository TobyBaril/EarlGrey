#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("-i", "--in_seq"), default=NA, type = "character", help="Path to input sequence (required)"),
  make_option(c("-d", "--directory"), default=NA, type = "character", help="Path to directory (required)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# check variables provided
if(is.na(opt$in_seq)){
  stop("Path to in sequence must be supplied")
}
if(is.na(opt$directory)){
  stop("Path to directory must be supplied")
}

opt$out_seq <- sub(".*/", "", opt$in_seq)
opt$rps_table <- paste0(opt$directory, "/chimeras/", opt$out_seq, ".rps.out")

# make empty variable function
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(BSgenome))

# read in fasta
rm_seq_in <- readDNAStringSet(paste0(opt$in_seq))
names(rm_seq_in) <- sub(" .*", "", names(rm_seq_in))
rm_seq_info <- tibble(seqnames = names(rm_seq_in), width = width(rm_seq_in))

compiled_acceptable <- tibble()
truly_chimeric_ranges <- GRanges()
questionable <- tibble()
no_domains_seq <- DNAStringSet()

# read in acceptable domains, rbind additionally discovered
additional_domains <- read_tsv("data/additional_domains.tsv", show_col_types = FALSE) %>%
  dplyr::select(ref, abbrev)
acceptable_domains <- read_tsv("data/acceptable_domains.tsv", show_col_types = FALSE) %>%
  rbind(additional_domains)
unacceptable_domains <- read_tsv("data/unacceptable_domains.tsv", show_col_types = FALSE)

if(file.size(opt$rps_table)==0){
  writeXStringSet(rm_seq_in, paste0(opt$directory, "/chimeras/clean_", opt$out_seq))
  quit()
}

# read rps blast out
rps_blast_out <- read_tsv(file = opt$rps_table,
                                           col_names = c("seqnames", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                                         "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle"),
                                           show_col_types = FALSE) %>%
                   tidyr::separate(stitle, into = c("ref", "abbrev", "full"), sep = ", ", extra = "drop")

rps_blast_out <- rps_blast_out %>% dplyr::filter(length >= 0.5*slen)

# select acceptable and non acceptable
contains_acceptable <- rps_blast_out[rps_blast_out$ref %in% acceptable_domains$ref,]
contains_unacceptable <- rps_blast_out[!rps_blast_out$ref %in% acceptable_domains$ref,]

# Filter repeats containing only suitable domains
completely_acceptable <- rps_blast_out %>%
  filter(seqnames %in% contains_acceptable$seqnames,
         !seqnames %in% contains_unacceptable$seqnames) %>%
  arrange(seqnames, qstart)

# Filter repeats containing no suitable domains
questionable <- rps_blast_out %>%
  filter(!seqnames %in% contains_acceptable$seqnames,
         seqnames %in% contains_unacceptable$seqnames) %>%
  arrange(seqnames, qstart)

# Filter repeats containing both suitable and unsuitable domains
chimeric <- rps_blast_out %>%
  filter(seqnames %in% contains_acceptable$seqnames,
         seqnames %in% contains_unacceptable$seqnames) %>%
  arrange(seqnames, qstart) %>%
  mutate(unacceptable = ifelse(ref %in% acceptable_domains$ref, "false", "true"))

# Determine regions with acceptable domains overlapping "unacceptable"
chimeric_ranges <- chimeric %>%
  dplyr::mutate(start = ifelse(qstart < qend, qstart, qend),
                end = ifelse(qstart > qend, qstart, qend),
                strand = ifelse(qstart < qend, "+", "-")) %>%
  dplyr::select(seqnames, start, end, strand, abbrev, unacceptable) %>%
  plyranges::as_granges()

chimeric_ranges_unacceptable_false <- chimeric_ranges[chimeric_ranges$unacceptable == "false",]
chimeric_ranges_unacceptable_true <- chimeric_ranges[chimeric_ranges$unacceptable == "true",]

# find regions which don't overlap with acceptable domains
truly_chimeric_ranges <- filter_by_non_overlaps(chimeric_ranges_unacceptable_true, chimeric_ranges_unacceptable_false)

truly_chimeric <- chimeric[chimeric$seqnames %in% seqnames(truly_chimeric_ranges),]
false_positive_chimeric <- chimeric[!chimeric$seqnames %in% seqnames(truly_chimeric_ranges),] %>%
  dplyr::select(-unacceptable)

# remove repeats unacceptable domains
unacceptable_chimeric <- truly_chimeric %>%
  filter(ref %in% unacceptable_domains$ref)
unacceptable_chimeric <- truly_chimeric %>%
  filter(seqnames %in% unacceptable_chimeric$seqnames) %>%
  dplyr::select(-unacceptable)
questionable <- rbind(questionable, unacceptable_chimeric)

truly_chimeric <- truly_chimeric %>%
  filter(!seqnames %in% unacceptable_chimeric$seqnames)

# flip if domain with highest bitscore is in reverse strand
compiled_acceptable <- rbind(completely_acceptable, false_positive_chimeric) %>%
  dplyr::mutate(strand = ifelse(qstart < qend, "+", "-"), start = 1, end = qlen) %>%
  dplyr::group_by(seqnames) %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  plyranges::as_granges()

# identify sequences with no domains
no_domains_seq <- rm_seq_in[!names(rm_seq_in) %in% c(chimeric$seqnames, completely_acceptable$seqnames, questionable$seqnames)]

## add step to combine data
# write to file (check if any filtered, if not write all in to output)
completely_acceptable_seq <- getSeq(rm_seq_in, compiled_acceptable)
names(completely_acceptable_seq) <- seqnames(compiled_acceptable)
writeXStringSet(c(completely_acceptable_seq, no_domains_seq), paste0(opt$directory, "/chimeras/clean_", opt$out_seq))

# in writing to file flip if appropriate
if(nrow(truly_chimeric) > 0){
  
  truly_chimeric_ranges <- truly_chimeric %>%
    dplyr::mutate(strand = ifelse(qstart < qend, "+", "-"), start = 1, end = qlen) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::arrange(-bitscore) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    plyranges::as_granges()
  
  chimeric_seq <- getSeq(rm_seq_in, truly_chimeric_ranges)
  names(chimeric_seq) <- seqnames(truly_chimeric_ranges)
  writeXStringSet(chimeric_seq, paste0(opt$directory, "/chimeras/chimeric_", opt$out_seq))

}

if(nrow(questionable) > 0){
  
  questionable_ranges <- questionable %>%
    dplyr::mutate(strand = ifelse(qstart < qend, "+", "-"), start = 1, end = qlen) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::arrange(-bitscore) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    plyranges::as_granges()
  
  questionable_seq <- getSeq(rm_seq_in, questionable_ranges)
  names(questionable_seq) <- seqnames(questionable_ranges)
  writeXStringSet(questionable_seq, paste0(opt$directory, "/chimeras/questionable_", opt$out_seq))
  
}
