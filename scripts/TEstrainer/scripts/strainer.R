#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("-i", "--in_seq"), default=NA, type = "character", help="Path to input sequence (required)"),
  make_option(c("-p", "--plot_chimeras"), type = "logical", help="Set to TRUE to make plots of chimeras", default = TRUE),
  make_option(c("-d", "--directory"), default=NA, type = "character", help="Path to directory (required)",)
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
suppressPackageStartupMessages(library(grid))

# read in fasta
rm_seq_in <- readDNAStringSet(paste0(opt$in_seq))
names(rm_seq_in) <- sub(" .*", "", names(rm_seq_in))
rm_seq_info <- tibble(seqnames = names(rm_seq_in), width = width(rm_seq_in))

compiled_acceptable <- tibble()
truly_chimeric_ranges <- GRanges()
questionable <- tibble()
no_domains_seq <- DNAStringSet()

# read in acceptable domains and exceptional domains and manually added additional domains
exceptional_domains <- read_tsv("data/exceptional_domains.tsv", show_col_types = FALSE) %>%
  dplyr::select(ref, class) %>% dplyr::rename(excep_class = class)
additional_domains <- read_tsv("data/additional_domains.tsv", show_col_types = FALSE) %>%
  dplyr::select(ref, class) %>% dplyr::rename(excep_class = class)
exceptional_domains <- rbind(exceptional_domains, additional_domains)

acceptable_domains <- read_tsv("data/acceptable_domains.tsv", show_col_types = FALSE)

if(file.size(opt$rps_table)==0){
  writeXStringSet(rm_seq_in, paste0(opt$directory, "/chimeras/clean_", opt$out_seq))
  quit()
}

# read rps blast out
rps_blast_out <- read_tsv(file = opt$rps_table,
                                           col_names = c("seqnames", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                                         "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle"),
                                           show_col_types = FALSE) %>%
                   tidyr::separate(stitle, into = c("ref", "abbrev", "full"), sep = ", ", extra = "drop") %>%
  dplyr::filter(length >= 0.5*slen)

rps_blast_out$n <- 1:nrow(rps_blast_out)

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

# check if domains fit exceptional circumstances (grepl used to allow for the various ERVs)
excep_check <- chimeric %>%
  inner_join(exceptional_domains)

for( i in 1:nrow(excep_check)){
  if(grepl(pattern = excep_check$excep_class[i], x = excep_check$seqnames[i])){
    excep_check$unacceptable[i] <- "false"
  } else {
    excep_check$unacceptable[i] <- "true"
  }
}

# select single for each excep check, ensuring acceptable retained if exists
excep_check <- excep_check %>%
  group_by(n) %>%
  dplyr::arrange(unacceptable) %>%
  dplyr::slice(1) %>%
  dplyr::select(-excep_class)

# add excep_checked back into fold
chimeric <- chimeric %>%
  dplyr::filter(!n %in% excep_check$n) %>%
  base::rbind(excep_check) %>%
  arrange(n)

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
# plotting section for manual consideration of truly_chimeric elements
if(opt$plot == TRUE & nrow(truly_chimeric) > 0){
  
  if(!file.exists(paste0(opt$directory, "/chimeras/plots/"))){
    dir.create(paste0(opt$directory, "/chimeras/plots/"))
  }
  
  for(i in 1:length(base::unique(truly_chimeric$seqnames))){
    
    to_plot <- truly_chimeric[truly_chimeric$seqnames == base::unique(truly_chimeric$seqnames)[i],] %>% arrange(qstart) %>%
      mutate(missing = ifelse(ref %in% acceptable_domains$ref, FALSE, TRUE))
    to_plot$y = 1:nrow(to_plot)
    to_plot$text_x = (to_plot$qstart + to_plot$qend)/2
    ggplot(to_plot) +
      geom_rect(aes(xmin = qstart, xmax = qend, ymin = y-1, ymax = y, fill = missing)) +
      geom_segment(aes(x = qstart, xend = qend, y = y, yend = y),
                   arrow = arrow(length = unit(10,"points"))) +
      geom_text(aes(x = text_x, y = y-0.5, label = abbrev)) +
      scale_x_continuous(limits = c(0,to_plot$qlen[1]), name = "Coordinates (bp)", expand = c(0,0)) +
      scale_y_continuous(name = NULL, breaks = NULL) +
      scale_fill_discrete(name="Atypical of TEs?") +
      ggtitle(to_plot$seqnames[1]) +
      theme_bw()
    
    ggsave(filename = paste0(opt$directory, "/chimeras/plots/", gsub("/", "_", base::unique(truly_chimeric$seqnames)[i]), ".svg"), device = "svg")
    
  }
  
}
