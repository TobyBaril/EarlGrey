#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("-i", "--in_seq"), default=NA, type = "character", help="Input sequence (required)"),
  make_option(c("-d", "--directory"), type="character", default=NULL, help="Path to data directory (required)", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

# check variables provided
if(is.na(opt$in_seq)){
  stop("Path to in sequence must be supplied")
}
if(is.na(opt$directory)){
  stop("Path to in sequence must be supplied")
}

opt$out_seq <- sub(".*/", "", opt$in_seq)

# make empty variable function
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(BSgenome))

# read in library and make bed like table
in_seq <- readDNAStringSet(opt$in_seq)
in_seq_tbl <- tibble(seqnames = names(in_seq), og_width = width(in_seq)) %>%
  dplyr::mutate(draft_seqnames = sub("#.*", "", seqnames))

# read in, rearrange and calculate percent tandem repeats for SA-SSR data
# check sassr empty, if so create empty tbl, else analyse data
if(file.size(paste0(opt$directory, "/trf/", opt$out_seq, ".sassr")) == 0){
  
  sassr_calc <- tibble(seqnames = character(), sassr_perc_tr = double())
  
  sassr_select <- tibble(seqnames = character(), start = double(), end = double(),
                         period = double(), count = double(), og_width = integer(),
                         ssr = character(), ssr_width = double(), package = "sassr")
  
} else {

  sassr <- read_tsv(paste0(opt$directory, "/trf/", opt$out_seq, ".sassr"), 
                    skip = 1, col_names = c("seqnames", "ssr", "count", "start"), show_col_types = F) %>%
    dplyr::mutate(ssr = ifelse(is.na(ssr), "NA", ssr),
                  period = as.double(width(ssr))) %>%
    dplyr::mutate(ssr_width = count*period, end = start + ssr_width, start = start +1) %>%
                  mutate(draft_seqnames = sub("#.*", "", seqnames)) %>%
    inner_join(in_seq_tbl, by = "seqnames") %>%
    arrange(seqnames) %>%
    filter(count > 2)
  
  sassr_calc <- as_tibble(reduce(as_granges(sassr))) %>%
    dplyr::select(-strand) %>%
    dplyr::mutate(seqnames = as.character(seqnames)) %>%
    group_by(seqnames) %>%
    mutate(total_width = sum(width)) %>%
    ungroup() %>%
    inner_join(in_seq_tbl, by = "seqnames") %>%
    dplyr::mutate(sassr_perc_tr = 100*total_width/og_width) %>%
    dplyr::select(seqnames, sassr_perc_tr) %>%
    base::unique()
  
  sassr_select <- sassr %>% select(seqnames, start, end, period, count, og_width, ssr, ssr_width) %>%
    mutate(package = "sassr")

}

# check trf empty, if so create empty tbl, else analyse data
# read in, rearrange and calculate percent tandem repeats for TRF data
if(file.size(paste0(opt$directory, "/trf/", opt$out_seq, ".trf")) == 0){
  trf_select <- tibble(seqnames = character(), start = double(), end = double(),
                         period = double(), count = double(), og_width = integer(),
                         ssr = character(), ssr_width = double(), package = "trf")
  
  trf_calc <- tibble(seqnames = character(), trf_perc_tr = double())
  
} else {
  trf <- read_tsv(paste0(opt$directory, "/trf/", opt$out_seq, ".trf"),
                  col_names = c("draft_seqnames", "start", "end", "period", "count", "ssr"), show_col_types = F) %>%
    mutate(ssr = ifelse(is.na(ssr), "NA", ssr),
           draft_seqnames = sub("@", "", sub("#.*", "", draft_seqnames))) %>%
    dplyr::mutate(ssr_width = end - start + 1) %>%
    inner_join(in_seq_tbl, by = "draft_seqnames")
  
  trf_select <- trf  %>%
    filter(count > 2) %>%
    dplyr::select(seqnames, start, end)
  
  trf_calc <- as_tibble(reduce(as_granges(trf_select))) %>%
    dplyr::select(-strand) %>%
    dplyr::mutate(seqnames = as.character(seqnames)) %>%
    group_by(seqnames) %>%
    mutate(total_width = sum(width)) %>%
    ungroup() %>%
    inner_join(in_seq_tbl, by = "seqnames") %>%
    dplyr::mutate(trf_perc_tr = 100*total_width/og_width) %>%
    dplyr::select(seqnames, trf_perc_tr) %>%
    base::unique()
  
  trf_select <- trf %>% select(seqnames, start, end, period, count, og_width, ssr, ssr_width) %>%
    mutate(package = "trf")
  
}

# read in, rearrange and calculate percent tandem repeats for MREPS data
if(file.size(paste0(opt$directory, "/trf/", opt$out_seq, ".mreps")) == 0){
  mreps_select <- tibble(seqnames = character(), start = double(), end = double(),
                         period = double(), count = double(), og_width = integer(),
                         ssr = character(), ssr_width = double(), package = "mreps")
  
  mreps_calc <- tibble(seqnames = character(), mreps_perc_tr = double())

} else {
  mreps <- read_tsv(paste0(opt$directory, "/trf/", opt$out_seq, ".mreps"),
                    col_names = c("draft_seqnames", "start", "end", "ssr_width", "period", "count", "error", "sequence"), show_col_types = F) %>%
    mutate(ssr = substr(x = sequence, start = 0, stop = period)) %>%
    inner_join(in_seq_tbl, by = "draft_seqnames")
  
  mreps_select <- mreps %>%
    filter(count > 2) %>%
    dplyr::select(seqnames, start, end)
  
  mreps_calc <- as_tibble(reduce(as_granges(mreps_select))) %>%
    dplyr::select(-strand) %>%
    dplyr::mutate(seqnames = as.character(seqnames)) %>%
    group_by(seqnames) %>%
    mutate(total_width = sum(width)) %>%
    ungroup() %>%
    inner_join(in_seq_tbl, by = "seqnames") %>%
    dplyr::mutate(mreps_perc_tr = 100*total_width/og_width) %>%
    dplyr::select(seqnames, mreps_perc_tr) %>%
    base::unique()
  
  mreps_select <- mreps %>% select(seqnames, start, end, period, count, og_width, ssr, ssr_width) %>%
    mutate(package = "mreps")
}

# Compile data
compiled_tr <- rbind(rbind(trf_select, mreps_select), sassr_select)
stats_tr <- full_join(sassr_calc, full_join(mreps_calc, trf_calc, by = "seqnames"), by = "seqnames")

# if no data write to og seq to file as nonsatellite and exit
if(nrow(compiled_tr) == 0){
  writeXStringSet(c(trimmed_seq, untouched_seq), paste0(opt$directory, "/trf/", opt$out_seq, ".nonsatellite"))
  stop()
}

# Identifying satellite/simple repeats
over50tr <- stats_tr %>%
  filter(sassr_perc_tr >= 50 | mreps_perc_tr >= 50 | trf_perc_tr >= 50) %>%
  mutate(sassr_perc_tr = ifelse(is.na(sassr_perc_tr), 0, sassr_perc_tr),
         mreps_perc_tr = ifelse(is.na(mreps_perc_tr), 0, mreps_perc_tr),
         trf_perc_tr = ifelse(is.na(trf_perc_tr), 0, trf_perc_tr))

satellites <- compiled_tr[compiled_tr$seqnames %in% over50tr$seqnames,] %>%
  dplyr::group_by(seqnames) %>%
  dplyr::arrange(-ssr_width, period) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(perc_tr = 100*ssr_width/og_width) %>%
  inner_join(over50tr, by = "seqnames") %>%
  dplyr::mutate(max_perc_tr = pmax(sassr_perc_tr, mreps_perc_tr, trf_perc_tr)) %>%
  dplyr::arrange(max_perc_tr)

# for macrosatellites trim to get single copy
macrosatellites <- satellites %>%
  filter(max_perc_tr > 90, count >= 2) %>%
  filter(period > 200) %>%
  mutate(start = start + period,
         end = start + period) %>%
  as_granges()
macrosatellites_seq <- getSeq(in_seq, macrosatellites)
names(macrosatellites_seq) <- sub("#.*", "#Satellite", seqnames(macrosatellites))

# leave micro to minisatellite repeats intact
other_satellites <- satellites %>%
  dplyr::filter(!seqnames %in% seqnames(macrosatellites)) %>%
  dplyr::mutate(start = 1, end = og_width) %>%
  as_granges()
other_satellites_seq <- getSeq(in_seq, other_satellites)
names(other_satellites_seq) <- sub("#.*", "#Satellite", seqnames(other_satellites))

# Combine satellite repeats and write to file
writeXStringSet(c(other_satellites_seq, macrosatellites_seq), paste0(opt$directory, "/trf/", opt$out_seq, ".satellites"))

# filtering for trimming
trim_3 <- compiled_tr %>%
  mutate(tr_len = period*count) %>%
  filter(count >=3, !seqnames %in% over50tr$seqnames, period < 20) %>%
  filter(og_width-end <= period) %>%
  dplyr::group_by(seqnames) %>%
  arrange(-end, period) %>%
  dplyr::slice(1)
trim_5 <- compiled_tr %>%
  mutate(tr_len = period*count) %>%
  filter(count >=3, !seqnames %in% over50tr$seqnames, period < 20) %>%
  filter(start <= period) %>%
  dplyr::group_by(seqnames) %>%
  arrange(-end, period) %>%
  dplyr::slice(1)

# Determine thos to be trimmed both ends
trim_both_5 <- trim_5[trim_5$seqnames %in% trim_3$seqnames,] %>%
  dplyr::mutate(start = ifelse(period < 4, end - 8, end - (2*period))) %>%
  dplyr::mutate(start = ifelse(start < 1, 1, start)) %>%
  dplyr::select(seqnames, start)
trim_both_3 <- trim_3[trim_3$seqnames %in% trim_5$seqnames,] %>%
  dplyr::mutate(end = ifelse(period < 4, start + 8, start + (2*period))) %>%
  dplyr::mutate(end = ifelse(end > og_width, og_width, end)) %>%
  dplyr::select(seqnames, end)
trim_both <- inner_join(trim_both_5, trim_both_3, by = "seqnames")

# determine those to be trimmed one end
trim_only_5 <- trim_5[!trim_5$seqnames %in% trim_3$seqnames,] %>%
  dplyr::mutate(start = ifelse(period < 4, end - 8, end - (2*period))) %>%
  dplyr::mutate(start = ifelse(start < 1, 1, start), end = og_width) %>%
  dplyr::select(seqnames, start, end)
trim_only_3 <- trim_3[!trim_3$seqnames %in% trim_5$seqnames,] %>%
  dplyr::mutate(end = ifelse(period < 4, start + 8, start + (2*period))) %>%
  dplyr::mutate(end = ifelse(end > og_width, og_width, end), start = 1) %>%
  dplyr::select(seqnames, start, end)

to_trim <- as_granges(rbind(trim_both, rbind(trim_only_5, trim_only_3)))
trimmed_seq <- getSeq(in_seq, to_trim)
names(trimmed_seq) <- seqnames(to_trim)

# determine untouched sequences
untouched_seq <-
  in_seq[!sub("#.*", "", names(in_seq)) %in% sub("#.*", "", c(names(trimmed_seq), names(macrosatellites_seq), names(other_satellites_seq))),]

# combine untouched and trimmed sequences and write to file
writeXStringSet(c(trimmed_seq, untouched_seq), paste0(opt$directory, "/trf/", opt$out_seq, ".nonsatellite"))
