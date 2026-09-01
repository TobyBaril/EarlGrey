#!/usr/bin/Rscript

library(ORFik)
library(tidyverse)
library(plyranges)
library(BSgenome)

# read in sequences
repbase_dna <- Biostrings::readDNAStringSet("~/Databases/Repbase/RepBase5May2021.fasta")
dfam_dna <- Biostrings::readDNAStringSet("seq/Dfam/Dfam_3.6_curatedonly.fasta")

# get info about sequences
repbase_seq_info <- tibble(seqnames = names(repbase_dna), width = width(repbase_dna))
repbase_seq_info <- repbase_seq_info %>%
  tidyr::separate(col = seqnames, into = c("seqnames", "class", "species"), sep = "\t") %>%
  dplyr::select(-species) %>%
  filter(!class %in% c("Unknown", "Satellite","Simple Repeat"))

# convert dfam classifications to Repbase
dfam_seq_info <- tibble(seqnames = names(dfam_dna), width = width(dfam_dna)) %>%
  tidyr::separate(col = seqnames, into = c("seqnames", "class"), sep = "#") %>%
  mutate(class = ifelse(grepl("ERV", class), "ERV", class)) %>%
  dplyr::mutate(class = case_when(
    grepl("SINE", class) ~ "SINE",
    grepl("/", class) ~ sub(".*/", "", class),
    TRUE ~ class)
  ) %>%
  mutate(
    class = ifelse(class %in% c("TcMar", "TcMar-ISRm11", "TcMar-m44", "TcMar-Pogo", "TcMar-Fot1", "TcMar-Tc2", "TcMar-Mariner", "TcMar-Tc1", "TcMar-Tc4", "TcMar-Tigger"), "Mariner/Tc1", class),
    class = ifelse(class %in% c("hAT-hAT5", "hAT-hAT19", "hAT-hATx", "hAT-hATm", "hAT-Tag1", "hAT-Ac", "hAT-Blackjack", "hAT-Tip100", "hAT-Charlie"), "hAT", class),
    class = ifelse(class %in% c("Pao"), "BEL", class),
    class = ifelse(class %in% c("CMC-Transib"), "Transib", class),
    class = ifelse(class %in% c("PiggyBac"), "piggyBac", class),
    class = ifelse(class %in% c("Naiad", "Chlamys"), "Penelope", class),
    class = ifelse(class %in% c("CMC-Chapaev", "CMC-Chapaev-3", "CMC-EnSpm"), "EnSpm/CACTA", class),
    class = ifelse(class %in% c("Academ-1"), "Academ", class),
    class = ifelse(class %in% c("Ngaro"), "DIRS", class),
    class = ifelse(class %in% c("PIF-ISL2EU"), "ISL2EU", class),
    class = ifelse(class %in% c("R2-Hero"), "Hero", class),
    class = ifelse(class %in% c("CR1-Zenon"), "CR1", class),
    class = ifelse(class %in% c("endogenous retrovirus", "Endogenous Retrovirus", "ERV-1_PM-I", "ERV-2_PM-I"), "ERV", class),
    class = ifelse(class %in% c("LTR retrotransposon", "LTR Retrotransposon"), "LTR", class),
    class = ifelse(class %in% c("RTE-RTE", "RTE-BovB"), "RTE", class)
  ) %>%
  filter(!seqnames %in% repbase_seq_info$seqnames)

dfam_dna <- dfam_dna[sub("#.*", "", names(dfam_dna)) %in% dfam_seq_info$seqnames,]
repbase_dna <- c(repbase_dna, dfam_dna)
repbase_seq_info_class_n <- as_tibble(as.data.frame(table(repbase_seq_info$class))) %>%
  dplyr::filter(Freq > 5) %>%
  dplyr::mutate(Var1 = as.character(Var1)) %>%
  dplyr::rename(class = Var1)
repbase_seq_info <- rbind(repbase_seq_info, dfam_seq_info) %>%
  filter(!class %in% c("Unknown", "Satellite", "subtelomeric", "Y-chromosome", "centromeric", "tRNA", "snRNA", "rRNA", "tRNA", "Multicopy gene"),
         class %in% repbase_seq_info_class_n$class)
repbase_seq_info$name = 1:nrow(repbase_seq_info)

repbase_dna <- repbase_dna[sub("#.*", "", sub("\t.*", "", names(repbase_dna))) %in% repbase_seq_info$seqnames,]
names(repbase_dna) <- sub("#.*", "", sub("\t.*", "", names(repbase_dna)))
writeXStringSet(repbase_dna, "predata/repbase_dfam_compiled.fasta")

#### RUN RPSTBLASTN AGAINST CDD ####


# read in repbase rpstblastn
repbase_rps_out <- readr::read_tsv(file = "predata/repbase_dfam_compiled.rps.out",
                                   col_names = c("seqnames", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                                 "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle"),
                                   show_col_types = F) %>%
  dplyr::filter(length >= 0.5*slen, evalue <= 0.01)

missing <- repbase_rps_out %>% filter(!seqnames %in% repbase_seq_info$seqnames)

repbase_rps_out <- tidyr::separate(data = repbase_rps_out, col = stitle, into = c("ref", "abbrev", "full"), sep = ", ")
repbase_rps_out <- dplyr::inner_join(repbase_rps_out, repbase_seq_info) %>%
  filter(class != "Multicopy gene")
domain_info <- repbase_rps_out %>% dplyr::select(ref, abbrev, full) %>% base::unique()

# get info about repbase classes
repbase_class_info <- dplyr::as_tibble(BiocGenerics::as.data.frame(table(repbase_seq_info$class)))
repbase_class_info <- dplyr::mutate(.data = repbase_class_info, Var1 = as.character(Var1))
repbase_class_info <- dplyr::rename(.data = repbase_class_info, class = Var1)
repbase_class_info_min <- repbase_class_info[repbase_class_info$Freq < 10,]
repbase_class_info_max <- repbase_class_info[repbase_class_info$Freq >= 10,]
repbase_class_info_max <- dplyr::arrange(.data = repbase_class_info_max, -Freq)

# loop over to get info
common_domains <- tibble()
for(i in 1:nrow(repbase_class_info_max)){
  holder <- repbase_rps_out[repbase_rps_out$class == repbase_class_info_max$class[i],] %>%
    group_by(seqnames, abbrev) %>%
    arrange(-bitscore) %>%
    dplyr::slice(1) %>% # ensure one representative for each domain for each repeat
    ungroup()
  if(nrow(holder >= 1)){
    if(repbase_class_info_max$Freq[i] > 299){
    holder_tbl <- as_tibble(as.data.frame(table(holder$ref))) %>%
      mutate(Var1 = as.character(Var1),
             perc = Freq/repbase_class_info_max$Freq[i]) %>%
      dplyr::rename(ref = Var1) %>%
      filter(perc >= 0.01, Freq >= 3) %>%
      inner_join(domain_info) %>%
      mutate(class = repbase_class_info_max$class[i]) %>%
      dplyr::arrange(-Freq)
    } else {
      holder_tbl <- as_tibble(as.data.frame(table(holder$ref))) %>%
        mutate(Var1 = as.character(Var1),
               perc = Freq/repbase_class_info_max$Freq[i]) %>%
        dplyr::rename(ref = Var1) %>%
        filter(Freq >= 3) %>%
        inner_join(domain_info) %>%
        mutate(class = repbase_class_info_max$class[i]) %>%
        dplyr::arrange(-Freq)
    }
    common_domains <- rbind(common_domains, holder_tbl)
  }
}

# write to file
base::unique(common_domains %>% select(ref, abbrev)) %>%
  dplyr::filter(!startsWith(tolower(abbrev), "ig"),
                !startsWith(tolower(abbrev), "znmc"),
                abbrev != "V-set") %>%
  write_tsv("data/acceptable_domains.tsv")

base::unique(common_domains %>% select(class, ref, abbrev)) %>%
  dplyr::filter(startsWith(tolower(abbrev), "ig") | startsWith(tolower(abbrev), "znmc") | abbrev == "V-set") %>%
  write_tsv("data/exceptional_domains.tsv")
