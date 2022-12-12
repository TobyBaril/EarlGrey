library(tidyverse)

in_seq <- "seq/Chiloscyllium_punctatum-families.fa"
data_dir <- "TS_Chiloscyllium_punctatum-families.fa_5438/"
seq_name <- sub(".*\\/", "", in_seq)

og_seq<-Biostrings::readDNAStringSet(in_seq)
og_tibble <- tibble(seqnames = names(og_seq), width = width(og_seq)) %>%
  mutate(raw_names = sub("#.*", "", seqnames))

chimeric_seq<-Biostrings::readDNAStringSet(paste0(data_dir, "/chimeras/chimeric_", seq_name))
unacceptable_seq <-Biostrings::readDNAStringSet(paste0(data_dir, "/chimeras/questionable_", seq_name))
chimeric_seq <- c(chimeric_seq, unacceptable_seq)
chimeric_tibble <- tibble(seqnames = names(chimeric_seq), width = width(chimeric_seq)) %>%
  mutate(raw_names = sub("#.*", "", seqnames))

rps_blast_out <- read_tsv(file = paste0(data_dir, "chimeras/", seq_name, ".rps.out"),
                          col_names = c("seqnames", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                        "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle"),
                          show_col_types = FALSE) %>%
  tidyr::separate(stitle, into = c("ref", "abbrev", "full"), sep = ", ", extra = "drop") %>%
  dplyr::filter(length >= 0.5*slen)

acceptable_domains <- rbind(read_tsv("data/acceptable_domains.tsv"), read_tsv("data/additional_domains.tsv", col_select = c(ref, abbrev)))
unacceptable_domains <- rbind(read_tsv("data/unacceptable_domains.tsv", col_select = c(ref, abbrev)))

chimeric_tibble <- rps_blast_out[sub("#.*", "", rps_blast_out$seqnames) %in% c(chimeric_tibble$raw_names),]

unacceptable <- chimeric_tibble[!chimeric_tibble$ref %in% acceptable_domains$ref,]$seqnames

chimeric_tibble <- rps_blast_out[sub("#.*", "", rps_blast_out$seqnames) %in% sub("#.*", "", unacceptable),]

chimeric_tibble <- chimeric_tibble %>%
  mutate(missing = ifelse(ref %in% acceptable_domains$ref, FALSE, TRUE))
chimeric_names <- base::unique(chimeric_tibble$seqnames)
chimeric_names

i=4

to_plot <- chimeric_tibble[chimeric_tibble$seqnames == chimeric_names[i],] %>% arrange(qstart)
to_plot$y = 1:nrow(to_plot)
to_plot$text_x = (to_plot$qstart + to_plot$qend)/2
ggplot(to_plot) +
  geom_rect(aes(xmin = qstart, xmax = qend, ymin = y-1, ymax = y, fill = missing)) +
  geom_text(aes(x = text_x, y = y-0.5, label = abbrev)) +
  scale_x_continuous(limits = c(0,to_plot$qlen[1])) +
  ggtitle(to_plot$seqnames[1])
to_plot %>%
  filter(!ref %in% unacceptable_domains$ref, !ref %in% acceptable_domains$ref) %>%
  dplyr::select(ref, abbrev) %>%
  base::unique() %>%
  mutate(organism = "Chiloscyllium_punctatum") %>%
  write_tsv("new_unacceptable_domains.tsv")
