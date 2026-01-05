# load libraries
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(tidyverse))

# set options

options(scipen = 100, stringsAsFactors = FALSE)

#get inputs

args <- commandArgs()
# print(args)
gff.in <- args[6]
gff.out <- args[7]

# read gff

input <- read.gff(gff.in)

# cut overlapping regions in half
input %<>%
  arrange(seqid, start) %>%
  mutate(new.start = case_when(seqid == lag(seqid) & start < lag(end) & end > lag(end) ~ as.integer((start + ((lag(end) - start)/2)) + 1),
                               seqid == lag(seqid) & start == lag(end) ~ as.integer(start + 1),
                               .default = start),
         new.end = case_when(seqid == lead(seqid) & end > lead(start) & end < lead(end) ~ as.integer((end - (end - lead(start))/2)),
                             seqid == lead(seqid) & end == lead(start) ~ as.integer(end),
                             .default = end)) %>%
  mutate(start = new.start,
         end = new.end) %>%
  select(-c(new.start, new.end))

# identify nested and unnested TEs

# nest detect function
detect_nesting <- function(df) {
  df %>%
    arrange(seqid, start) %>%
    mutate(
      nested = case_when(
        seqid == lag(seqid) &
          start > lag(start) &
          end   < lag(end) ~ "FULLY_NESTED",
        TRUE ~ "NOT_NESTED"
      )
    )
}

# iterate to detect all nestings
all_nested <- list()

iteration <- 1
current <- input

repeat {
  print(paste0("Searching for Nested TEs. Iteration:", iteration))
  annotated <- detect_nesting(current)
  
  newly_nested <- annotated %>%
    filter(nested == "FULLY_NESTED") %>%
    mutate(nesting_round = iteration)
  
  if (nrow(newly_nested) == 0) break
  
  all_nested[[iteration]] <- newly_nested
  
  current <- annotated %>%
    filter(nested == "NOT_NESTED") %>%
    select(-nested)
  
  iteration <- iteration + 1
}

nested_only <- bind_rows(all_nested)
unnested_only <- current

# recombine to make GFF
nested_only <- nested_only %>%
  mutate(attributes = paste0(attributes, ";nested=", nested, "-round", nesting_round)) %>%
  select(-c(nested, nesting_round))

output.gff <- bind_rows(nested_only, unnested_only) %>%
  arrange(seqid, start)

# Write Table
write.table(output.gff, gff.out, sep = "\t", quote = F, row.names = F, col.names = F)
