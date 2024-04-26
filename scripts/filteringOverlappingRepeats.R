# load libraries
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(tidyverse))

# set options

options(scipen = 100, stringsAsFactors = FALSE)

#get inputs

args <- commandArgs()
print(args)
gff.in <- args[6]
gff.out <- args[7]

# read gff

input <- read.gff(gff.in)
#input <- read.gff("/home/toby/projects/earlGreyREwrite/TESTENV/autoTest_RM_EarlGrey/autoTest_RM_mergedRepeats/looseMerge/autoTest_RM.rmerge.gff.sorted")
#input <- read.gff("/home/toby/projects/earlGreyREwrite/TESTENV/autoTest_newOverlapFilter_EarlGrey/autoTest_newOverlapFilter_mergedRepeats/looseMerge/autoTest_newOverlapFilter.rmerge.gff.sorted")

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

# Write Table
write.table(input, gff.out, sep = "\t", quote = F, row.names = F, col.names = F)


