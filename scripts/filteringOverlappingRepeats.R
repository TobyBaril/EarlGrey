# load libraries

#library.path <- .libPaths()
#library(GenomicFeatures, lib.loc = library.path)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(ape)
library(tidyverse)

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

# sort table by scaffold and start

input <- input %>% 
  arrange(seqid, start)

# cut overlapping regions in half

for (i in 2:length(input$seqid)) {
  if (input$seqid[i] == input$seqid[i-1]) {
    if (input$start[i] < input$end[i-1]) {
      ovr <- (input$end[i-1] - input$start[i]) / 2
      input$start[i] <- as.integer(input$start[i] + ovr)
      input$end[i-1] <- as.integer(input$end[i-1] - ovr)
      if (input$start[i] == input$end[i-1]) {
        input$start[i] <- as.integer(input$start[i] + 1)
      }
    }
  }
}

# Write Table

write.table(input, gff.out, sep = "\t", quote = F, row.names = F, col.names = F)


