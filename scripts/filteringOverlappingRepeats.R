# load libraries

#library.path <- .libPaths()
#library(GenomicFeatures, lib.loc = library.path)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(ape)

# set options

options(scipen = 100, stringsAsFactors = FALSE)

#get inputs

args <- commandArgs()
print(args)
gff.in <- args[6]
gff.out <- args[7]

# read gff

input <- read.gff(gff.in)
#input <- read.gff("/media/toby/projectDrive/earlGrey/analysis/simulatedData/simulatingSequencesRoundTwo/benchmarkingAnalysis/overlappingAnnotations/simulatedDataset.400mb_out_sequence_nest.fasta.mod.EDTA.TEanno.gff3.sorted")

# add ID column

input$id <- paste(input[, 3], input[, 9], sep = "---")

# convert to GRanges object

gr <- with(input, GRanges(
  seqnames = input[, 1],
  IRanges(start = input[, 4], end = input[, 5]),
  id = input[, 10]))

# find overlaps

hits <- findOverlaps(gr, gr, minoverlap = 1)

# remove self hits

hits <- hits[queryHits(hits) != subjectHits(hits)]

# Determine features that are shorter than the overlapping feature

mcols(hits)$queryWidth = width(gr[queryHits(hits)]);
mcols(hits)$subjectWidth = width(gr[subjectHits(hits)]);
mcols(hits)$hit <- ifelse(
  mcols(hits)$queryWidth < mcols(hits)$subjectWidth, 
  queryHits(hits), 
  subjectHits(hits));

# Remove those shorter overlapping features
gr.final <- gr[-unique(mcols(hits)$hit)]
output <- input[input$id %in% gr.final$id,1:9]

# Write Table

write.table(output, gff.out, sep = "\t", quote = F, row.names = F, col.names = F)


