library(tidyverse)
library(plyr)
library(dplyr)
library(magrittr)
library(data.table)
library(doMC)
library(GenomicRanges)
options(scipen = 100, stringsAsFactors = FALSE)

  #####

args <- commandArgs()
print(args)
file <- args[6]
output <- args[7]
genomeSize <- args[8]
output2 <- args[9]
filtBed <- args[10]
filtSum <- args[11]
            
# Step 1 - Read in rmerge table

repeats <- fread(file)
colnames(repeats) <- c("Scaffold", "ReMa", "Repeat", "Start", "End", "Score", "Strand", "Dot", "Desc")

# Step 2 - Calculate size of repeats
  
repeats$size <- abs((repeats$End - repeats$Start)) + 1
repeats$Start <- as.double(repeats$Start) - 1
repeats$End <- as.double(repeats$End)

# Step 3 - Rectify overlaps

# Make GenomicRanges Matrix
#gr <- makeGRangesFromDataFrame(repeats,
#                               keep.extra.columns = TRUE,
#                               seqinfo = NULL,
#                               seqnames.field = "Scaffold",
#                               start.field = as.character("Start"),
#                               end.field = as.character("End"),
#                               strand.field = "Strand")

# Find overlaps
#hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
# Remove self-overlaps as query and subject are the same
#hits <- hits[queryHits(hits) != subjectHits(hits)]

# Determine features that are shorter than the overlapping feature
#mcols(hits)$querysize = size(gr[queryHits(hits)])
#mcols(hits)$subjectsize = size(gr[subjectHits(hits)])
#mcols(hits)$hit <- ifelse(
#  mcols(hits)$querysize < mcols(hits)$subjectsize,
#  queryHits(hits),
#  subjectHits(hits))

# Remove shorter overlapping features
#gr.final <- gr[-unique(mcols(hits)$hit)]

# Check there are no overlaps
#test <- findOverlaps(gr.final, gr.final, ignore.strand = TRUE)
# Remove self-overlaps and there should be none left
#test <- test[queryHits(test) != subjectHits(test)]
# There are now no overlaps
#rm(test)
      
# Turn back into dataframe
#mergedRepeats <- as.data.frame(gr.final)
#mergedRepeatsOut <- mergedRepeats[,c(1:3,7,8,5)]
#write.table(mergedRepeatsOut, file = output, quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)

mergedRepeats <- repeats
mergedRepeatsOut <- repeats[,c(1,4,5,3,6,7)]
write.table(mergedRepeatsOut, file = output, quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)
  
# Generate RM style summary table for large Repeats

sine <- mergedRepeats %>% filter(str_detect(Repeat, "SINE")) %>% distinct()
line <- mergedRepeats %>% filter(str_detect(Repeat, "LINE")) %>% distinct()
ltr <- mergedRepeats %>% filter(str_detect(Repeat, "LTR")) %>% distinct()
dna <- mergedRepeats %>% filter(str_detect(Repeat, "DNA")) %>% distinct()
rc <- mergedRepeats %>% filter(str_detect(Repeat, "RC")) %>% distinct()
unknown <- mergedRepeats %>% filter(str_detect(Repeat, "Unknown")) %>% distinct()
other <- mergedRepeats %>% filter(str_detect(Repeat, "ARTEFACT|Low|Retrop|Satellite|Simple|RNA")) %>% distinct()
other <- other[! other$Repeat == "SINE",]

# Total bp coverage of each class
sineCov <- sum(as.numeric(as.character(sine$size)))
lineCov <- sum(as.numeric(as.character(line$size)))
ltrCov <- sum(as.numeric(as.character(ltr$size)))
dnaCov <- sum(as.numeric(as.character(dna$size)))
rcCov <- sum(as.numeric(as.character(rc$size)))
unknownCov <- sum(as.numeric(as.character(unknown$size)))
otherCov <- sum(as.numeric(as.character(other$size)))
totalCov <- (sineCov + lineCov + ltrCov + dnaCov + rcCov + unknownCov) 
unmaskCov <- (as.numeric(genomeSize) - as.numeric(totalCov))

# Number of elements
sineTal <- as.numeric(tally(sine))
lineTal <- as.numeric(tally(line))
ltrTal <- as.numeric(tally(ltr))
dnaTal <- as.numeric(tally(dna))
rcTal <- as.numeric(tally(rc))
unknownTal <- as.numeric(tally(unknown))
otherTal <- as.numeric(tally(other))
totalTal <- (sineTal + lineTal + ltrTal + dnaTal + rcTal + unknownTal)

# Calculate seq percentage
genomeSize <- as.numeric(genomeSize)
sinePer <- (sineCov / genomeSize) * 100
linePer <- (lineCov / genomeSize) * 100
ltrPer <- (ltrCov / genomeSize) * 100
dnaPer <- (dnaCov / genomeSize) * 100
rcPer <- (rcCov / genomeSize) * 100
unknownPer <- (unknownCov / genomeSize) * 100
otherPer <- (otherCov / genomeSize) * 100
totalPer <- (sinePer + linePer + ltrPer + dnaPer + rcPer + unknownPer)
unmaskPer <- (100 - totalPer)

# Make Table  
eleTyp <- c("SINE", "LINE", "LTR", "DNA", "RC", "Unknown", "Total Interspersed Repeat", "Other", "Unmasked")
numEle <- c(sineTal, lineTal, ltrTal, dnaTal, rcTal, unknownTal, totalTal, otherTal, as.character("notApplicable"))
lenOcc <- c(sineCov, lineCov, ltrCov, dnaCov, rcCov, unknownCov, totalCov, otherCov, unmaskCov)
perSeq <- c(sinePer, linePer, ltrPer, dnaPer, rcPer, unknownPer, totalPer, otherPer, unmaskPer)
size <- c(rep(genomeSize, 9))

summary <- data.frame(eleTyp, numEle, lenOcc, perSeq, size)
colnames(summary) <- c("Family", "Number of Elements", "Length Occupied", "Percentage of Sequence", "GenomeSize")

    # Print table
write.table(summary, output2, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

filteredRepeatsOut <- mergedRepeatsOut
filteredRepeatsOut$length <- abs(filteredRepeatsOut$End - filteredRepeatsOut$Start)
filteredRepeatsOut <- filteredRepeatsOut[filteredRepeatsOut$length > 100,]
filteredRepeatsOut2 <- filteredRepeatsOut[,1:6]
write.table(filteredRepeatsOut2, file = filtBed, quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)


# Generate RM style summary table for large Repeats

sine <- filteredRepeatsOut %>% filter(str_detect(Repeat, "SINE")) %>% distinct()
line <- filteredRepeatsOut %>% filter(str_detect(Repeat, "LINE")) %>% distinct()
ltr <- filteredRepeatsOut %>% filter(str_detect(Repeat, "LTR")) %>% distinct()
dna <- filteredRepeatsOut %>% filter(str_detect(Repeat, "DNA")) %>% distinct()
rc <- filteredRepeatsOut %>% filter(str_detect(Repeat, "RC")) %>% distinct()
unknown <- filteredRepeatsOut %>% filter(str_detect(Repeat, "Unknown")) %>% distinct()
other <- filteredRepeatsOut %>% filter(str_detect(Repeat, "ARTEFACT|Low|Retrop|Satellite|Simple|RNA")) %>% distinct()
other <- other[! other$Repeat == "SINE",]

# Total bp coverage of each class
sineCov <- sum(as.numeric(as.character(sine$length)))
lineCov <- sum(as.numeric(as.character(line$length)))
ltrCov <- sum(as.numeric(as.character(ltr$length)))
dnaCov <- sum(as.numeric(as.character(dna$length)))
rcCov <- sum(as.numeric(as.character(rc$length)))
unknownCov <- sum(as.numeric(as.character(unknown$length)))
otherCov <- sum(as.numeric(as.character(other$length)))
totalCov <- (sineCov + lineCov + ltrCov + dnaCov + unknownCov) 
unmaskCov <- (as.numeric(genomeSize) - as.numeric(totalCov))

# Number of elements
sineTal <- as.numeric(tally(sine))
lineTal <- as.numeric(tally(line))
ltrTal <- as.numeric(tally(ltr))
dnaTal <- as.numeric(tally(dna))
rcTal <- as.numeric(tally(rc))
unknownTal <- as.numeric(tally(unknown))
otherTal <- as.numeric(tally(other))
totalTal <- (sineTal + lineTal + ltrTal + dnaTal + rcTal + unknownTal)

# Calculate seq percentage
genomeSize <- as.numeric(genomeSize)
sinePer <- (sineCov / genomeSize) * 100
linePer <- (lineCov / genomeSize) * 100
ltrPer <- (ltrCov / genomeSize) * 100
dnaPer <- (dnaCov / genomeSize) * 100
rcPer <- (rcCov / genomeSize) * 100
unknownPer <- (unknownCov / genomeSize) * 100
otherPer <- (otherCov / genomeSize) * 100
totalPer <- (sinePer + linePer + ltrPer + dnaPer + rcPer + unknownPer)
unmaskPer <- (100 - totalPer)

# Make Table  
eleTyp <- c("SINE", "LINE", "LTR", "DNA", "RC", "Unknown", "Total Interspersed Repeat", "Other", "Unmasked")
numEle <- c(sineTal, lineTal, ltrTal, dnaTal, rcTal, unknownTal, totalTal, otherTal, as.character("notApplicable"))
lenOcc <- c(sineCov, lineCov, ltrCov, dnaCov, rcCov, unknownCov, totalCov, otherCov, unmaskCov)
perSeq <- c(sinePer, linePer, ltrPer, dnaPer, rcPer, unknownPer, totalPer, otherPer, unmaskPer)
size <- c(rep(genomeSize, 9))

summary2 <- data.frame(eleTyp, numEle, lenOcc, perSeq, size)
colnames(summary2) <- c("Family", "Number of Elements", "Length Occupied", "Percentage of Sequence", "GenomeSize")

# Print table
write.table(summary2, file = filtSum, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
