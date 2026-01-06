suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(data.table))
options(scipen = 100, stringsAsFactors = FALSE)

#####

args <- commandArgs()
# print(args)
file <- args[6]
output <- args[7]
genomeSize <- args[8]
output2 <- args[9]
filtBed <- args[10]
filtSum <- args[11]
lowend <- args[12]

# Step 1 - Read in rmerge table

repeats <- fread(file)
colnames(repeats) <- c("Scaffold", "ReMa", "Repeat", "Start", "End", "Score", "Strand", "Dot", "Desc")

# Step 2 - Calculate size of repeats

# calculate size of all repeats
repeats$size <- abs((repeats$End - repeats$Start)) + 1
repeats$Start <- as.double(repeats$Start) - 1
repeats$End <- as.double(repeats$End)

mergedRepeats <- repeats
mergedRepeatsOut <- repeats[,c(1,4,5,3,6,7)]
write.table(mergedRepeatsOut, file = output, quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)

# Generate RM style summary table for all unnested repeats, and nested ones too

sine <- mergedRepeats %>% filter(str_detect(Repeat, "SINE")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
line <- mergedRepeats %>% filter(str_detect(Repeat, "LINE")) %>% filter(! str_detect(Repeat, "Penelope")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
ltr <- mergedRepeats %>% filter(str_detect(Repeat, "LTR")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
dna <- mergedRepeats %>% filter(str_detect(Repeat, "DNA")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
rc <- mergedRepeats %>% filter(str_detect(Repeat, "RC")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
ple <- mergedRepeats %>% filter(str_detect(Repeat, "Penelope")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
unknown <- mergedRepeats %>% filter(str_detect(Repeat, "Unknown")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
other <- mergedRepeats %>% filter(str_detect(Repeat, "ARTEFACT|Low|Retrop|Satellite|Simple|RNA")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
other <- other[! other$Repeat == "SINE",]

sine.nest <- mergedRepeats %>% filter(str_detect(Repeat, "SINE")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
line.nest <- mergedRepeats %>% filter(str_detect(Repeat, "LINE")) %>% filter(! str_detect(Repeat, "Penelope")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
ltr.nest <- mergedRepeats %>% filter(str_detect(Repeat, "LTR")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
dna.nest <- mergedRepeats %>% filter(str_detect(Repeat, "DNA")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
rc.nest <- mergedRepeats %>% filter(str_detect(Repeat, "RC")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
ple.nest <- mergedRepeats %>% filter(str_detect(Repeat, "Penelope")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
unknown.nest <- mergedRepeats %>% filter(str_detect(Repeat, "Unknown")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
other.nest <- mergedRepeats %>% filter(str_detect(Repeat, "ARTEFACT|Low|Retrop|Satellite|Simple|RNA")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
other.nest <- other.nest[! other.nest$Repeat == "SINE",]

# Total bp coverage of each class
sineCov <- sum(as.numeric(as.character(sine$size)))
lineCov <- sum(as.numeric(as.character(line$size)))
ltrCov <- sum(as.numeric(as.character(ltr$size)))
dnaCov <- sum(as.numeric(as.character(dna$size)))
rcCov <- sum(as.numeric(as.character(rc$size)))
pleCov <- sum(as.numeric(as.character(ple$size)))
unknownCov <- sum(as.numeric(as.character(unknown$size)))
otherCov <- sum(as.numeric(as.character(other$size)))
totalCov <- (sineCov + lineCov + ltrCov + dnaCov + rcCov + pleCov + unknownCov) 
unmaskCov <- (as.numeric(genomeSize) - as.numeric(totalCov))

sineCov.nest <- sum(as.numeric(as.character(sine.nest$size)))
lineCov.nest <- sum(as.numeric(as.character(line.nest$size)))
ltrCov.nest <- sum(as.numeric(as.character(ltr.nest$size)))
dnaCov.nest <- sum(as.numeric(as.character(dna.nest$size)))
rcCov.nest <- sum(as.numeric(as.character(rc.nest$size)))
pleCov.nest <- sum(as.numeric(as.character(ple.nest$size)))
unknownCov.nest <- sum(as.numeric(as.character(unknown.nest$size)))
otherCov.nest <- sum(as.numeric(as.character(other.nest$size)))

# Number of elements
sineTal <- as.numeric(tally(sine))
lineTal <- as.numeric(tally(line))
ltrTal <- as.numeric(tally(ltr))
dnaTal <- as.numeric(tally(dna))
rcTal <- as.numeric(tally(rc))
pleTal <- as.numeric(tally(ple))
unknownTal <- as.numeric(tally(unknown))
otherTal <- as.numeric(tally(other))
totalTal <- (sineTal + lineTal + ltrTal + dnaTal + rcTal + pleTal + unknownTal)

sineTal.nest <- as.numeric(tally(sine.nest))
lineTal.nest <- as.numeric(tally(line.nest))
ltrTal.nest <- as.numeric(tally(ltr.nest))
dnaTal.nest <- as.numeric(tally(dna.nest))
rcTal.nest <- as.numeric(tally(rc.nest))
pleTal.nest <- as.numeric(tally(ple.nest))
unknownTal.nest <- as.numeric(tally(unknown.nest))
otherTal.nest <- as.numeric(tally(other.nest))

# Calculate seq percentage
genomeSize <- as.numeric(genomeSize)
sinePer <- (sineCov / genomeSize) * 100
linePer <- (lineCov / genomeSize) * 100
ltrPer <- (ltrCov / genomeSize) * 100
dnaPer <- (dnaCov / genomeSize) * 100
rcPer <- (rcCov / genomeSize) * 100
plePer <- (pleCov / genomeSize) * 100
unknownPer <- (unknownCov / genomeSize) * 100
otherPer <- (otherCov / genomeSize) * 100
totalPer <- (sinePer + linePer + ltrPer + dnaPer + rcPer + plePer + unknownPer)
unmaskPer <- (100 - totalPer)

sinePer.nest <- (sineCov.nest / genomeSize) * 100
linePer.nest <- (lineCov.nest / genomeSize) * 100
ltrPer.nest <- (ltrCov.nest / genomeSize) * 100
dnaPer.nest <- (dnaCov.nest / genomeSize) * 100
rcPer.nest <- (rcCov.nest / genomeSize) * 100
plePer.nest <- (pleCov.nest / genomeSize) * 100
unknownPer.nest <- (unknownCov.nest / genomeSize) * 100
otherPer.nest <- (otherCov.nest / genomeSize) * 100

# Make Table  
eleTyp <- c("SINE", "LINE", "LTR", "DNA", "RC", "PLE", "Unknown", "Total Interspersed Repeat", "Other", "Unmasked",
            "SINE-nested", "LINE-nested", "LTR-nested", "DNA-nested", "RC-nested", "PLE-nested", "Unknown-nested", "Other-nested")
numEle <- c(sineTal, lineTal, ltrTal, dnaTal, rcTal, pleTal, unknownTal, totalTal, otherTal, as.character("notApplicable"),
            sineTal.nest, lineTal.nest, ltrTal.nest, dnaTal.nest, rcTal.nest, pleTal.nest, unknownTal.nest, otherTal.nest)
lenOcc <- c(sineCov, lineCov, ltrCov, dnaCov, rcCov, pleCov, unknownCov, totalCov, otherCov, unmaskCov,
            sineCov.nest, lineCov.nest, ltrCov.nest, dnaCov.nest, rcCov.nest, pleCov.nest, unknownCov.nest, otherCov.nest)
perSeq <- c(sinePer, linePer, ltrPer, dnaPer, rcPer, plePer, unknownPer, totalPer, otherPer, unmaskPer,
            sinePer.nest, linePer.nest, ltrPer.nest, dnaPer.nest, rcPer.nest, plePer.nest, unknownPer.nest, otherPer.nest)
size <- c(rep(genomeSize, 18))

summary <- data.frame(eleTyp, numEle, lenOcc, perSeq, size)
colnames(summary) <- c("Family", "Number of Elements", "Length Occupied", "Percentage of Sequence", "GenomeSize")

# Print table
write.table(summary, output2, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

filteredRepeatsOut <- mergedRepeatsOut
filteredRepeatsOut$length <- abs(filteredRepeatsOut$End - filteredRepeatsOut$Start)
if (lowend == "yes") {
  filteredRepeatsOut <- filteredRepeatsOut[filteredRepeatsOut$length > 100,]
}

# if end coordinate is before start, switch
filteredRepeatsOut2 <- filteredRepeatsOut[,1:6] %>%
  mutate(nStart = ifelse(End < Start, End, Start),
         nEnd = ifelse(End < Start, Start, End)) %>%
  mutate(Start = nStart,
         End = nEnd) %>%
  select(! c(nStart, nEnd))

write.table(filteredRepeatsOut2, file = filtBed, quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)


# Generate RM style summary table for large Repeats

sine <- mergedRepeats %>% filter(str_detect(Repeat, "SINE")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
line <- mergedRepeats %>% filter(str_detect(Repeat, "LINE")) %>% filter(! str_detect(Repeat, "Penelope")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
ltr <- mergedRepeats %>% filter(str_detect(Repeat, "LTR")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
dna <- mergedRepeats %>% filter(str_detect(Repeat, "DNA")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
rc <- mergedRepeats %>% filter(str_detect(Repeat, "RC")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
ple <- mergedRepeats %>% filter(str_detect(Repeat, "Penelope")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
unknown <- mergedRepeats %>% filter(str_detect(Repeat, "Unknown")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
other <- mergedRepeats %>% filter(str_detect(Repeat, "ARTEFACT|Low|Retrop|Satellite|Simple|RNA")) %>% filter(! str_detect(Desc, "FULLY_NESTED")) %>% distinct()
other <- other[! other$Repeat == "SINE",]

sine.nest <- mergedRepeats %>% filter(str_detect(Repeat, "SINE")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
line.nest <- mergedRepeats %>% filter(str_detect(Repeat, "LINE")) %>% filter(str_detect(Repeat, "Penelope")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
ltr.nest <- mergedRepeats %>% filter(str_detect(Repeat, "LTR")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
dna.nest <- mergedRepeats %>% filter(str_detect(Repeat, "DNA")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
rc.nest <- mergedRepeats %>% filter(str_detect(Repeat, "RC")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
ple.nest <- mergedRepeats %>% filter(str_detect(Repeat, "Penelope")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
unknown.nest <- mergedRepeats %>% filter(str_detect(Repeat, "Unknown")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
other.nest <- mergedRepeats %>% filter(str_detect(Repeat, "ARTEFACT|Low|Retrop|Satellite|Simple|RNA")) %>% filter(str_detect(Desc, "FULLY_NESTED")) %>% distinct()
other.nest <- other.nest[! other.nest$Repeat == "SINE",]

# Total bp coverage of each class
sineCov <- sum(as.numeric(as.character(sine$size)))
lineCov <- sum(as.numeric(as.character(line$size)))
ltrCov <- sum(as.numeric(as.character(ltr$size)))
dnaCov <- sum(as.numeric(as.character(dna$size)))
rcCov <- sum(as.numeric(as.character(rc$size)))
pleCov <- sum(as.numeric(as.character(ple$size)))
unknownCov <- sum(as.numeric(as.character(unknown$size)))
otherCov <- sum(as.numeric(as.character(other$size)))
totalCov <- (sineCov + lineCov + ltrCov + dnaCov + rcCov + pleCov + unknownCov) 
unmaskCov <- (as.numeric(genomeSize) - as.numeric(totalCov))

sineCov.nest <- sum(as.numeric(as.character(sine.nest$size)))
lineCov.nest <- sum(as.numeric(as.character(line.nest$size)))
ltrCov.nest <- sum(as.numeric(as.character(ltr.nest$size)))
dnaCov.nest <- sum(as.numeric(as.character(dna.nest$size)))
rcCov.nest <- sum(as.numeric(as.character(rc.nest$size)))
pleCov.nest <- sum(as.numeric(as.character(ple.nest$size)))
unknownCov.nest <- sum(as.numeric(as.character(unknown.nest$size)))
otherCov.nest <- sum(as.numeric(as.character(other.nest$size)))

# Number of elements
sineTal <- as.numeric(tally(sine))
lineTal <- as.numeric(tally(line))
ltrTal <- as.numeric(tally(ltr))
dnaTal <- as.numeric(tally(dna))
rcTal <- as.numeric(tally(rc))
pleTal <- as.numeric(tally(ple))
unknownTal <- as.numeric(tally(unknown))
otherTal <- as.numeric(tally(other))
totalTal <- (sineTal + lineTal + ltrTal + dnaTal + rcTal + pleTal + unknownTal)

sineTal.nest <- as.numeric(tally(sine.nest))
lineTal.nest <- as.numeric(tally(line.nest))
ltrTal.nest <- as.numeric(tally(ltr.nest))
dnaTal.nest <- as.numeric(tally(dna.nest))
rcTal.nest <- as.numeric(tally(rc.nest))
pleTal.nest <- as.numeric(tally(ple.nest))
unknownTal.nest <- as.numeric(tally(unknown.nest))
otherTal.nest <- as.numeric(tally(other.nest))

# Calculate seq percentage
genomeSize <- as.numeric(genomeSize)
sinePer <- (sineCov / genomeSize) * 100
linePer <- (lineCov / genomeSize) * 100
ltrPer <- (ltrCov / genomeSize) * 100
dnaPer <- (dnaCov / genomeSize) * 100
rcPer <- (rcCov / genomeSize) * 100
plePer <- (pleCov / genomeSize) * 100
unknownPer <- (unknownCov / genomeSize) * 100
otherPer <- (otherCov / genomeSize) * 100
totalPer <- (sinePer + linePer + ltrPer + dnaPer + rcPer + plePer + unknownPer)
unmaskPer <- (100 - totalPer)

sinePer.nest <- (sineCov.nest / genomeSize) * 100
linePer.nest <- (lineCov.nest / genomeSize) * 100
ltrPer.nest <- (ltrCov.nest / genomeSize) * 100
dnaPer.nest <- (dnaCov.nest / genomeSize) * 100
rcPer.nest <- (rcCov.nest / genomeSize) * 100
plePer.nest <- (pleCov.nest / genomeSize) * 100
unknownPer.nest <- (unknownCov.nest / genomeSize) * 100
otherPer.nest <- (otherCov.nest / genomeSize) * 100

# Make Table  
eleTyp <- c("SINE", "LINE", "LTR", "DNA", "RC", "PLE", "Unknown", "Total Interspersed Repeat", "Other", "Unmasked",
            "SINE-nested", "LINE-nested", "LTR-nested", "DNA-nested", "RC-nested", "PLE-nested", "Unknown-nested", "Other-nested")
numEle <- c(sineTal, lineTal, ltrTal, dnaTal, rcTal, pleTal, unknownTal, totalTal, otherTal, as.character("notApplicable"),
            sineTal.nest, lineTal.nest, ltrTal.nest, dnaTal.nest, rcTal.nest, pleTal.nest, unknownTal.nest, otherTal.nest)
lenOcc <- c(sineCov, lineCov, ltrCov, dnaCov, rcCov, pleCov, unknownCov, totalCov, otherCov, unmaskCov,
            sineCov.nest, lineCov.nest, ltrCov.nest, dnaCov.nest, rcCov.nest, pleCov.nest, unknownCov.nest, otherCov.nest)
perSeq <- c(sinePer, linePer, ltrPer, dnaPer, rcPer, plePer, unknownPer, totalPer, otherPer, unmaskPer,
            sinePer.nest, linePer.nest, ltrPer.nest, dnaPer.nest, rcPer.nest, plePer.nest, unknownPer.nest, otherPer.nest)
size <- c(rep(genomeSize, 18))

summary2 <- data.frame(eleTyp, numEle, lenOcc, perSeq, size)
colnames(summary2) <- c("Family", "Number of Elements", "Length Occupied", "Percentage of Sequence", "GenomeSize")

# Print table
write.table(summary2, file = filtSum, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
