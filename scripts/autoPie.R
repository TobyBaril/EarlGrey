# load libraries

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(kableExtra))

# set options

options(scipen = 100, stringsAsFactors = FALSE)

# get inputs


args <- commandArgs()
# print(args)
inBed <- args[6]
inGff <- args[7]
gen <- args[8]
savePie <- args[9]
saveTab <- args[10]
## bed file
## genome size
## pieName

# read table

input <- read.table(inBed, header = FALSE, sep = "\t")
input2 <- read.table(inGff, header = FALSE, sep = "\t")
gen <- as.numeric(gen)

# set colnames

colnames(input) <- c("scaf", "start", "end", "classif", "score", "strand")
colnames(input2) <- c("scaf", "method", "classif", "start", "end", "score", "strand", "dot", "attributes")

# simple classification

input$tclassif <- input$classif
input$tclassif <- gsub("^DNA.*", "DNA", input$tclassif)
input$tclassif <- gsub("^RC.*", "Rolling Circle", input$tclassif)
input$tclassif <- gsub(".*Penelope|^PLE.*", "Penelope", input$tclassif)
input$tclassif <- gsub("^LINE.*", "LINE", input$tclassif)
input$tclassif <- gsub("^SINE.*", "SINE", input$tclassif)
input$tclassif <- gsub("^LTR.*", "LTR", input$tclassif)
input$tclassif <- gsub("^Unknown.*|Retroposon.*|Unspecified.*", "Unclassified", input$tclassif)
input$tclassif <- gsub(".*RNA.*|^Satellite.*|^Simple_repeat.*|^Low_complexity.*|^ARTEFACT.*|^repeat.*|^Other.*", "Other (Simple Repeat, Microsatellite, RNA)", input$tclassif)

input2$tclassif <- input2$classif
input2$tclassif <- gsub("^DNA.*", "DNA", input2$tclassif)
input2$tclassif <- gsub("^RC.*", "Rolling Circle", input2$tclassif)
input2$tclassif <- gsub(".*Penelope|^PLE.*", "Penelope", input2$tclassif)
input2$tclassif <- gsub("^LINE.*", "LINE", input2$tclassif)
input2$tclassif <- gsub("^SINE.*", "SINE", input2$tclassif)
input2$tclassif <- gsub("^LTR.*", "LTR", input2$tclassif)
input2$tclassif <- gsub("^Unknown.*|Retroposon.*|Unspecified.*", "Unclassified", input2$tclassif)
input2$tclassif <- gsub(".*RNA.*|^Satellite.*|^Simple_repeat.*|^Low_complexity.*|^ARTEFACT.*|^repeat.*|^Other.*", "Other (Simple Repeat, Microsatellite, RNA)", input2$tclassif)


# set colours

# https://oobrien.com/2012/01/tube-colours/

# DNA = #E32017
# Rolling Circle = #EE7C0E
# Penelope = #7156A5
# LINE = #0098D4
# SINE = #9B0056
# LTR = #00782A
# Other = #F3A9BB
# Unclas = #A0A5A9
# Non-Repeat = #000000

colours <- as.data.frame(matrix(c("DNA", "#E32017",
                                  "Rolling Circle", "#EE7C0E",
                                  "Penelope", "#7156A5",
                                  "LINE", "#0098D4",
                                  "SINE", "#9B0056",
                                  "LTR", "#00782A",
                                  "Other (Simple Repeat, Microsatellite, RNA)", "#F3A9BB",
                                  "Unclassified", "#A0A5A9",
                                  "Non-Repeat", "#000000"), ncol = 2, byrow = TRUE))

colnames(colours) <- c("tclassif", "colour")
col <- colours$colour
names(col) <- colours$tclassif

# summary table

pieSum <- input %>% group_by(tclassif) %>% summarise(cov = sum(abs(end - start)))
pieCou <- input %>% group_by(tclassif) %>% tally()
colnames(pieCou) <- c("tclassif", "count")
pieSum <- merge(pieSum, pieCou)

# add non-repeat
nr <- data.frame("Non-Repeat", gen - sum(abs(input$end - input$start)))
colnames(nr) <- c("tclassif", "cov")
pieSum <- dplyr::bind_rows(pieSum, nr)

# set factor order

pieSum$tclassif %<>% as.factor() %>% ordered(levels = c("DNA",
                                                        "Rolling Circle",
                                                        "Penelope",
                                                        "LINE",
                                                        "SINE",
                                                        "LTR",
                                                        "Other (Simple Repeat, Microsatellite, RNA)",
                                                        "Unclassified",
                                                        "Non-Repeat"))

# get proportions

pieSum$proportion <- pieSum$cov / gen
pieSum$percentage <- (pieSum$cov / gen) * 100
pieSum$genome_size <- gen

# plot pie

plot <- ggplot(pieSum, aes(x = as.numeric(gen/2), y = proportion, fill = tclassif, width = as.numeric(gen))) +
        geom_bar(stat = "identity", position = "fill") +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = col) +
        theme_classic() +
        theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
              strip.text = element_text(size = 16, face = "italic"), legend.text = element_text(size = 16), legend.title = element_blank())

ggsave(savePie,
       plot = plot,
       scale = 1,
       width = 297,
       height = 210,
       units = "mm",
       dpi = 300,
       limitsize = FALSE)

# add number of classifications for each of the main (ie how many families)

input2$family <- gsub(".*Name=", "", input2$attributes)
input2$family <- gsub(";.*", "", input2$family)
input2$family <- toupper(input2$family)
classCount <- input2 %>% group_by(tclassif, family) %>% tally(name = "Number_of_Copies") %>% group_by(tclassif) %>% tally(name = "Number_of_Distinct_Classifications")
pieSum <- left_join(pieSum, classCount) %>%
  ungroup() %>%
  select(-c(proportion))

pieSum$tclassif %<>% as.factor() %>% ordered(levels = c("DNA",
                                                        "Rolling Circle",
                                                        "Penelope",
                                                        "LINE",
                                                        "SINE",
                                                        "LTR",
                                                        "Other (Simple Repeat, Microsatellite, RNA)",
                                                        "Unclassified",
                                                        "Non-Repeat"))

pieSum <- pieSum %>%
  arrange(tclassif)

colnames(pieSum) <- c("TE Classification", "Coverage (bp)", "Copy Number", "% Genome Coverage", "Genome Size", "TE Family Count")
      
# generate summary of TE family abundance

familyAbundance <- input2 %>%
        mutate(name = paste0(family, "#", classif)) %>%
        group_by(name) %>%
        summarise(coverage = sum((end - start) + 1)) %>% 
        arrange(-coverage)

familyTally <- input2 %>% 
        mutate(name = paste0(family, "#", classif)) %>%
        group_by(name) %>%
        tally(name = "copy_number")

familyAbundance <- merge(familyAbundance, familyTally) %>%
        arrange(-coverage)

colnames(familyAbundance) <- c("TE Family", "Coverage (bp)", "Copy Number")

saveFam <- gsub("highLevelCount", "familyLevelCount", saveTab)

write.table(pieSum, saveTab, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(familyAbundance, saveFam, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# add a user-readable version of the final table
saveKable <- gsub(".txt", ".kable", saveTab)
kable(pieSum, format = "pipe") %>%
  writeLines(saveKable)

saveKableFam <- gsub(".txt", ".kable", saveFam)
kable(familyAbundance, format = "pipe") %>%
  writeLines(saveKableFam)
