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

# read table

input <- read.table(inBed, header = TRUE, sep = "\t")
input2 <- read.table(inGff, header = FALSE, sep = "\t")
gen <- as.numeric(gen)

# set colnames
colnames(input2) <- c("scaf", "method", "classif", "start", "end", "score", "strand", "dot", "attributes")

# simple classification
input$tclassif <- input$Family
input$tclassif <- gsub("^DNA.*", "DNA", input$tclassif)
input$tclassif <- gsub("^RC.*", "Rolling Circle", input$tclassif)
input$tclassif <- gsub(".*Penelope|^PLE.*", "Penelope", input$tclassif)
input$tclassif <- gsub("^LINE.*", "LINE", input$tclassif)
input$tclassif <- gsub("^SINE.*", "SINE", input$tclassif)
input$tclassif <- gsub("^LTR.*", "LTR", input$tclassif)
input$tclassif <- gsub("^Unknown.*|Retroposon.*|Unspecified.*|^unknown.*", "Unclassified", input$tclassif)
input$tclassif <- gsub(".*RNA.*|^Satellite.*|^Simple_repeat.*|^Low_complexity.*|^ARTEFACT.*|^repeat.*|^Other.*|^Microsatellite.*", "Other (Simple Repeat, Microsatellite, RNA)", input$tclassif)
input <- input %>%
  mutate(tclassif = case_when(input$Family %like% "-nested" ~ paste0(tclassif, "-nested"),
                              .default = tclassif))
input$tclassif <- gsub("Unmasked", "Non-Repeat", input$tclassif)

input2$tclassif <- input2$classif
input2$tclassif <- gsub("^DNA.*", "DNA", input2$tclassif)
input2$tclassif <- gsub("^RC.*", "Rolling Circle", input2$tclassif)
input2$tclassif <- gsub(".*Penelope|^PLE.*", "Penelope", input2$tclassif)
input2$tclassif <- gsub("^LINE.*", "LINE", input2$tclassif)
input2$tclassif <- gsub("^SINE.*", "SINE", input2$tclassif)
input2$tclassif <- gsub("^LTR.*", "LTR", input2$tclassif)
input2$tclassif <- gsub("^Unknown.*|Retroposon.*|Unspecified.*|^unknown.*", "Unclassified", input2$tclassif)
input2$tclassif <- gsub(".*RNA.*|^Satellite.*|^Simple_repeat.*|^Low_complexity.*|^ARTEFACT.*|^repeat.*|^Other.*|^Microsatellite.*", "Other (Simple Repeat, Microsatellite, RNA)", input2$tclassif)
input2 <- input2 %>%
  mutate(tclassif = case_when(attributes %like% "FULLY_NESTED" ~ paste0(tclassif, "-nested"),
                              .default = tclassif))


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

# use summary table generated previously...
pieCou <- input %>% 
  select(tclassif, Number.of.Elements)
colnames(pieCou) <- c("tclassif", "count")
pieSum <- input

# set factor order
pieSum$tclassif %<>% as.factor() %>% ordered(levels = c("DNA",
                                                        "DNA-nested",
                                                        "Rolling Circle",
                                                        "Rolling Circle-nested",
                                                        "Penelope",
                                                        "Penelope-nested",
                                                        "LINE",
                                                        "LINE-nested",
                                                        "SINE",
                                                        "SINE-nested",
                                                        "LTR",
                                                        "LTR-nested",
                                                        "Other (Simple Repeat, Microsatellite, RNA)",
                                                        "Other (Simple Repeat, Microsatellite, RNA)-nested",
                                                        "Unclassified",
                                                        "Unclassified-nested",
                                                        "Total Interspersed Repeat",
                                                        "Non-Repeat"))

# get proportions

pieSum$proportion <- pieSum$Percentage.of.Sequence / 100
pieSum$percentage <- pieSum$Percentage.of.Sequence
pieSum$genome_size <- pieSum$GenomeSize

# plot pie

plot <- ggplot(pieSum %>%
                 filter(! tclassif %like% "-nested",
                        ! tclassif %in% c("Total Interspersed Repeat")), aes(x = as.numeric(gen/2), y = proportion, fill = tclassif, width = as.numeric(gen))) +
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



# generate summary of TE family abundance

familyAbundance <- input2 %>%
  mutate(family = gsub(".*Name=", "", gsub(";TSTART.*", "", attributes))) %>%
  mutate(name = paste0(family, "#", classif)) %>%
  group_by(name) %>%
  summarise(coverage = sum((end - start) + 1)) %>% 
  arrange(-coverage)

familyTally <- input2 %>% 
  mutate(family = gsub(".*Name=", "", gsub(";TSTART.*", "", attributes))) %>%
  mutate(name = paste0(family, "#", classif)) %>%
  group_by(name) %>%
  tally(name = "copy_number")

familyAbundance <- merge(familyAbundance, familyTally) %>%
  arrange(-coverage)

colnames(familyAbundance) <- c("TE Family", "Coverage (bp)", "Copy Number")

saveFam <- gsub("highLevelCount", "familyLevelCount", saveTab)

classTally <- input2 %>% 
  group_by(tclassif) %>%
  tally(name = "Copy.Number") %>%
  rename("tclassif" = "Family")

classFamilyCount <- input2 %>%
  mutate(family = gsub(".*Name=", "", gsub(";TSTART.*", "", attributes))) %>%
  select(family, tclassif) %>%
  distinct() %>%
  group_by(tclassif) %>%
  tally(name = "Family.Count") %>%
  rename("tclassif" = "Family")

pieSum <- pieSum %>%
  ungroup() %>%
  select(-c(Family)) %>%
  rename("tclassif" = "Family") %>%
  left_join(classTally) %>%
  left_join(classFamilyCount)

# finalise table
pieSum$Family %<>% as.factor() %>% ordered(levels = c("DNA",
                                                      "DNA-nested",
                                                      "Rolling Circle",
                                                      "Rolling Circle-nested",
                                                      "Penelope",
                                                      "Penelope-nested",
                                                      "LINE",
                                                      "LINE-nested",
                                                      "SINE",
                                                      "SINE-nested",
                                                      "LTR",
                                                      "LTR-nested",
                                                      "Other (Simple Repeat, Microsatellite, RNA)",
                                                      "Other (Simple Repeat, Microsatellite, RNA)-nested",
                                                      "Unclassified",
                                                      "Unclassified-nested",
                                                      "Total Interspersed Repeat",
                                                      "Non-Repeat"))

pieSum <- pieSum %>%
  arrange(Family) %>%
  select(-c(percentage, genome_size, proportion, Copy.Number)) %>%
  relocate(Family, Length.Occupied, Number.of.Elements, Percentage.of.Sequence, GenomeSize, Family.Count)

colnames(pieSum) <- c("TE Classification", "Coverage (bp)", "Copy Number", "% Genome Coverage", "Genome Size", "TE Family Count")

write.table(pieSum, saveTab, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(familyAbundance, saveFam, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# add a user-readable version of the final table
saveKable <- gsub(".txt", ".kable", saveTab)
kable(pieSum, format = "pipe") %>%
  writeLines(saveKable)

saveKableFam <- gsub(".txt", ".kable", saveFam)
kable(familyAbundance, format = "pipe") %>%
  writeLines(saveKableFam)
