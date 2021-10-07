# load libraries

library(tidyverse)
library(data.table)

# set options

options(scipen = 100, stringsAsFactors = FALSE)

# get inputs


args <- commandArgs()
print(args)
inBed <- args[6]
gen <- args[7]
savePie <- args[8]
saveTab <- args[9]
## bed file
## genome size
## pieName

# read table

#inBed <- "~/projects/butterfly/importantFiles/filterRep/brenthisIno.filteredRepeats.bed"
input <- read.table(inBed, header = FALSE, sep = "\t")
gen <- as.numeric(gen)
#gen <- 406860652

# set colnames

colnames(input) <- c("scaf", "start", "end", "classif", "score", "strand")

# simple classification

input$tclassif <- input$classif
input$tclassif <- gsub("^DNA.*", "DNA", input$tclassif)
input$tclassif <- gsub("^RC.*", "Rolling Circle", input$tclassif)
input$tclassif <- gsub(".*Penelope", "Penelope", input$tclassif)
input$tclassif <- gsub("^LINE.*", "LINE", input$tclassif)
input$tclassif <- gsub("^SINE.*", "SINE", input$tclassif)
input$tclassif <- gsub("^LTR.*", "LTR", input$tclassif)
input$tclassif <- gsub("^Unknown.*|Retroposon.*|Unspecified.*", "Unclassified", input$tclassif)
input$tclassif <- gsub(".*RNA.*|^Satellite.*|^Simple_repeat.*|^Low_complexity.*|^ARTEFACT.*|^repeat.*|^Other.*", "Other (Simple Repeat, Microsatellite, RNA)", input$tclassif)

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

pieSum$percentage <- pieSum$cov / gen
pieSum$gen <- gen

# plot pie

plot <- ggplot(pieSum, aes(x = as.numeric(gen/2), y = percentage, fill = tclassif, width = as.numeric(gen))) +
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

write.table(pieSum, saveTab, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
