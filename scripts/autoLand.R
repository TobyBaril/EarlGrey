# load libraries

library(tidyverse)
library(data.table)
library(magrittr)

# set options

options(scipen = 100, stringsAsFactors = FALSE)

# get inputs

args <- commandArgs()
print(args)
inDist <- args[6]
gen <- args[7]
spec <- args[8]
saveLand <- args[9]
## bed file
## genome size
## landscape name

# test
#inDist = "~/projects/monarch/datingInsertions/danausPlexippus.divsum.unaltered"
#gen = 248676414
#spec = "test"

# read table

dist <- as.data.frame(t(fread(inDist, skip = "Div", header = TRUE)))
dist$species <- spec
dist$type <- rownames(dist)

# fix shitty columns

dist <- dist[! dist$type == "X",]
colnames(dist) <- c(dist[1,1:71], "species", "type")
dist <- dist[! dist$type == "Div",]
rownames(dist) <- 1:nrow(dist)

# get sidple TE classification levels

dist$classif <- dist$type
dist$classif <- gsub("^DNA.*", "DNA", dist$classif)
dist$classif <- gsub("^LTR.*", "LTR", dist$classif)
dist$classif <- gsub(".*Penelope.*", "Penelope", dist$classif)
dist$classif <- gsub("^LINE.*", "LINE", dist$classif)
dist$classif <- gsub("^SINE.*", "SINE", dist$classif)
dist$classif <- gsub("^RC.*", "Rolling Circle", dist$classif)
dist$classif <- gsub("^Unknown.*|Retroposon.*|Unspecified.*", "Unclassified", dist$classif)
dist$classif <- gsub(".*RNA.*|^Satellite.*|^Simple_repeat.*|^Low_complexity.*|^ARTEFACT.*|^Other.*|^Segmenta.*", "Other (Simple Repeat, Microsatellite, RNA)", dist$classif)


# long format

piv <- dist[,c(1:72, 74)] %>% pivot_longer(!c("species", "classif"), names_to = "Divergence", values_to = "count") %>% group_by(species, classif, Divergence) %>% summarise(Count = sum(count))

piv$Divergence <- as.numeric(piv$Divergence)

colourOrder <- c("DNA",
                 "Rolling Circle",
                 "Penelope",
                 "LTR",
                 "LINE",
                 "SINE",
                 "Other (Simple Repeat, Microsatellite, RNA)",
                 "Unclassified")

piv$classif %<>% as.factor() %>% ordered(levels = colourOrder)

# add genome size and scale number of elements

genomeSizes <- as.data.frame(matrix(c(spec, gen),
                                    byrow = TRUE, ncol = 2))

colnames(genomeSizes) <- c("species", "genomeSize")

piv <- merge(piv, genomeSizes)

piv$genomeSize <- as.numeric(piv$genomeSize)

piv$proportion <- piv$Count / piv$genomeSize

for (i in 1:nrow(piv)) {
  if (piv$Count[i] < 0) {
    piv$Count[i] <- 0
  }
}

for (i in 1:length(piv$proportion)){
  if (piv$proportion[i] < 0) {
    piv$proportion[i] <- 0
  }
  if (i %% 100 == 0) { #if x is divisible by 100
    cat(sprintf("Completed: %s\r", i)) #print Completed: %s\r
  }
}


# summarise

piv <- piv %>% distinct()
piv2 <- piv %>% group_by(species, Divergence, classif) %>% summarise(proportion = sum(proportion))

# set plot limits

limits <- piv2 %>% group_by(Divergence) %>% summarise(perc = sum(proportion * 100))
ylimit <- max(limits$perc) * 1.1

# plot

plot <- ggplot(piv2[! piv2$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
  geom_bar(stat = "identity") +
  ylim(0, ylimit) +
  scale_x_reverse() +
  scale_fill_manual(values = c("#E32017", "#EE7C0E", "#7156A5", "#00782A", "#0098D4", "#9B0056","#A0A5A9")) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.y = element_text(face = "italic", angle = 0),
        axis.text.y = element_text(size = 8)) +
  labs(y = "Percentage of Genome", x = "Kimura Distance (CpG Adjusted)", fill = "")

ggsave(saveLand,
       plot = plot,
       scale = 1,
       width = 297,
       height = 210,
       units = "mm",
       dpi = 300,
       limitsize = FALSE)