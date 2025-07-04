suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(magrittr))

args <- commandArgs()
# print(args)
bedIn <- args[6]
gffIn <- args[7]
gffOut <- args[8]

# test
#bedIn <- "~/projects/monarch/sortingSINEs/repeatMasker/merge/danausPlexippus.filteredRepeats.bed"
#gffIn <- "/home/toby/projects/monarch/sortingSINEs/repeatMasker/merge/remerge/danausPlexippus.rmerge.gff.filtered"
#gffOut <- "/home/toby/projects/monarch/sortingSINEs/repeatMasker/merge/remerge/danausPlexippus.filteredRepeats.gff"

# step 1 - read in tables

tab1 <- read.table(bedIn)
tab2 <- read.table(gffIn, sep = "\t")

tab1$id <- paste(tab1$V1, tab1$V2, tab1$V3, tab1$V4, tab1$V5, sep = "_")
tab2$id <- paste(tab2$V1, (tab2$V4 - 1), tab2$V5, tab2$V3, tab2$V6, sep = "_")

tab2 <- tab2 %>% distinct()

fil <- tab2[tab2$id %in% tab1$id,]
fil <- fil[,c(1:9)]
fil$V8 <- "." #set phase to be Geneious-Friendly

for (i in 1:length(fil$V1)) {
  if (fil$V5[i] < fil$V4[i]) {
    start <- fil$V5[i]
    end <- fil$V4[i]
    fil$V4[i] <- start
    fil$V5[i] <- end
  }
}

fil <- fil %>%
  mutate(V2 = case_when(V9 %like% "^ID=HLE*" ~ "HELIANO",
                        .default = "Earl_Grey"),
         V9 = toupper(V9),
         ID = paste0(sub(".*;ID=", "ID=", sub(";SHORTTE.*$", "", V9)), "_", row_number()),
         Name = paste0(sub(".*;ID=", "Name=", sub(";SHORTTE.*$", "", V9))),
         V9 = paste0(ID, ";", Name, ";", sub(";ID=.*", "", V9), ";", sub(".*SHORTTE", "SHORTTE", V9))) %>%
  select(-c(ID, Name))

write.table(fil, gffOut, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
