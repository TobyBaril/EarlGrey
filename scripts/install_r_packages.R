# Install R Packages

install.packages(c("BiocManager","optparse", "ape"), repos="https://www.stats.bris.ac.uk/R/")
BiocManager::install(c("plyranges", "BSgenome"), update = TRUE, ask = FALSE)
