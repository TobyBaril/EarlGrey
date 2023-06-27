# Install R Packages

install.packages(c("BiocManager","optparse", "ape", "optparse"), repos="https://www.stats.bris.ac.uk/R/")
BiocManager::install(c("plyranges", "BSgenome"), update = FALSE)
