# Install R Packages
#
# NOTE: This script is only for manual / non-pixi (non-conda) installs.
# If you provision EarlGrey with pixi, these R packages are already installed
# as conda packages by `pixi install` (declared in pixi/pixi.toml as
# r-optparse, r-ape, bioconductor-plyranges, bioconductor-bsgenome, plus the
# rest of the plotting stack) and pinned in pixi/pixi.lock. In that case you do
# NOT need to run this script.

install.packages(c("BiocManager","optparse", "ape", "optparse"), repos="https://www.stats.bris.ac.uk/R/")
BiocManager::install(c("plyranges", "BSgenome"), update = FALSE)
