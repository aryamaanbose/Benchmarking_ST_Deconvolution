# Load the renv package and restore the environment from the renv.lock file
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::restore()

# Install the 'remotes' package if it's not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Define GitHub packages with their repositories
github_packages <- c(
  "YingMa0107/CARD", 
  "xuranw/MuSiC",
  "SONGDONGYUAN1994/scDesign3",
  "IOBR/IOBR"
)

# Install GitHub packages
for (repo in github_packages) {
  remotes::install_github(repo, dependencies = TRUE)
}


#renv::init()

# List of CRAN and Bioconductor packages with specific versions
packages <- c(
  "dplyr" = "1.1.4",
  "Seurat" = "5.1.0",
  "patchwork" = "1.2.0",
  "ggplot2" = "3.5.1",
  "imager" = "1.0.2",
  "readr" = "2.1.5",
  "SpatialExperiment" = "1.10.0",
  "magick" = "2.8.3",
  "cowplot" = "1.1.3",
  "jpeg" = "0.1-10",
  "data.table" = "1.15.4",
  "tidyverse" = "2.0.0",
  "matrixStats" = "1.3.0",
  "scran" = "1.28.2",
  "gridExtra" = "2.3",
  "reshape2" = "1.4.4",
  "here" = "1.0.1",
  "mclust" = "6.1.1",
  "TOAST" = "1.14.0",
  "spdep" = "1.3-5"
)

# Install CRAN and Bioconductor packages with specific versions
for (pkg in names(packages)) {
  version <- packages[pkg]
  renv::install(sprintf("%s@%s", pkg, version))
}

# Take a snapshot of the environment
#renv::snapshot()

# Load the packages
lapply(names(packages), library, character.only = TRUE)
