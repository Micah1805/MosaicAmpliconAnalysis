# install.R — installs required packages for MosaicAmpliconAnalysis

install.packages(c("tidyverse", "phyloseq", "vegan", "data.table", "ggplot2", "ape"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("microbiome","dada2"))
