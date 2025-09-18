# install.R — installs required packages for MosaicAmpliconAnalysis

install.packages(c("tidyverse", "dplyr", "plyr, "tydyr, "vegan", "data.table", "ggplot2", "ape", "readxl", "xlsx", "stringr", "grid", "magrittr", "paletteer", "RColorBrewer", "zoo", "microeco"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("microbiome","dada2","Phyloseq", "Biostrings","shortRead"))

