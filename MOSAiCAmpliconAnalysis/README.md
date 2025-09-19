# MosaicAmpliconAnalysis

This repository contains a structured pipeline for the analysis of amplicon sequencing data from the MOSAiC expedition. 
It has been refactored for clarity, reproducibility, and scientific publication standards.

---

## 🗂 Repository Structure

```
MosaicAmpliconAnalysis/
├── data/             # Input data or instructions to download datasets
  ├── 16S/            # Output of 16S dada2 analysis
    ├── Fastqs/       # Place to put 16S fastqs downloaded from ENA
  ├── 18S/            # Output of 18S dada2 analysis
    ├── Fastqs/       # Place to put 18S fastqs downloaded from ENA
├── output/           # Output directory for generated results
  ├── Microeco/       # Microeco statistics output
├── scripts/          # Analysis scripts (ordered and modular)
  ├── run_all.R       # Master script for full pipeline execution
├── figures/          # Final plots and visualizations
├── requirements.txt  # R dependencies
├── install.R         # Script to install requirements
├── LICENSE           # Licensing information
├── CITATION.cff      # Citation information
└── README.md         # Project overview and instructions
```

---

## 🔧 Requirements

# R version: >= 4.2

# Required CRAN packages
tidyverse
dplyr
plyr
tidyr
vegan
data.table
ggplot2
ape
readxl
xlsx
stringr
grid
magrittr
paletteer
RColorBrewer
zoo

microeco
# Packages required for microeco see: https://chiliubio.github.io/microeco_tutorial/intro.html#dependence

# Bioconductor package
BiocManager
dada2
microbiome
Biostrings
shortRead
phyloseq

You can install the dependencies via:

```r
install.packages(c("tidyverse", "dplyr", "plyr, "tydyr, "vegan", "data.table", "ggplot2", "ape", "readxl", "xlsx", "stringr", "grid", "magrittr", "paletteer", "RColorBrewer", "zoo", "microeco"))
```

```r
#Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("microbiome","dada2","Phyloseq", "Biostrings","shortRead"))
```

---

## ▶️ Usage

To reproduce the analysis pipeline:

1. Place raw input data into the `data/` directory.
2. Run scripts in the defined order from the `scripts/` folder, or use the master script:

```bash
Rscript scripts/run_all.R
```

Individual scripts:
```bash
Rscript scripts/01_generate_asv_tables_prokaryotes.R
Rscript scripts/02_generate_asv_tables_eukaryotes.R
...
```

Each script has header comments describing:
- Purpose
- Input files
- Output files

---

## 📊 Output

The `output/` and `figures/` directories will be populated as the scripts run.
All results used in the publication will be generated there.

---

## 🧾 Citation

If you use this repository, please cite:

> [Insert Paper Title]  
> [Authors]  
> [Journal, Year]  
> DOI: [Insert DOI here]

---

## 📄 License

This project is licensed under the Creative Commons Attribution 4.0 International (CC-BY 4.0).
