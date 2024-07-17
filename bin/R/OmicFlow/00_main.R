# Load libraries ---------------------------------------------------------------
library("dplyr")
library("ggplot2")
library("tidyr")
library("dplyr")
library("biomformat")
library("phyloseq")
library("ape")
library("Biostrings")
library("patchwork")
library("vegan")
library("optparse")

# Load Functions ---------------------------------------------------------------
# Fetch current Rscript path
Rscript <- sub("--file=", "", commandArgs()[4])
current_path <- sub(basename(Rscript), "", normalizePath(Rscript))

# Data wrangling functions
source(file = paste0(current_path, "99_project-functions.R"))

# Utils
sourceDir(path = paste0(current_path, "utils"))

# Parse command line
data_00 <- parse_commandline()
outfile_path <- normalizePath(data_00$outDir)
filter_taxa_tab <- normalizePath(paste0(current_path, "../documents/database_filter"))

# Fetch additional files from preprocessing pipeline ---------------------------
preprocess_stats <- readr::read_delim(paste0(outfile_path, "/read_stats.txt"))
shannon_file <- readr::read_csv(file = paste0(outfile_path, "/alpha_rarefaction/shannon.csv"))

# Taxa to be visualized:
taxa_names <- c("Phylum", "Family", "Genus")

# Run all scripts --------------------------------------------------------------
source(file = paste0(current_path, "01_load.R"))
source(file = paste0(current_path, "02_clean.R"))
source(file = paste0(current_path, "03_data_visualization.R"))
source(file = paste0(current_path, "04_Model-UniFrac-PCoA.R"))
source(file = paste0(current_path, "05_Model-RDA.R"))
rmarkdown::render(input = paste0(current_path, "../documents/report.Rmd"),
                  output_file = paste0(outfile_path, "/report.html"))
