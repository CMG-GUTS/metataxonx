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
#outfile_path <- normalizePath(data_00$outDir)
outfile_path <- normalizePath(getwd())
filter_taxa_tab <- normalizePath(paste0(current_path, "../documents/database_filter"))


# Fetch additional files from preprocessing pipeline ---------------------------
read_stats_path <- paste0(outfile_path, "/read_stats.txt")
shannon_path <- paste0(outfile_path, "/alpha_rarefaction/shannon.csv")
sankeyplot <- paste0(outfile_path, "/biotaviz_sankey_prepfile-AverageAllSamples.png")
if (file.exists(read_stats_path)) {
    preprocess_stats <- readr::read_delim(read_stats_path)
} else {
    preprocess_stats <- NULL
}

if (file.exists(shannon_path)) {
    shannon_file <- readr::read_csv(file = shannon_path)
} else {
    shannon_file <- NULL
}

if (file.exists(sankeyplot)) {
    sankeyplot <- sankeyplot
} else {
    sankeyplot <- NULL
}

# collect untrimmed html page
if (file.exists(paste0(outfile_path, "/untrimmed/multiqc_report.html"))) {
    untrimmed_multiqc = paste0(outfile_path, "/untrimmed/multiqc_report.html")
} else {
    untrimmed_multiqc <- NULL
}
# collect trimmed html page
if (file.exists(paste0(outfile_path, "/trimmed/multiqc_report.html"))) {
    trimmed_multiqc = paste0(outfile_path, "/trimmed/multiqc_report.html")
} else {
    trimmed_multiqc <- NULL
}
# collect pear html page
if (file.exists(paste0(outfile_path, "/multiqc_pear/multiqc_report.html"))) {
    pear_multiqc = paste0(outfile_path, "/multiqc_pear/multiqc_report.html")
} else {
    pear_multiqc <- NULL
}

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
