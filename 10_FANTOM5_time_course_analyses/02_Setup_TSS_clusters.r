#!/usr/bin/env Rscript
# filepath: /mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/09_FANTOM5_time_course_analyses/01_Setup_TSS_data.R
# Setup TSS counts from CAGE data
# Author: Martin Loza
# Date: 25/12/30
# Description: Transform the CAGE bed files to TSS counts

# Change R language to English
Sys.setenv(LANGUAGE = "en")

# Init
suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(CAGEr)
    library(rtracklayer)
    library(ggplot2)
})

# Local variables 
seed = 777
date = "251230"
# Select dataset to analyze
sel_dataset <- "Monocyte-derived" # Options: "Lymphatic", "Aortic" "MCF7" "Monocyte-derived"  "Saos-2" 
# selected genome
genome = 'BSgenome.Hsapiens.UCSC.hg19'
# Define colors for strand plots
red = "#E41A1C"
blue = "#090a0bff"
# Define colors for gene types
green = "#4DAF4A"
purple = "#984EA3"
text_size = 18
width = 18.6
dot_size = 4
line_size = 1.5
dpi = 300

# Install the genome package if not already installed
if (!requireNamespace(genome, quietly = TRUE)) {
    BiocManager::install(genome)
}
# Load the genome package
library(genome, character.only = TRUE)

fantom_dir = "/mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/Data/FANTOM5/"
ctss_bed_dir = paste0(fantom_dir, "bed_files/")
ctss_clusters_out_dir = paste0(fantom_dir, "CTSS_clusters/")
plots_dir = paste0(ctss_clusters_out_dir, "Plots_preprocessing/")

# Create output directories if they don't exist
if (!dir.exists(ctss_clusters_out_dir)) {
    dir.create(ctss_clusters_out_dir, recursive = TRUE)
}
if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
}

# Load the hg19 gene annotation from UCSC
url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ensGene.gtf.gz"
hg19_transcripts <- import.gff(url, format = "gtf")

# ============================================================================
# Load and setup
# ============================================================================

# Load the metadata
metadata = read.csv(paste0(fantom_dir, "FANTOM5_timecourse_metadata_selected_series_251219.tsv"), sep = "\t", header = TRUE)
# Get the file names, they are listed in the file_name column but includes the http path. 
# We need to remove everything before the last /
metadata$file_name = sapply(metadata$file_name, function(x) {
    parts = str_split(x, "/")[[1]]
    return(parts[length(parts)])
})
# When processing the bed files we change the file names to include only the first word and the time point.
# So we need to update the file_name column accordingly. 
# We will take the first word of the data_name and append the time point.
metadata$data_name_first_word = sapply(metadata$data_name, function(x) {
    parts = str_split(x, " ")[[1]]
    return(parts[1])
})
metadata$data_name_simple = paste0(metadata$data_name_first_word, "_", metadata$time_point)
# Create underscore version for CAGE analysis (replace dashes with underscores)
metadata$data_name_underscore = gsub("-", "_", metadata$data_name_simple)
# Reorder columns
metadata <- metadata %>% select(data_name, data_name_simple, data_name_underscore, time_point, time_sequence, file_name, everything())
# Check the result
head(metadata, 2)

# ============================================================================
# Focus on one dataset
# ============================================================================

table(metadata$data_name_first_word)


metadata_sel = metadata %>% filter(data_name_first_word == sel_dataset) 

### TESTS #####
# For testing, take only the first two time points
# metadata_sel <- metadata_sel[1:2, ]


ordered_time_points = metadata_sel$time_point
print(ordered_time_points)

# Get the available files in the directory
files = list.files(ctss_bed_dir, full.names = TRUE)
# Select the files related to the selected dataset
sel_files <- files[str_detect(files, sel_dataset)] 
# Order them by time point
sel_file = c()
for (tp in ordered_time_points) {
    file_tp = sel_files[str_detect(sel_files, tp)]
    sel_file = c(sel_file, file_tp)
}
print(sel_file)

# ============================================================================
# Create and process CAGEr object
# ============================================================================

# create CAGEr object
cage <- CAGEexp(
    genomeName = genome,
    inputFiles = sel_file,
    inputFilesType = "ctss",
    sampleLabels = metadata_sel$data_name_underscore
)

# Load the CAGE files
cage <- getCTSS(cage)

# Merge samples into one 
cage_merged <- mergeSamples(cage, 
                    mergeIndex = 1:length(metadata_sel$data_name_underscore),
                    mergedSampleLabels = metadata_sel$data_name_underscore)
# Annotate the merged CAGEr object
cage_merged <- annotateCTSS(cage_merged, hg19_transcripts)

# ============================================================================
# Normalization
# ============================================================================

p <- plotReverseCumulatives(cage_merged, fitInRange = c(5, 1000))
print(p)
# Save plot
ggsave(paste0(plots_dir, sel_dataset, "_reverse_cumulatives.png"), plot = p, width = 10, height = 8, dpi = dpi)

# Normalization
cage_merged <- normalizeTagCount(cage_merged, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)

# ============================================================================
# Filter and cluster
# ============================================================================

# Filter low expression CTSS
min_tpm <- 0.1
cage_merged <- filterLowExpCTSS(cage_merged, thresholdIsTpm = TRUE, nrPassThreshold = 2, threshold = min_tpm)

# Cluster CTSS into tag clusters
cage_merged <- distclu(cage_merged, maxDist = 20, keepSingletonsAbove = min_tpm)

# Calculate cumulative CTSS distribution and quantile positions
cage_merged <- cumulativeCTSSdistribution(cage_merged, clusters = "tagClusters", useMulticore = TRUE)
cage_merged <- quantilePositions(cage_merged, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

# Plot interquantile width
p <- plotInterquantileWidth(cage_merged, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
print(p)
# Save plot
ggsave(paste0(plots_dir, sel_dataset, "_interquantile_width.png"), plot = p, width = 10, height = 8, dpi = dpi)

# Aggregate tag clusters
cage_merged <- aggregateTagClusters(cage_merged, tpmThreshold = 0.1, qLow = 0.1, qUp = 0.9, maxDist = 100)
print(cage_merged$outOfClusters / cage_merged$librarySizes * 100)

# Loop over the samples to extract the consensus clusters
consensus_list <- list()
for (i in 1:length(metadata_sel$data_name_underscore)) {
    consensus_list[[i]] <- consensusClustersGR(cage_merged, sample = metadata_sel$data_name_underscore[i]) %>% as.data.frame()
}

# Confirm that all the consensus clusters df have the same number of rows
print(sapply(consensus_list, nrow))

# Check first consensus cluster
print(head(consensus_list[[1]]))

# Select columns of interest and rename the expression column
sel_cols <- c("seqnames", "start", "end", "width", "strand", "tpm")
for (i in 1:length(consensus_list)) {
    consensus_list[[i]] <- consensus_list[[i]] %>% 
        select(all_of(sel_cols)) %>%
        dplyr::rename(!!paste0("tpm_", metadata_sel$data_name_underscore[i]) := "tpm")
}

# Merge all consensus clusters into a single data frame
consensus_merged <- consensus_list[[1]]
for (i in 2:length(consensus_list)) {
    consensus_merged <- consensus_merged %>%
        left_join(consensus_list[[i]], by = c("seqnames", "start", "end", "width", "strand"))
}

print(head(consensus_merged))
print(dim(consensus_merged))

# Save the final consensus merged dataframe
write.table(consensus_merged, paste0(ctss_clusters_out_dir, sel_dataset, "_consensus_clusters.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste0("\nConsensus clusters saved to: ", ctss_clusters_out_dir, sel_dataset, "_consensus_clusters.tsv\n"))

cat("\nAnalysis complete!\n")