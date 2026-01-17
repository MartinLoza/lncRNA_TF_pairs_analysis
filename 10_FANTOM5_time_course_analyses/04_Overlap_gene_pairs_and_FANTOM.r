# FANTOM5 time course analysis
# Overlap gene pairs and FANTOM TSSs
# Author: Martin Loza
# Date: 25/12/30

# Change R language to English
Sys.setenv(LANGUAGE = "en")

# Init
suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(GenomicRanges)
})

# Local variables 
date = "251230"

# Define datasets to process
datasets <- c("Lymphatic", "Aortic", "MCF7", "Monocyte-derived", "Saos-2")

gene_pairs_file = "/mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/09_FANTOM5_time_course_analyses/Results/Ranked_lncRNA_TF_gene_pairs_hg19_TSS_251230.tsv"
fantom_tss_dir = "/mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/Data/FANTOM5/CTSS_clusters/"
out_dir = "/mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/09_FANTOM5_time_course_analyses/Results/Gene_pairs_FANTOM_expression/"

# Local Functions
CreateGRangesObject <- function(ranges = NULL,
                seqnames = "chr", start = "start", end = "end",
                expand = NULL, metadata = NULL, retain_md = FALSE){
  gr_object <- GRanges(seqnames = ranges[[seqnames]], 
             ranges = IRanges(start = ranges[[start]], end = ranges[[end]]))
  
  if(!is.null(expand)){
    gr_object <- gr_object + expand
  }

  if(!is.null(metadata)){
    mcols(gr_object) <- metadata
  } else if (retain_md) {
    mcols(gr_object) <- ranges[, !(names(ranges) %in% c(seqnames, start, end))]
  }

  return(gr_object)
}

# Load the gene pairs annotated with hg19 TSS information 
cat("Loading gene pairs data...\n")
gene_pairs_ranked <- read.table(gene_pairs_file, header = TRUE, sep = "\t", 
                            stringsAsFactors = FALSE, check.names = FALSE)
cat("Number of gene pairs: ", nrow(gene_pairs_ranked), "\n\n")

# Create unique TSS from gene pairs
gene_tss_hg19 <- rbind(
  gene_pairs_ranked %>%
    dplyr::select(seqnames = lncrna_tss_hg19_chromosome, 
                  start = lncrna_tss_hg19_start, 
                  end = lncrna_tss_hg19_end, 
                  gene_tss_id = lncRNA_tss_id),
  gene_pairs_ranked %>%
    dplyr::select(seqnames = tf_tss_hg19_chromosome, 
                  start = tf_tss_hg19_start, 
                  end = tf_tss_hg19_end, 
                  gene_tss_id = TF_tss_id)
) %>%
    distinct()

# Create GenomicRanges object for gene TSS 
gene_tss_gr <- CreateGRangesObject(ranges = gene_tss_hg19,
                                 seqnames = "seqnames",
                                 start = "start",
                                 end = "end",
                                 metadata = gene_tss_hg19["gene_tss_id"],
                                 retain_md = TRUE)

# Loop over datasets
for(sel_dataset in datasets) {
  cat("========================================\n")
  cat("Processing dataset: ", sel_dataset, "\n")
  cat("========================================\n")
  
  # Load the TSS clusters from FANTOM5 for the selected dataset
  fantom_tss_file <- paste0(fantom_tss_dir, sel_dataset, "_consensus_clusters.tsv")
  
  if(!file.exists(fantom_tss_file)) {
    cat("Warning: File not found - ", fantom_tss_file, "\n")
    cat("Skipping dataset: ", sel_dataset, "\n\n")
    next
  }
  
  fantom_tss <- read.table(fantom_tss_file, header = TRUE, sep = "\t", 
                              stringsAsFactors = FALSE, check.names = FALSE)
  cat("Number of FANTOM5 TSS clusters: ", nrow(fantom_tss), "\n")
  
  # Create cluster ID
  fantom_tss <- fantom_tss %>%
      mutate(cluster_id = paste0(seqnames, ":", start, "-", end))
  
  # Create GenomicRanges object for FANTOM TSS
  fantom_tss_gr <- CreateGRangesObject(ranges = fantom_tss,
                                     seqnames = "seqnames",
                                     start = "start",
                                     end = "end",
                                     metadata = fantom_tss["cluster_id"],
                                     retain_md = TRUE)
  
  # Overlap gene TSS with FANTOM5 TSS clusters
  overlaps <- findOverlaps(query = gene_tss_gr,
                           subject = fantom_tss_gr,
                           ignore.strand = TRUE)
  cat("Number of overlaps: ", length(overlaps), "\n")
  
  # Extract metadata from overlaps
  gene_tss_ids <- gene_tss_gr$gene_tss_id[queryHits(overlaps)]
  cluster_ids <- fantom_tss_gr$cluster_id[subjectHits(overlaps)]
  
  # Create mapping dataframe
  tss_to_cluster <- data.frame(
      gene_tss_id = gene_tss_ids,
      cluster_id = cluster_ids,
      stringsAsFactors = FALSE
  )
  
  # Get TPM columns from fantom_tss
  tpm_cols <- grep("tpm", names(fantom_tss), value = TRUE, ignore.case = TRUE)
  
  # Merge expression data with the mapping
  tss_expression <- tss_to_cluster %>%
      left_join(fantom_tss %>% dplyr::select(cluster_id, all_of(tpm_cols)), 
                by = "cluster_id")
  
  # Handle duplicates (average expression values and concatenate cluster_ids)
  tss_expression <- tss_expression %>%
      group_by(gene_tss_id) %>%
      summarise(cluster_id = paste(cluster_id, collapse = ","),
                across(all_of(tpm_cols), ~ mean(.x, na.rm = TRUE)),
                .groups = "drop")
  
  # Transfer to gene_pairs_ranked - lncRNA expression
  lncrna_expr <- tss_expression %>%
      dplyr::select(gene_tss_id, cluster_id, everything())
  names(lncrna_expr)[-1] <- paste0("lncRNA_", names(lncrna_expr)[-1])
  
  gene_pairs_expression <- gene_pairs_ranked %>%
      left_join(lncrna_expr, by = c("lncRNA_tss_id" = "gene_tss_id"))
  
  # Transfer to gene_pairs_ranked - TF expression
  tf_expr <- tss_expression %>%
      dplyr::select(gene_tss_id, cluster_id, everything())
  names(tf_expr)[-1] <- paste0("TF_", names(tf_expr)[-1])
  
  gene_pairs_expression <- gene_pairs_expression %>%
      left_join(tf_expr, by = c("TF_tss_id" = "gene_tss_id"))
  
  # Fill NAs with zeros
  # gene_pairs_expression[is.na(gene_pairs_expression)] <- 0
  
  cat("Number of columns in output: ", ncol(gene_pairs_expression), "\n")
  
  # Save the final table
  output_file <- paste0(out_dir, sel_dataset, "_", date, ".tsv")
  write.table(gene_pairs_expression, file = output_file, sep = "\t", 
              quote = FALSE, row.names = FALSE)
  cat("Saved: ", output_file, "\n\n")
}

cat("========================================\n")
cat("Processing complete!\n")
cat("========================================\n")
