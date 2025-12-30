#!/usr/bin/env Rscript
# Coexpression analysis for FANTOM5 time course data
# Author: Martin Loza
# Date: 25/12/30
# Description: Perform coexpression analysis for lncRNA-TF pairs across FANTOM5 datasets

# Change R language to English
Sys.setenv(LANGUAGE = "en")

# Init
suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(patchwork)
})

# Local variables 
seed = 777
date = "251230"
# Available datasets (using underscore naming convention)
datasets <- c("Lymphatic", "Aortic", "MCF7", "Monocyte-derived", "Saos-2")

# Define colors for plots
red = "#E41A1C"
blue = "#090a0bff"
green = "#4DAF4A"
purple = "#984EA3"
text_size = 18
width = 18.6
dot_size = 4
line_size = 1.5
dpi = 300

# Directories
main_dir = "/mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/09_FANTOM5_time_course_analyses/Results/"
gene_pairs_dir = paste0(main_dir, "Gene_pairs_FANTOM_expression/")
out_dir = paste0(main_dir, "Coexpression_results/")
plots_dir = paste0(main_dir, "Coexpression_results/Plots/")

# Create output directories if they don't exist
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}
if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
}

# ============================================================================
# Function to perform coexpression analysis for one dataset
# ============================================================================

perform_coexpression_analysis <- function(dataset_name) {
    
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("Processing dataset:", dataset_name, "\n")
    cat(rep("=", 80), "\n\n", sep = "")
    
    # Load the pre-annotated gene pairs with expression data
    gene_pairs_file <- paste0(gene_pairs_dir, dataset_name, "_251230.tsv")
    
    if (!file.exists(gene_pairs_file)) {
        cat("WARNING: Gene pairs file not found for", dataset_name, "\n")
        cat("Expected file:", gene_pairs_file, "\n")
        return(NULL)
    }
    
    gene_pairs_expression <- read.table(gene_pairs_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cat("Loaded", nrow(gene_pairs_expression), "gene pairs with expression data\n")
    
    if (nrow(gene_pairs_expression) == 0) {
        cat("WARNING: No gene pairs with expression data for", dataset_name, "\n")
        return(NULL)
    }
    # print(head(gene_pairs_expression,3))
    
    # Extract expression matrices for lncRNA and TF genes
    lncRNA_expression <- gene_pairs_expression %>%
        select(starts_with("lncRNA_tpm_")) %>%
        as.matrix()
    
    TF_expression <- gene_pairs_expression %>%
        select(starts_with("TF_tpm_")) %>%
        as.matrix()
    
    # Set row names to gene pair IDs
    rownames(lncRNA_expression) <- gene_pairs_expression$gene_pair_id
    rownames(TF_expression) <- gene_pairs_expression$gene_pair_id
    
    cat("Expression matrices dimensions:", dim(lncRNA_expression), "\n")
    
    # Calculate correlations
    n_pairs <- nrow(gene_pairs_expression)
    coexpression_results <- data.frame(
        gene_pair_id = gene_pairs_expression$gene_pair_id,
        lncRNA_gene_id = gene_pairs_expression$lncRNA_tss_id,
        TF_gene_id = gene_pairs_expression$TF_tss_id,
        pearson_correlation = numeric(n_pairs),
        pearson_pvalue = numeric(n_pairs),
        spearman_correlation = numeric(n_pairs),
        spearman_pvalue = numeric(n_pairs),
        n_samples = integer(n_pairs)
    )
    
    cat("Calculating correlations...\n")
    for (i in 1:n_pairs) {
        lncRNA_expr <- as.numeric(lncRNA_expression[i, ])
        TF_expr <- as.numeric(TF_expression[i, ])
        
        # Check if either gene has all NA values
        if (all(is.na(lncRNA_expr)) || all(is.na(TF_expr))) {
            # Set correlation and p-values to NA
            coexpression_results$pearson_correlation[i] <- NA
            coexpression_results$pearson_pvalue[i] <- NA
            coexpression_results$spearman_correlation[i] <- NA
            coexpression_results$spearman_pvalue[i] <- NA
            coexpression_results$n_samples[i] <- 0
            next
        }

        # # Check if either gene has all Zero values
        # if (sum(is.na(lncRNA_expr)) || all(is.na(TF_expr))) {
        #     # Set correlation and p-values to NA
        #     coexpression_results$pearson_correlation[i] <- NA
        #     coexpression_results$pearson_pvalue[i] <- NA
        #     coexpression_results$spearman_correlation[i] <- NA
        #     coexpression_results$spearman_pvalue[i] <- NA
        #     coexpression_results$n_samples[i] <- 0
        #     next
        # }
        
        # Count valid samples (both genes have non-NA values)
        valid_samples <- !is.na(lncRNA_expr) & !is.na(TF_expr)
        n_valid <- sum(valid_samples)
        
        # Check if we have enough valid samples for correlation (minimum 3)
        if (n_valid < 3) {
            coexpression_results$pearson_correlation[i] <- NA
            coexpression_results$pearson_pvalue[i] <- NA
            coexpression_results$spearman_correlation[i] <- NA
            coexpression_results$spearman_pvalue[i] <- NA
            coexpression_results$n_samples[i] <- n_valid
            next
        }
        
        # Pearson correlation
        pearson_test <- cor.test(lncRNA_expr, TF_expr, method = "pearson", use = "complete.obs")
        coexpression_results$pearson_correlation[i] <- pearson_test$estimate
        coexpression_results$pearson_pvalue[i] <- pearson_test$p.value
        
        # Spearman correlation
        spearman_test <- cor.test(lncRNA_expr, TF_expr, method = "spearman", exact = FALSE, use = "complete.obs")
        coexpression_results$spearman_correlation[i] <- spearman_test$estimate
        coexpression_results$spearman_pvalue[i] <- spearman_test$p.value
        
        coexpression_results$n_samples[i] <- n_valid
    }
    
    # Apply FDR correction (only on non-NA p-values)
    coexpression_results$pearson_fdr <- p.adjust(coexpression_results$pearson_pvalue, method = "fdr")
    coexpression_results$spearman_fdr <- p.adjust(coexpression_results$spearman_pvalue, method = "fdr")
    
    # Add significance flags (handle NA values)
    coexpression_results$pearson_significant_0.05 <- ifelse(is.na(coexpression_results$pearson_fdr), 
                                                              FALSE, 
                                                              coexpression_results$pearson_fdr < 0.05)
    coexpression_results$pearson_significant_0.01 <- ifelse(is.na(coexpression_results$pearson_fdr), 
                                                              FALSE, 
                                                              coexpression_results$pearson_fdr < 0.01)
    coexpression_results$spearman_significant_0.05 <- ifelse(is.na(coexpression_results$spearman_fdr), 
                                                               FALSE, 
                                                               coexpression_results$spearman_fdr < 0.05)
    coexpression_results$spearman_significant_0.01 <- ifelse(is.na(coexpression_results$spearman_fdr), 
                                                               FALSE, 
                                                               coexpression_results$spearman_fdr < 0.01)
    
    # Print summary statistics
    cat("\n=== Coexpression Analysis Summary ===\n")
    cat("Total gene pairs analyzed:", nrow(coexpression_results), "\n")
    cat("Gene pairs with valid correlations:", sum(!is.na(coexpression_results$pearson_correlation)), "\n")
    cat("Gene pairs with NA correlations:", sum(is.na(coexpression_results$pearson_correlation)), "\n\n")
    
    cat("--- Pearson Correlation ---\n")
    cat("Mean correlation:", round(mean(coexpression_results$pearson_correlation, na.rm = TRUE), 3), "\n")
    cat("Median correlation:", round(median(coexpression_results$pearson_correlation, na.rm = TRUE), 3), "\n")
    cat("Significant pairs (FDR < 0.05):", sum(coexpression_results$pearson_significant_0.05, na.rm = TRUE), "\n")
    cat("Significant pairs (FDR < 0.01):", sum(coexpression_results$pearson_significant_0.01, na.rm = TRUE), "\n")
    cat("Positively correlated (r > 0):", sum(coexpression_results$pearson_correlation > 0, na.rm = TRUE), "\n")
    cat("Negatively correlated (r < 0):", sum(coexpression_results$pearson_correlation < 0, na.rm = TRUE), "\n\n")
    
    cat("--- Spearman Correlation ---\n")
    cat("Mean correlation:", round(mean(coexpression_results$spearman_correlation, na.rm = TRUE), 3), "\n")
    cat("Median correlation:", round(median(coexpression_results$spearman_correlation, na.rm = TRUE), 3), "\n")
    cat("Significant pairs (FDR < 0.05):", sum(coexpression_results$spearman_significant_0.05, na.rm = TRUE), "\n")
    cat("Significant pairs (FDR < 0.01):", sum(coexpression_results$spearman_significant_0.01, na.rm = TRUE), "\n")
    cat("Positively correlated (rho > 0):", sum(coexpression_results$spearman_correlation > 0, na.rm = TRUE), "\n")
    cat("Negatively correlated (rho < 0):", sum(coexpression_results$spearman_correlation < 0, na.rm = TRUE), "\n")
    
    # Save results (keep NA values in the output)
    sel_cols_save <- c('gene_pair_id', 'lncRNA_gene_id', 'TF_gene_id', 'n_samples',
                       'pearson_correlation', 'pearson_pvalue', 'pearson_fdr',
                       'spearman_correlation', 'spearman_pvalue', 'spearman_fdr')
    tmp_data_save <- coexpression_results %>%
        select(all_of(sel_cols_save))
    
    output_file <- paste0(out_dir, dataset_name, "_coexpression_results_", date, ".tsv")
    write.table(tmp_data_save, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\nResults saved to:", output_file, "\n")
    
    return(coexpression_results)
}

# ============================================================================
# Function to create plots for one dataset
# ============================================================================

create_plots <- function(dataset_name, coexpression_results) {
    
    if (is.null(coexpression_results)) {
        return()
    }
    
    cat("\nCreating plots for", dataset_name, "...\n")
    
    # Create dataset-specific plot directory
    dataset_plot_dir <- paste0(plots_dir, dataset_name, "/")
    if (!dir.exists(dataset_plot_dir)) {
        dir.create(dataset_plot_dir, recursive = TRUE)
    }
    
    # 1. Histogram plots
    p1 <- ggplot(coexpression_results, aes(x = pearson_correlation)) +
        geom_histogram(bins = 50, fill = "white", color = "black") +
        geom_vline(xintercept = 0, color = red, linetype = "dashed", linewidth = 1) +
        labs(x = "Pearson Correlation", y = "Count") +
        theme_minimal() +
        theme(text = element_text(size = text_size),
              strip.text = element_text(size = text_size + 1),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
    
    p2 <- ggplot(coexpression_results, aes(x = spearman_correlation)) +
        geom_histogram(bins = 50, fill = "white", color = "black") +
        geom_vline(xintercept = 0, color = red, linetype = "dashed", linewidth = 1) +
        labs(x = "Spearman Correlation", y = "Count") +
        theme_minimal() +
        theme(text = element_text(size = text_size),
              strip.text = element_text(size = text_size + 1),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
    
    p_histogram <- p1 + p2 + plot_annotation(
        title = gsub("_", "-", dataset_name),
        theme = theme(plot.title = element_text(hjust = 0.5, size = text_size))
    )
    
    ggsave(filename = paste0(dataset_plot_dir, dataset_name, "_histogram_", date, ".pdf"),
           plot = p_histogram, width = width * 0.6, height = 5.0, units = "in", dpi = dpi)
    
    # 2. Volcano plots
    plot_data <- coexpression_results %>%
        mutate(
            pearson_log_fdr = -log10(pearson_fdr),
            spearman_log_fdr = -log10(spearman_fdr),
            pearson_color = ifelse(pearson_significant_0.05, "Significant", "Not significant"),
            spearman_color = ifelse(spearman_significant_0.05, "Significant", "Not significant")
        ) %>%
        filter(!is.na(pearson_log_fdr) & !is.na(spearman_log_fdr))
    
    color_palette <- c("Significant" = red, "Not significant" = "grey")
    
    p1 <- ggplot(plot_data, aes(x = pearson_correlation, y = pearson_log_fdr, color = pearson_color)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        scale_color_manual(values = color_palette, name = "Significance") +
        labs(x = "Pearson Correlation", y = "-log10(FDR)") +
        theme_minimal() +
        theme(text = element_text(size = text_size),
              strip.text = element_text(size = text_size + 1),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
    
    p2 <- ggplot(plot_data, aes(x = spearman_correlation, y = spearman_log_fdr, color = spearman_color)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        scale_color_manual(values = color_palette, name = "Significance") +
        labs(x = "Spearman Correlation", y = "-log10(FDR)") +
        theme_minimal() +
        theme(text = element_text(size = text_size),
              strip.text = element_text(size = text_size + 1),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
    
    p_volcano <- p1 + p2 + plot_layout(guides = "collect") +
        plot_annotation(
            title = gsub("_", "-", dataset_name),
            theme = theme(plot.title = element_text(hjust = 0.5, size = text_size))
        ) & theme(legend.position = "bottom")
    
    ggsave(filename = paste0(dataset_plot_dir, dataset_name, "_volcano_", date, ".pdf"),
           plot = p_volcano, width = width * 0.6, height = 5.0, units = "in", dpi = dpi)
    
    # 3. Method comparison plot
    method_cor_test <- cor.test(coexpression_results$pearson_correlation,
                                coexpression_results$spearman_correlation,
                                use = "complete.obs")
    method_cor <- method_cor_test$estimate
    method_pval <- method_cor_test$p.value
    
    p_comparison_methods <- ggplot(coexpression_results, aes(x = pearson_correlation, y = spearman_correlation)) +
        geom_point(alpha = 0.3, size = 1.5) +
        geom_abline(intercept = 0, slope = 1, color = red, linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
        geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
        annotate("text", x = -Inf, y = Inf,
                 label = paste0("Correlation: ", round(method_cor, 3), "\nP-value: ", format(method_pval, scientific = TRUE, digits = 3)),
                 hjust = -0.1, vjust = 1.5, size = 5) +
        labs(x = "Pearson Correlation", y = "Spearman Correlation",
             title = gsub("_", "-", dataset_name)) +
        theme_minimal() +
        theme(text = element_text(size = text_size),
              strip.text = element_text(size = text_size + 1),
              axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5, size = text_size))
    
    ggsave(filename = paste0(dataset_plot_dir, dataset_name, "_comparison_methods_", date, ".pdf"),
           plot = p_comparison_methods, width = width * 0.3, height = 5.0, units = "in", dpi = dpi)
    
    cat("Plots saved to:", dataset_plot_dir, "\n")
}

# ============================================================================
# Process all datasets
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("Starting coexpression analysis for all FANTOM5 datasets\n")
cat(rep("=", 80), "\n", sep = "")

results_list <- list()

for (dataset in datasets) {
    results <- perform_coexpression_analysis(dataset)
    if (!is.null(results)) {
        results_list[[dataset]] <- results
        create_plots(dataset, results)
    }
}

cat("\n", rep("=", 80), "\n", sep = "")
cat("Analysis complete for all datasets!\n")
cat("Results saved to:", out_dir, "\n")
cat("Plots saved to:", plots_dir, "\n")
cat(rep("=", 80), "\n", sep = "")