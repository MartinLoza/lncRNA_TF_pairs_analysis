#!/bin/bash

# This script will convert BED files to FASTA and find enriched motifs using a control
# Since we are using a docker version of meme, we need to make sure we are in the right directory. 
# The Docker's container should be initialized in the following directory:
# /mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis

input_dir="04_Figure_window_size_and_GO/06_Species_comparison_lncRNA_TF/Results/"
output_dir="04_Figure_window_size_and_GO/06_Species_comparison_lncRNA_TF/Results/AME/Human_specific_vs_Conserved_TSS/"

motif_file="Data/JASPAR/JASPAR2022_CORE_vertebrates_redundant_v2.meme"
genome_file="Data/UCSC_hg38/hg38.fa"
# query_sample="Conserved_TSS_regions_260114"
# control_sample="Human_specific_TSS_regions_260114"
query_sample="Human_specific_TSS_regions_260114"
control_sample="Conserved_TSS_regions_260114"

# Make output directory if it doesn't exist
mkdir -p "$output_dir"

# Convert BED files to FASTA files
echo "Converting BED files to FASTA..."

# Convert query sample
if [ -f "${input_dir}${query_sample}.bed" ]; then
    bed2fasta -s -o "${output_dir}${query_sample}.fa" "${input_dir}${query_sample}.bed" "${genome_file}"
    echo "Converted ${query_sample}.bed to ${query_sample}.fa"
else
    echo "File ${input_dir}${query_sample}.bed does not exist"
    exit 1
fi

# Convert control sample
if [ -f "${input_dir}${control_sample}.bed" ]; then
    bed2fasta -s -o "${output_dir}${control_sample}.fa" "${input_dir}${control_sample}.bed" "${genome_file}"
    echo "Converted ${control_sample}.bed to ${control_sample}.fa"
else
    echo "File ${input_dir}${control_sample}.bed does not exist"
    exit 1
fi

# Run AME motif enrichment analysis
echo "Running AME motif enrichment analysis..."
ame --verbose 2 --oc "$output_dir" \
    --scoring avg \
    --method fisher \
    --hit-lo-fraction 0.25 \
    --evalue-report-threshold 10.0 \
    --control "${output_dir}${control_sample}.fa" \
    "${output_dir}${query_sample}.fa" \
    "${motif_file}"

echo "AME analysis complete. Results in ${output_dir}"
