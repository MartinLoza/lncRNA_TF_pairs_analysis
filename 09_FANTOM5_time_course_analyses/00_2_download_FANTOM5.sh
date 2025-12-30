#!/bin/bash
#
# Simple FANTOM5 download script using URL list file
#

# Configuration
URL_LIST_FILE="/mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/Data/FANTOM5/FANTOM5_download_urls_251219.txt"
DOWNLOAD_DIR="/mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/Data/FANTOM5/CTSS_files/"

# Create download directory
echo "Creating download directory: $DOWNLOAD_DIR"
mkdir -p "$DOWNLOAD_DIR"
cd "$DOWNLOAD_DIR"

# Check if URL list exists
if [ ! -f "$URL_LIST_FILE" ]; then
    echo "Error: URL list file not found: $URL_LIST_FILE"
    echo "Please run the R notebook to generate the URL list first."
    exit 1
fi

# Count total URLs
total_urls=$(wc -l < "$URL_LIST_FILE")

echo "========================================="
echo "FANTOM5 Simple Download"
echo "========================================="
echo "URL list file: $URL_LIST_FILE"
echo "Download directory: $DOWNLOAD_DIR"
echo "Total URLs: $total_urls"
echo "========================================="

# Download each URL
counter=0
while IFS= read -r url; do
    ((counter++))
    
    # Extract filename from URL
    filename=$(basename "$url" | sed 's/%2520/ /g' | sed 's/%252c/,/g' | sed 's/%2528/(/g' | sed 's/%2529/)/g')
    
    echo "[$counter/$total_urls] $filename"
    
    if [ -f "$filename" ]; then
        echo "  ✓ Already exists"
    else
        if wget -c -q --timeout=60 "$url" -O "$filename"; then
            if [ -s "$filename" ]; then
                size=$(stat -c%s "$filename" 2>/dev/null || echo "unknown")
                echo "  ✓ Downloaded ($size bytes)"
            else
                echo "  ✗ Empty file"
                rm -f "$filename"
            fi
        else
            echo "  ✗ Failed"
        fi
    fi
    
    # Small delay
    sleep 0.5
done < "$URL_LIST_FILE"

echo "========================================="
echo "Download completed!"
downloaded=$(find . -name "*.ctss.bed.gz" -o -name "*.nobarcode.bam" | wc -l)
echo "Downloaded files: $downloaded / $total_urls"
echo "========================================="