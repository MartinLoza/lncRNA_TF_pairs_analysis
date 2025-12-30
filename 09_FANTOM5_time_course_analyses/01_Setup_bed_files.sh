#!/bin/bash
#
# Process FANTOM5 downloaded files: decompress and rename
#

# Configuration
DOWNLOAD_DIR="/mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/Data/FANTOM5/CTSS_files/"
PROCESSED_DIR="/mnt/c/Users/Marti/Documents/Projects/lncRNA_TF_pairs_analysis/Data/FANTOM5/bed_files/"

# Create processed directory
echo "Creating processed directory: $PROCESSED_DIR"
mkdir -p "$PROCESSED_DIR"

# Check if download directory exists
if [ ! -d "$DOWNLOAD_DIR" ]; then
    echo "Error: Download directory not found: $DOWNLOAD_DIR"
    exit 1
fi

cd "$DOWNLOAD_DIR"

# Count total files
total_files=$(find . -name "*.ctss.bed.gz" | wc -l)

echo "========================================="
echo "FANTOM5 File Processing"
echo "========================================="
echo "Source directory: $DOWNLOAD_DIR"
echo "Output directory: $PROCESSED_DIR"
echo "Total .gz files: $total_files"
echo "========================================="

# Process each .gz file
counter=0
for gzfile in *.ctss.bed.gz; do
    [ -f "$gzfile" ] || continue
    
    ((counter++))
    
    echo "[$counter/$total_files] Processing: $gzfile"
    
    # Extract first word and time point
    first_word=$(echo "$gzfile" | awk '{print $1}')
    
    # Extract time point - try multiple patterns
    # Patterns: 00hr30min, 00hr, day1, day2, etc.
    timepoint=$(echo "$gzfile" | grep -oP '\d{2}hr\d{2}min|\d{2}hr|day\d+' | head -1)
    
    if [ -n "$timepoint" ]; then
        new_name="${first_word}_${timepoint}.bed"
        
        # Check if already processed
        if [ -f "$PROCESSED_DIR/$new_name" ]; then
            echo "  ✓ Already processed: $new_name"
        else
            echo "  → Decompressing to: $new_name"
            
            # Decompress and reorganize columns (chr, start, strand, score)
            gunzip -c "$gzfile" | awk 'BEGIN{OFS="\t"} {print $1, $2, $6, $5}' > "$PROCESSED_DIR/$new_name"
            
            if [ -s "$PROCESSED_DIR/$new_name" ]; then
                size=$(stat -c%s "$PROCESSED_DIR/$new_name" 2>/dev/null || echo "unknown")
                echo "  ✓ Success ($size bytes)"
            else
                echo "  ✗ Failed"
                rm -f "$PROCESSED_DIR/$new_name"
            fi
        fi
    else
        echo "  ✗ No timepoint found in filename"
    fi
done

echo "========================================="
echo "Processing completed!"
processed=$(find "$PROCESSED_DIR" -name "*.bed" | wc -l)
echo "Processed BED files: $processed"
echo "Location: $PROCESSED_DIR"
echo "========================================="