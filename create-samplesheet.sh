#!/bin/bash

# generate_samplesheet.sh
# Usage: ./generate_samplesheet.sh <directory_path> [output_file]

set -euo pipefail

# Default values
DIRECTORY="${1:-dev_results/ORFQuant/}"
OUTPUT_FILE="${2:-samplesheet.csv}"

# Check if directory exists
if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist" >&2
    exit 1
fi

# Function to extract run name from filename
extract_run_name() {
    local filename="$1"
    # Remove the ncORFs_ prefix and file extension
    local run_name=$(echo "$filename" | sed 's/^ncORFs_//' | sed 's/\.[^.]*$//')
    echo "$run_name"
}

# Function to determine file type from extension
get_file_type() {
    local filename="$1"
    local extension="${filename##*.}"
    
    case "$extension" in
        bed) echo "bed" ;;
        bedgraph|bg) echo "bedgraph" ;;
        bigbed|bb) echo "bigbed" ;;
        bigwig|bw) echo "bigwig" ;;
        *) echo "unknown" ;;
    esac
}

# Create the header
echo "Run,FileType,Path" > "$OUTPUT_FILE"

# Process files in the directory
find "$DIRECTORY" -type f \( -name "*.bed" -o -name "*.bedgraph" -o -name "*.bg" -o -name "*.bigbed" -o -name "*.bb" -o -name "*.bigwig" -o -name "*.bw" \) | sort | while read -r filepath; do
    filename=$(basename "$filepath")
    run_name=$(extract_run_name "$filename")
    file_type=$(get_file_type "$filename")
    
    # Convert relative path to absolute path
    abs_path=$(realpath "$filepath")
    
    # Skip if file type is unknown
    if [[ "$file_type" != "unknown" ]]; then
        echo "$run_name,$file_type,$abs_path" >> "$OUTPUT_FILE"
    else
        echo "Warning: Skipping file with unknown type: $filename" >&2
    fi
done

echo "Samplesheet generated: $OUTPUT_FILE"
echo "Found $(tail -n +2 "$OUTPUT_FILE" | wc -l) files"

# Show a preview of the generated samplesheet
echo ""
echo "Preview of generated samplesheet:"
echo "=================================="
head -10 "$OUTPUT_FILE"
if [[ $(wc -l < "$OUTPUT_FILE") -gt 10 ]]; then
    echo "... (showing first 10 lines)"
fi