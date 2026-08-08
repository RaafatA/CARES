#!/bin/bash

# Exit if any command fails or an unset variable is used
set -euo pipefail

# 1. Define your paths here
TARGET_DIR="$HOME/projects/tomato/Acc_list"
CSV_FILE="Study_Layout.csv"

# Check if the target directory exists
if [[ ! -d "$TARGET_DIR" ]]; then
    echo "Error: Directory $TARGET_DIR does not exist!"
    exit 1
fi

# Check if the CSV layout file exists
if [[ ! -f "$CSV_FILE" ]]; then
    echo "Error: $CSV_FILE not found!"
    exit 1
fi

echo "Starting batch processing using layout lookup..."
echo "Reading files from: $TARGET_DIR"
echo "----------------------------------------"

# 2. Loop through all .txt files in the specific directory
for file in "$TARGET_DIR"/*.txt; do
    
    # Safety check: if no .txt files exist, exit the loop
    [[ -e "$file" ]] || continue
    
    # Extract just the filename from the path (e.g., GSE109672.txt)
    filename=$(basename "$file")
    
    # Extract the study name by removing the ".txt" extension (e.g., GSE109672)
    study="${filename%.txt}"
    
    # Look up the layout from the CSV file
    layout=$(awk -F',' -v s="$study" '$1==s {print $2}' "$CSV_FILE" | tr -d '\r' | tr '[:lower:]' '[:upper:]')
    
    if [[ -z "$layout" ]]; then
        echo "  [WARNING] No layout found for $study in $CSV_FILE. Skipping..."
        echo "----------------------------------------"
        continue
    fi
    
    echo "Processing Study: $study | Layout: $layout"
    
    # 3. Create the study folder inside the target directory
    mkdir -p "$TARGET_DIR/$study"
    
    # Loop through each sample accession inside the current .txt file
    while read -r sample || [[ -n "$sample" ]]; do
        
        # Skip empty lines or headers
        [[ -z "$sample" || "$sample" == \#* ]] && continue
        
        # Download based on the library layout using fasterq-dump
        # Files will be saved directly into the newly created study folder
        if [[ "$layout" == *"PAIRED"* ]]; then
            echo "    -> Downloading $sample (Paired-End) into $study/"
            fasterq-dump --split-files -O "$TARGET_DIR/$study" "$sample"
            
        elif [[ "$layout" == *"SINGLE"* ]]; then
            echo "    -> Downloading $sample (Single-End) into $study/"
            fasterq-dump -O "$TARGET_DIR/$study" "$sample"
            
        else
            echo "    -> [ERROR] Unknown layout ($layout) for $sample. Skipping."
        fi
        
    done < "$file"
    
    echo "Completed processing for $study"
    echo "----------------------------------------"
    
done

echo "All studies processed successfully!"
