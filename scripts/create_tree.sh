#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_folder> <output_folder> <threads>"
    exit 1
fi

# Assign input arguments to variables
input_folder="$1"
output_folder="$2"
threads="$3"

# Check if the input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Error: Input folder not found: $input_folder"
    exit 1
fi

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Iterate through each FASTA file in the input folder
for file in "$input_folder"/*.fasta; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        filename_no_ext="${filename%.*}"
        output_treefile="$output_folder/${filename_no_ext}.treefile"

        echo "Processing: $filename"

        # Perform phylogenetic tree construction using FastTree
        if ! FastTree "$file" > "$output_treefile"; then
            echo "Error: Failed to construct phylogenetic tree for $filename"
            continue
        fi

        echo "Phylogenetic tree constructed for $filename"
    fi
done

echo "All files processed."
