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

# Iterate through each file in the input folder
for file in "$input_folder"/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        output_file="$output_folder/$filename"
        temp_file="$output_folder/${filename}.temp.fasta"

        # Perform some operation on the input file and write the result to the output file
        echo "Processing: $filename"

        # Step 1: Run the Python script to select MSA
        if ! python scripts/select_msa.py -i "$file" -o "$temp_file" -n 25; then
            echo "Error: Failed to process $filename with select_msa.py"
            continue
        fi

        # Step 2: Run MAFFT on the temporary file
        if ! mafft --thread "$threads" "$temp_file" > "$output_file"; then
            echo "Error: Failed to process $filename with MAFFT"
            continue
        fi

        # Clean up temporary file
        rm "$temp_file"

        echo "Completed: $filename"
    fi
done

echo "All files processed."
