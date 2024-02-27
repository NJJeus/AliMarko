input_folder=$1
output_folder=$2
threads=$3

# Check if the input and output folders are provided
if [ -z "$input_folder" ] || [ -z "$output_folder" ]; then
    echo "Usage: $0 <input_folder> <output_folder>"
    exit 1
fi

# Check if the input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Input folder not found"
    exit 1
fi

# Check if the output folder exists, if not create it
if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi

# Iterate through each file in the input folder
for file in "$input_folder"/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        output_file="$output_folder/$filename"  # Path of the output file

        # Perform some operation on the input file and write the result to the output file
        # For demonstration purpose, let's just copy the input file to the output folder
        python scripts/select_msa.py -i $file -o $file.temp.fasta
        mafft $file.temp.fasta  > $output_file
 
        echo "Processed: $filename"
    fi
done
