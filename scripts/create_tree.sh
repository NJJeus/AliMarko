# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input_folder output_folder"
    exit 1
fi

input_folder=$1
output_folder=$2
threads=$3

# Check if the input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Input folder not found."
    exit 1
fi

# Check if the output folder exists, if not, create it
if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi

# Iterate through each file in the input folder
for file in "$input_folder"/*.fasta; do
    filename=$(basename "$file")
    filename_no_ext="${filename%.*}"

    # Perform multiple sequence alignment using MAFFT

    # Construct a phylogenetic tree using IQ-TREE
    #iqtree -s "$file" -fast -nt $threads -pre "$output_folder/${filename_no_ext}_tree" -bb 1000
    FastTree "$file" > "$output_folder"/$filename_no_ext.treefile 

    echo "Phylogenetic tree constructed for $filename"
done

echo "All files processed."