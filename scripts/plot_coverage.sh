#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -b|--bamfile <BAMFILE> -l|--listfile <LISTFILE> -f|--reference <REFERENCE> -o|--output_dir <OUTPUT_DIR> [--default]"
    exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b=*|--bamfile=*)
            BAMFILE="${1#*=}"
            ;;
        -l=*|--listfile=*)
            LISTFILE="${1#*=}"
            ;;
        -f=*|--reference=*)
            REFERENCE="${1#*=}"
            ;;
        -o=*|--output_dir=*)
            OUTPUT_DIR="${1#*=}"
            ;;
        --default)
            DEFAULT=YES
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$BAMFILE" || -z "$LISTFILE" || -z "$REFERENCE" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Extract base name of the BAM file
BAMNAME=$(basename "$BAMFILE" .bam)

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each line in the list file
while IFS=',' read -r species rname endpos; do
    # Create species-specific directory
    species_dir="${OUTPUT_DIR}/${species}"
    mkdir -p "$species_dir"

    # Run bamsnap command
    timeout 20 bamsnap \
        -draw coordinates bamplot coverage base \
        -bam "$BAMFILE" \
        -out "${species_dir}/${rname}.png" \
        -pos "${rname}:552-${endpos}" \
        -ref "$REFERENCE" \
        -width 2000 -height 1000 \
        -read_thickness 10 -read_gap_height 5 \
        -title_fontsize 25 \
        -bamplot coverage read \
        -draw coordinates bamplot -no_target_line \
        -title "${BAMNAME}_${species}_${rname}"
done < "$LISTFILE"
