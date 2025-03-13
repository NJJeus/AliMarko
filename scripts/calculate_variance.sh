#!/bin/bash

# Default threshold value
THRESHOLD=200

# Function to display usage
usage() {
    echo "Usage: $0 -b|--bamfile <bamfile> -l|--listfile <covfile> -f|--reference <reference> [-t|--threshold <threshold>] -o|--output_dir <output>"
    exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b=*|--bamfile=*)
            BAMFILE="${1#*=}"
            ;;
        -l=*|--listfile=*)
            COVFILE="${1#*=}"
            ;;
        -f=*|--reference=*)
            REFERENCE="${1#*=}"
            ;;
        -t=*|--threshold=*)
            THRESHOLD="${1#*=}"
            ;;
        -o=*|--output_dir=*)
            OUT="${1#*=}"
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$BAMFILE" || -z "$COVFILE" || -z "$REFERENCE" || -z "$OUT" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Extract base name of the BAM file
BAMNAME=$(basename "$BAMFILE" .bam)

# Initialize output file
> "$OUT"
echo "rname,snps,deep_sites" >> "$OUT"

# Process coverage file
tail -n +2 "$COVFILE" | awk -v threshold="$THRESHOLD" '$6 > threshold' | cut -f1 | while read -r line; do
    echo "Processing region: $line"
    
    # Count variants
    n_var=$(freebayes -f "$REFERENCE" "$BAMFILE" -p 1 -r "$line" \
        | bcftools filter -e "QUAL <= 20" | grep -vc '^#')
    
    # Calculate depth
    depth=$(samtools depth "$BAMFILE" -r "$line" | awk '$3 > 5' | wc -l)
    
    # Append results to output file
    echo "$line,$n_var,$depth" >> "$OUT"
done

echo "Processing complete. Results saved to $OUT."
