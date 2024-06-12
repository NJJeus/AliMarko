import argparse
import os
import pysam

def main():
    parser = argparse.ArgumentParser(description='Filter BAM file based on match ratio')
    parser.add_argument('-i', '--input', help='Input BAM file path', required=True)
    parser.add_argument('-o', '--output', help='Output BAM file path', required=True)
    args = parser.parse_args()

    input_bam = args.input
    output_bam = args.output

    if not os.path.exists(input_bam):
        print(f"Error: Input file '{input_bam}' does not exist.")
        return

    with pysam.AlignmentFile(input_bam, "rb") as bamfile, pysam.AlignmentFile(output_bam, "wb", template=bamfile) as output:
        for read in bamfile:
            cigar = read.cigarstring
            if cigar:
                matches = sum([int(length) for op, length in read.cigartuples if op == 0])  # Calculate total matches
                if matches >= (2/3) * len(read.query_sequence):  # Check if matches are at least 2/3 of read length
                    output.write(read)

if __name__ == "__main__":
    main()
