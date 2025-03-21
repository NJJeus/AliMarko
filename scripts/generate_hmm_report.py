import pandas as pd
import numpy as np
import argparse
import os
from Bio import SeqIO
from typing import List, Optional, Dict
import logging


# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Script to generate reports from the result file of the script 'analyze_seq.py'.")

    parser.add_argument('-i', '--input', type=str, required=True, help='Input file or regular expression of file')
    parser.add_argument('-i2', '--input2', type=str, help='Second input file or regular expression of file')
    parser.add_argument('-c', '--contigs_file', type=str, required=True, help='Contigs file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output folder')
    parser.add_argument('-m', '--hmm_info', type=str, required=True, help='Folder with HMM drawings')
    parser.add_argument('-t', '--hmm_out', type=str, help='One column CSV with models that have passed the threshold')
    parser.add_argument('-n', '--output_contigs', type=str, required=True, help='Output contigs that were matched by HMM')

    return parser.parse_args()


def check_file_exists(file_path: str, message: str) -> None:
    """Check if a file exists, otherwise log an error message and exit."""
    if not os.path.isfile(file_path):
        logging.error(message)
        exit(1)


def load_csv_file(file_path: str, index_col: int = 0) -> pd.DataFrame:
    """Load a CSV file into a DataFrame."""
    return pd.read_csv(file_path, index_col=index_col)


def export_seq_data(seq: str) -> List[str]:
    """Extract and format sequence data."""
    try:
        name, chain, part, length, length_contig, skew, _ = seq.split(';')
        chain = chain.replace('chain_', '')
        part = part.replace('part_', '')
        length = length.replace('len_', '')
        skew = skew.replace('skew_', '')
        length_contig = length_contig.replace('len_contig_', '')
        frame = int(chain.replace('+1', '1').replace('-1', '4')) + int(skew)
        return [name, part, length, frame, length_contig]
    except ValueError as e:
        logging.error(f"Error parsing sequence data: {e}")
        return [None] * 5  # Return a list of None values to maintain DataFrame structure


def get_coordinates(dataframe: pd.DataFrame) -> np.ndarray:
    """Calculate coordinates based on frame and part using vectorized operations."""
    frame_condition = dataframe['Frame'] <= 3
    from_coords = np.where(
        frame_condition,
        dataframe['From'] * 3 + dataframe['Part'],
        (dataframe['Lengh'] - dataframe['From'] * 3) + dataframe['Part']
    )
    to_coords = np.where(
        frame_condition,
        dataframe['To'] * 3 + dataframe['Part'],
        (dataframe['Lengh'] - dataframe['To'] * 3) + dataframe['Part']
    )
    return np.column_stack((from_coords, to_coords))


def create_report(input_df: pd.DataFrame, hmm_df: pd.DataFrame) -> pd.DataFrame:
    """Create a report by merging input data with HMM info."""
    report = input_df.merge(hmm_df, left_on='Query', right_on='Model_ID', how='left').drop('Model_ID', axis=1)
    report['Score_ratio'] = report.Score / report.Threshold
    report = report.sort_values('Score_ratio', ascending=False)

    try:
        seq_data = np.array(report.Name.apply(export_seq_data).to_list())
        report[['Name', 'Part', 'Lengh', 'Frame', 'Length_contig']] = seq_data
        report[['Lengh', 'From', 'To', 'Part', 'Frame', 'Length_contig']] = report[
            ['Lengh', 'From', 'To', 'Part', 'Frame', 'Length_contig']
        ].astype('int')
        report[['From', 'To']] = get_coordinates(report)
        report = report.drop('Part', axis=1)
    except Exception as e:
        logging.error(f"An error occurred while processing the report: {e}")

    return report


def main() -> None:
    """Main function to execute the script."""
    args = parse_arguments()

    # Validate input files
    check_file_exists(args.input, "An input file does not exist")
    check_file_exists(args.contigs_file, "An input contigs file does not exist")
    check_file_exists(args.hmm_info, "HMM info folder does not exist")

    # Load data
    input_file = load_csv_file(args.input)
    hmm_info = load_csv_file(args.hmm_info)
    contigs_file = args.contigs_file

    # Process input files
    report = create_report(input_file, hmm_info)

    if args.input2:
        check_file_exists(args.input2, "An input file 2 does not exist")
        input_file2 = load_csv_file(args.input2)
        report2 = create_report(input_file2, hmm_info)
        report = pd.concat([report, report2])

    # Filter and save reports
    report = report.query('Score_ratio > 0.5')

    if args.hmm_out:
        report[['Query']].drop_duplicates().to_csv(args.hmm_out, header=None, index=None)

    report.reset_index(drop=True).to_csv(args.output)

    # Extract and save selected contigs
    selected_contigs = report['Name'].unique().tolist()
    contigs = []
    with open(contigs_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in selected_contigs:
                contigs.append(record)

    SeqIO.write(contigs, args.output_contigs, 'fasta')
    logging.info("Script executed successfully.")


if __name__ == "__main__":
    main()

