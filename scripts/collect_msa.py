import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os
import time
import warnings
from typing import List, Dict, Tuple

# Suppress warnings
warnings.filterwarnings("ignore")


def validate_file(file_path: str, error_message: str) -> None:
    """Validate if a file exists."""
    if not os.path.isfile(file_path):
        raise FileNotFoundError(error_message)


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    description = "A script to process genomic data."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-g', '--genome_reference', type=str, required=True, help='Path to the genome reference file')
    parser.add_argument('-c', '--contig_fasta', type=str, required=True, help='Path to the contig FASTA file')
    parser.add_argument('-i', '--ictv_report', type=str, required=True, help='Path to the ICTV report file')
    parser.add_argument('-r', '--contig_report', type=str, required=True, help='Path to the contig report file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output folder path')
    parser.add_argument('-t', '--ictv_taxo', type=str, default='ictv_tables/ictv_taxo.xlsx', help='Path to the ICTV taxonomy file')

    return parser.parse_args()


def is_intersecting(segment1: List[int], segment2: List[int]) -> bool:
    """Check if two segments intersect."""
    segment1 = sorted(segment1)
    segment2 = sorted(segment2)
    return not (segment1[1] * 1.1 <= segment2[0] or segment1[0] * 1.1 >= segment2[1])


def remove_intersecting_segments(segments: List[List[int]]) -> Tuple[Dict[int, List[int]], List[int]]:
    """Remove intersecting segments from a list."""
    result = {segments[0][0]: segments[0][1]}
    saved_indices = [segments[0][0]]

    for i in range(1, len(segments)):
        current_index, current_segment = segments[i]
        is_intersects = any(is_intersecting(current_segment, result[saved_index]) for saved_index in saved_indices)

        if not is_intersects:
            result[current_index] = current_segment
            saved_indices.append(current_index)

    return result, saved_indices


def extract_translated_seq(match_row: pd.Series, sequence: SeqRecord) -> str:
    """Extract translated sequence based on match row."""
    From, To, Frame = match_row['From'], match_row['To'], match_row['Frame']
    if From > To:
        From, To = To, From

    shift_dict = {1: 0, 2: 1, 3: 2, 4: 0, 5: 1, 6: 2}
    shift = shift_dict[Frame]

    try:
        if Frame in [1, 2, 3]:
            translated_seq = sequence.seq[From + shift:To].translate()
        else:
            translated_seq = sequence.seq[From:To + 1].reverse_complement()[shift:].translate()
    except Exception as e:
        print(f"Error translating sequence: {e}")
        raise

    return translated_seq


def gb_str_to_list(genbank_str: str) -> List[str]:
    """Convert GENBANK accession string to a list of accessions."""
    if pd.isna(genbank_str):
        return pd.NA

    accessions = genbank_str.split(';')
    accessions = [acc.split(':')[-1].replace(' ', '') for acc in accessions]
    return [acc.split('(')[0] for acc in accessions]


def load_ictv_data(ictv_report_path: str, ictv_taxo_path: str) -> Tuple[pd.DataFrame, Dict]:
    """Load and preprocess ICTV data."""
    ictv_report = pd.read_csv(ictv_report_path, index_col=0)
    ictv_report['Name'] = ictv_report['Name'].apply(lambda x: x.split('.')[0])

    ictv_taxo = pd.read_excel(ictv_taxo_path).dropna(subset='Virus GENBANK accession')
    ictv_taxo['genbank_list'] = ictv_taxo['Virus GENBANK accession'].apply(gb_str_to_list)
    ictv_taxo['Isolate_id'] = ictv_taxo['Species'] + '{' + ictv_taxo['Sort'].astype('str') + '_' + ictv_taxo['Isolate Sort'].astype('str') + '}' + '|' + ictv_taxo['Genus']

    # Create taxo_index
    taxo_index = pd.DataFrame({
        'ictv_taxo_index': [item.Isolate_id for _, item in ictv_taxo.iterrows() for accession in item['genbank_list']],
        'accessions': [accession for _, item in ictv_taxo.iterrows() for accession in item['genbank_list']]
    }).drop_duplicates(subset='accessions').set_index('accessions').to_dict(orient='index')

    return ictv_report, taxo_index


def process_contig(name_c_contig: str, contigs: pd.DataFrame, ictv_report: pd.DataFrame, taxo_index: Dict, genome_reference: str, contigs_fasta: str, output_folder: str) -> None:
    """Process a single contig."""
    print(f'{name_c_contig} ', time.strftime('%H:%M:%S'))

    original_array = contigs.query(f"Name == '{name_c_contig}'").reset_index()[['index', 'From', 'To']].to_numpy()
    reshaped_array = [[i[0], [i[1], i[2]]] for i in original_array]
    filtered_segments, saved_indices = remove_intersecting_segments(reshaped_array)

    contigs_sample = contigs.iloc[saved_indices] if saved_indices else contigs.head(3)
    contig_rows = contigs_sample[['From', 'To', 'Frame', 'Query']].to_dict(orient='index')
    hmms = list(contigs_sample.Query.unique())

    for cur_hmm in hmms[:3]:
        print(f'    {cur_hmm} ', time.strftime('%H:%M:%S'))
        ictv_names = ictv_report.query(f'Query == "{cur_hmm}"')[['Name', 'From', 'To', 'Frame']].reset_index()
        ictv_names2 = ictv_names['Name'].unique().tolist()
        genome_seqs = []

        with open(genome_reference, "r") as file:
            for sequence in SeqIO.parse(file, "fasta"):
                seq_id = sequence.id.split('.')[0]
                if seq_id not in ictv_names2:
                    continue

                ictv_matches = [ictv_match for _, ictv_match in ictv_names.iterrows() if str(ictv_match['Name']) == seq_id]
                for match_count, ictv_match in enumerate(ictv_matches):
                    translated_seq = extract_translated_seq(ictv_match, sequence)
                    output_seq = SeqRecord(translated_seq, id=str(taxo_index[seq_id]['ictv_taxo_index']).replace(' ', '_') + f"#{match_count}")
                    genome_seqs.append(output_seq)

        contig_seqs = []
        with open(contigs_fasta, "r") as file:
            for sequence in SeqIO.parse(file, "fasta"):
                if sequence.id != name_c_contig:
                    continue

                for k, values in contig_rows.items():
                    if values['Query'] != cur_hmm:
                        continue

                    From, To, Frame = values['From'], values['To'], values['Frame']
                    if From > To:
                        From, To = To, From

                    shift_dict = {1: 0, 2: 1, 3: 2, 4: 0, 5: 1, 6: 2}
                    shift = shift_dict[Frame]

                    if Frame in [1, 2, 3]:
                        aa_seq = SeqRecord(sequence.seq[From + shift:To].translate())
                    else:
                        aa_seq = SeqRecord(sequence.seq[From:To].reverse_complement()[shift:].translate())

                    aa_seq.id = f'CONTIG__{sequence.name}__{values["Query"]}__{k}'
                    contig_seqs.append(aa_seq)

        for c_seq in contig_seqs:
            output_name = f"{output_folder}/test_{c_seq.id.replace('CONTIG__', '')}.fasta"
            SeqIO.write([c_seq] + genome_seqs, output_name, 'fasta')


def main():
    start_time = time.process_time()
    args = parse_arguments()

    # Validate input files
    validate_file(args.genome_reference, "Genome reference file doesn't exist")
    validate_file(args.contig_fasta, "Contig FASTA file doesn't exist")
    validate_file(args.ictv_report, "ICTV report file doesn't exist")
    validate_file(args.contig_report, "Contig report file doesn't exist")
    validate_file(args.ictv_taxo, "ICTV taxonomy file doesn't exist")

    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)

    print("Start reading frames ", time.strftime('%H:%M:%S'))

    # Load and preprocess data
    ictv_report, taxo_index = load_ictv_data(args.ictv_report, args.ictv_taxo)

    # Load contigs
    contigs = pd.read_csv(args.contig_report, index_col=0).query('Lengh > 300')
    contigs_names = contigs['Name'].unique()

    # Process each contig
    for name_c_contig in contigs_names:
        process_contig(name_c_contig, contigs, ictv_report, taxo_index, args.genome_reference, args.contig_fasta, args.output)

    print(f"Total processing time: {time.process_time() - start_time} seconds")


if __name__ == "__main__":
    main()
