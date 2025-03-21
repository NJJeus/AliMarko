import os
import sys
import glob
import re
import gzip
import argparse
from typing import List, Dict, Tuple, Optional, Set
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import pyhmmer

# Constants
DEFAULT_THREADS = 1
DEFAULT_BATCH_SIZE = 1000
DEFAULT_THRESHOLD = 50
DEFAULT_TO_STOP = False

# Genetic code table
GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    description = "The script gets a fastq or a fasta file, translates it, and analyzes sequences with HMM profiles."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-f', '--file', type=str, required=True, help='A file or directory with fastq/fasta extension')
    parser.add_argument('-o', '--output', type=str, required=True, help='An output folder')
    parser.add_argument('-m', '--hmm', type=str, required=True, help='A file or directory with HMM profile')
    parser.add_argument('-t', '--threads', type=int, default=DEFAULT_THREADS, help='Amount of threads to use')
    parser.add_argument('-b', '--batch', type=int, default=DEFAULT_BATCH_SIZE, help='Batch size for sequence analysis')
    parser.add_argument('-s', '--to_stop', action='store_true', help='Translate sequences only to stop codon')
    parser.add_argument('-r', '--threshold', type=int, default=DEFAULT_THRESHOLD, help='Minimum amino acid sequence length for analysis')
    parser.add_argument('-l', '--positive_list', type=str, help='A CSV file with HMMs to use. It should have only one column')

    return parser.parse_args()

def validate_input_path(path: str, is_directory: bool = False) -> str:
    """Validate and ensure the existence of the input path (file or directory)."""
    if is_directory:
        # If it's a directory, ensure it exists or create it
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
        return path
    else:
        # If it's a file, ensure its parent directory exists or create it
        parent_dir = os.path.dirname(path)
        if parent_dir and not os.path.exists(parent_dir):
            os.makedirs(parent_dir, exist_ok=True)
        return path

def load_positive_list(file_path: Optional[str]) -> Set[str]:
    """Load a list of positive HMMs from a CSV file."""
    if file_path:
        return set(pd.read_csv(file_path, header=None)[0].tolist())
    return set()

def _translate_to_stop(seq: str, genetic_code: Dict[str, str], threshold: int, info: str) -> List[Tuple[str, str]]:
    """Translate a sequence until a stop codon is encountered and attach metadata."""
    proteins = []
    skews = {0: len(seq), 1: len(seq), 2: len(seq)}
    for match in re.finditer(r'TGA|TAA|TAG', seq):
        skew = match.start() % 3
        skews[skew] = min(skews[skew], match.start())

    for skew in [key for key, value in skews.items() if value >= threshold]:
        protein = []
        seq_trunc = seq[skew:]
        end = min(skews[skew], len(seq_trunc))

        for i in range(0, (end // 3) * 3, 3):
            aa = genetic_code.get(seq_trunc[i:i + 3], 'X')
            protein.append(aa)
        proteins.append(("".join(protein), f'{info};skew_{skew};'))
    return proteins

def _translate_non_stop(seq: str, genetic_code: Dict[str, str], threshold: int, info: str) -> List[Tuple[str, str]]:
    """Translate the entire sequence, including stop codons, and attach metadata."""
    proteins = []
    for skew in [0, 1, 2]:
        protein = []
        seq_trunc = seq[skew:]
        for i in range(0, (len(seq_trunc) // 3) * 3, 3):
            aa = genetic_code.get(seq_trunc[i:i + 3], 'X')
            protein.append(aa)
        proteins.append(("".join(protein), f'{info};skew_{skew};'))
    return proteins

def translate_sequence(seq: str, genetic_code: Dict[str, str], threshold: int, to_stop: bool, info: str) -> List[Tuple[str, str]]:
    """Translate a DNA sequence into a protein sequence and return it with metadata."""
    if to_stop:
        return _translate_to_stop(seq, genetic_code, threshold, info)
    else:
        return _translate_non_stop(seq, genetic_code, threshold, info)

def _split_sequence(seq: str, info: str, max_length: int = 300000) -> List[Tuple[str, str]]:
    """Split a sequence into smaller parts and attach metadata."""
    return [(seq[i:i + max_length], f'{info};part_{i};len_{len(seq[i:i + max_length])};len_contig_{len(seq)}') for i in range(0, len(seq), max_length)]

def analyze_sequences(hmms: List[pyhmmer.plan7.HMM], sequences: List[pyhmmer.easel.DigitalSequence], scores_global: np.ndarray, threads: int) -> np.ndarray:
    """Analyze sequences using HMM profiles."""
    hits = pyhmmer.hmmer.hmmsearch(hmms, sequences, cpus=threads)
    for hit in hits:
        for domain in hit:
            try:
                scores = np.array([(hit.query_name, domain.name, i.score, i.alignment.target_from, i.alignment.target_to) for i in domain.domains])
                scores_global = np.concatenate([scores_global, scores])
            except ValueError:
                pass
    return scores_global

def analyze_file(seq_file: str, hmms: List[pyhmmer.plan7.HMM], threshold: int, batch_size: int, file_type: str, to_stop: bool, threads: int) -> pd.DataFrame:
    """Analyze sequences in a file using HMM profiles."""
    sequences = []
    scores_global = np.empty([0, 5])

    for i, dna_record in enumerate(SeqIO.parse(seq_file, file_type)):
        # Split the sequence into smaller parts
        dna_seqs = _split_sequence(str(dna_record.seq), f'{dna_record.name};chain_+1') + _split_sequence(str(dna_record.seq.reverse_complement()), f'{dna_record.name};chain_-1')

        # Translate each DNA sequence into protein sequences
        aa_seqs = []
        for seq, info in dna_seqs:
            translated = translate_sequence(seq, GENETIC_CODE, threshold, to_stop, info)
            aa_seqs.extend(translated)

        # Debug: Print the structure of aa_seqs
        print(f"aa_seqs structure: {aa_seqs}")

        # Convert protein sequences to pyhmmer-compatible format
        digitized_seqs = [pyhmmer.easel.TextSequence(sequence=seq, name=bytes(info, 'utf-8')).digitize(hmms[0].alphabet) for seq, info in aa_seqs]
        sequences.extend(digitized_seqs)

        # Analyze sequences in batches
        if (i + 1) % batch_size == 0:
            scores_global = analyze_sequences(hmms, sequences, scores_global, threads)
            sequences.clear()

    # Analyze any remaining sequences
    if sequences:
        scores_global = analyze_sequences(hmms, sequences, scores_global, threads)
    # Create the DataFrame without specifying dtype
    result_frame = pd.DataFrame(scores_global, columns=['Query', 'Name', 'Score', 'From', 'To'])
    # Set the data types for each column using astype
    result_frame = result_frame.astype({
		'Query': 'str',        # Ensure Query is string
		'Name': 'str',         # Ensure Name is string
		'Score': 'float32',    # Ensure Score is 32-bit float
		'From': 'int32',       # Ensure From is 32-bit integer
		'To': 'int32'          })
    return result_frame


def load_hmms(hmm_path: str, positive_list: Optional[Set[str]] = None) -> List[pyhmmer.plan7.HMM]:
    """Load HMM profiles from a file or directory."""
    hmms = []
    if os.path.isfile(hmm_path):
        # Load a single HMM file
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            hmms.append(hmm_file.read())
    elif os.path.isdir(hmm_path):
        # Load all HMM files in the directory
        for hmm_file in glob.glob(os.path.join(hmm_path, '*')):
            try:
                with pyhmmer.plan7.HMMFile(hmm_file) as f:
                    hmms.append(f.read())
            except Exception as e:
                print(f"Error loading HMM file {hmm_file}: {e}")
    else:
        raise ValueError(f"Invalid HMM path: {hmm_path}")

    # Filter HMMs based on the positive list
    if positive_list:
        hmms = [hmm for hmm in hmms if str(hmm.name.decode("utf-8")) in positive_list]

    return hmms

def main():
    args = parse_arguments()

    # Validate input paths
    input_path = validate_input_path(args.file)
    hmm_path = validate_input_path(args.hmm)

    # Handle output path
    if os.path.isdir(args.output):
        # If output is a directory, ensure it exists
        output_path = validate_input_path(args.output, is_directory=True)
    else:
        # If output is a file, ensure its parent directory exists
        output_path = validate_input_path(args.output)

    # Load positive list if provided
    positive_list = load_positive_list(args.positive_list)

    # Load HMM profiles
    hmms = load_hmms(hmm_path, positive_list)

    # Analyze sequences
    results = {}
    for file in glob.glob(input_path):
        file_type = _get_file_type(file)
        file_name = os.path.basename(file).rsplit('.', 2)[0] if file.endswith('.gz') else os.path.basename(file).rsplit('.', 1)[0]
        print(f'Analyzing file: {file_name}')
        results[file] = analyze_file(file, hmms, args.threshold, args.batch, file_type, args.to_stop, args.threads)

        # Save results
        if os.path.isfile(args.file):
            results[file].to_csv(output_path)
        else:
            results[file].to_csv(os.path.join(output_path, f'{file_name}.csv'))

def _get_file_type(file_path: str) -> str:
    """Determine the file type based on the extension."""
    if file_path.endswith('.gz'):
        ext = file_path.split('.')[-2]
    else:
        ext = file_path.split('.')[-1]
    if ext in ['fasta', 'fna', 'ffn', 'faa', 'frn', 'fa']:
        return 'fasta'
    elif ext in ['fastq', 'fq']:
        return 'fastq'
    else:
        raise ValueError(f'Unsupported file type: {ext}')

if __name__ == '__main__':
    main()
