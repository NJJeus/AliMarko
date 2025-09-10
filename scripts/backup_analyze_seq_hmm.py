import os
import sys
import glob
import gzip
import re
import argparse
from typing import List, Tuple, Set, Optional, Dict, Any
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import pyhmmer

# Constants
BATCH_SIZE = 1000
THRESHOLD = 50
MAX_SEQ_LENGTH = 300000

# Genetic Code Table
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
    parser.add_argument('-t', '--threads', type=int, default=1, help='Amount of threads to use')
    parser.add_argument('-b', '--batch', type=int, default=BATCH_SIZE, help='Batch size for sequence analysis')
    parser.add_argument('-s', '--to_stop', action='store_true', help='Translate sequences only to stop codon')
    parser.add_argument('-r', '--threshold', type=int, default=THRESHOLD, help='Minimum amino acid sequence length for analysis')
    parser.add_argument('-l', '--positive_list', type=str, help='A CSV file with HMMs to use. It should have only one column')
    return parser.parse_args()

def validate_path(path: str) -> str:
    """Validate if the given path exists (file or directory)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Path does not exist: {path}")
    return path

def ensure_output_dir(output_dir: str) -> None:
    """Ensure the output directory exists. If not, create it."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def load_positive_list(file_path: Optional[str]) -> Set[str]:
    """Load a set of HMM names from a CSV file."""
    if not file_path:
        return set()
    return set(pd.read_csv(file_path, header=None)[0].tolist())

def load_hmms(hmm_path: str) -> List[Any]:
    """Load HMM profiles from a file or directory."""
    if os.path.isfile(hmm_path):
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            return [hmm_file.read()]
    elif os.path.isdir(hmm_path):
        hmms = []
        for hmm_file in glob.glob(os.path.join(hmm_path, '*')):
            with pyhmmer.plan7.HMMFile(hmm_file) as hmm_file_obj:
                hmms.append(hmm_file_obj.read())
        return hmms
    else:
        raise ValueError(f"Invalid HMM path: {hmm_path}")

def translate_sequence(seq: str, genetic_code: Dict[str, str], to_stop: bool = True, threshold: int = THRESHOLD) -> List[str]:
    """Translate a DNA sequence into a protein sequence."""
    proteins = []
    seq_len = len(seq)
    skews = {0: seq_len, 1: seq_len, 2: seq_len}

    if to_stop:
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
        proteins.append("".join(protein))
    return proteins

def analyze_sequences(hmms: List[Any], sequences: List[Any], threads: int) -> np.ndarray:
    """Analyze sequences using HMM profiles."""
    scores_global = np.empty([0, 5])
    hits = pyhmmer.hmmer.hmmsearch(hmms, sequences, cpus=threads)
    for hit in hits:
        for hit2 in hit:
            try:
                scores = np.array([(hit.query_name, hit2.name, i.score, i.alignment.target_from, i.alignment.target_to) for i in hit2.domains])
                scores_global = np.concatenate([scores_global, scores])
            except ValueError:
                pass
    return scores_global

def process_file(file_path: str, hmms: List[Any], threshold: int, batch: int, file_type: str, to_stop: bool, threads: int) -> pd.DataFrame:
    """Process a single file and analyze its sequences."""
    sequences = []
    scores_global = np.empty([0, 5])
    with gzip.open(file_path, "rt") if file_path.endswith('.gz') else open(file_path, "r") as file:
        for i, record in enumerate(SeqIO.parse(file, file_type)):
            dna_seqs = [str(record.seq), str(record.seq.reverse_complement())]
            for dna_seq in dna_seqs:
                aa_seqs = translate_sequence(dna_seq, GENETIC_CODE, to_stop=to_stop, threshold=threshold)
                # Use the alphabet from the first HMM in the list
                aa_seqs = [pyhmmer.easel.TextSequence(sequence=seq, name=bytes(f"{record.name};skew_{skew}", 'utf-8')).digitize(hmms[0].alphabet) for skew, seq in enumerate(aa_seqs)]
                sequences.extend(aa_seqs)
            if (i + 1) % batch == 0:
                scores_global = analyze_sequences(hmms, sequences, threads)
                sequences = []
    if sequences:
        scores_global = analyze_sequences(hmms, sequences, threads)
    return pd.DataFrame(scores_global, columns=['Query', 'Name', 'Score', 'From', 'To'])

def main():
    args = parse_arguments()
    validate_path(args.file)
    validate_path(args.hmm)
    ensure_output_dir(os.path.basename(os.path.dirname(args.output)))  # Ensure the output directory exists

    positive_list = load_positive_list(args.positive_list)
    hmms = load_hmms(args.hmm)
    if positive_list:
        hmms = [hmm for hmm in hmms if str(hmm.name.decode("utf-8")) in positive_list]

    results = {}
    for file_path in glob.glob(args.file):
        file_type = 'fasta' if any(file_path.endswith(ext) for ext in ['fasta', 'fna', 'ffn', 'faa', 'frn', 'fa']) else 'fastq'
        results[file_path] = process_file(file_path, hmms, args.threshold, args.batch, file_type, args.to_stop, args.threads)
        # Construct the output file path
        file_name = os.path.basename(os.path.dirname(file_path)) + '.csv'
        output_path = str(args.output)
        print(output_path)
        results[file_path].to_csv(output_path, index=False)

if __name__ == "__main__":
    main()
