import subprocess
import glob
from Bio import SeqIO
import pandas as pd
import argparse
import os

def parse_arguments():
    """Parse command-line arguments."""
    description = """The script gets a FASTA file with several sequences and selects n sequences that are the most similar to the first sequence using BLAST. The selected sequences are written to an output FASTA file. The script creates several temporary files in the folder of the output file."""
    
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output FASTA file')
    parser.add_argument('-n', '--n_seqs', type=int, default=30, help='Number of sequences to select (default: 30)')
    
    return parser.parse_args()

def create_blast_database(input_file, database_file):
    """Create a BLAST database from the input FASTA file."""
    make_db_cmd = f"makeblastdb -in {input_file} -dbtype prot -out {database_file}"
    subprocess.run(make_db_cmd, shell=True, check=True)

def run_blast(query_file, database_file, blast_results):
    """Run BLASTp and save the results."""
    blast_cmd = (
        f"blastp -query {query_file} -db {database_file} -out {blast_results} "
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' "
        "-gapopen 11 -gapextend 1 -evalue 100"
    )
    subprocess.run(blast_cmd, shell=True, check=True)

def read_blast_results(blast_results, n_seqs):
    """Read BLAST results and select top n sequences."""
    column_names = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    
    df = pd.read_csv(blast_results, sep='\t', header=None, names=column_names)
    df = df.sort_values(by='bitscore', ascending=False).drop_duplicates(subset='sseqid')
    
    return df.head(n_seqs + 1)['sseqid'].tolist()

def write_selected_sequences(input_file, output_file, selected_ids):
    """Write selected sequences to the output FASTA file."""
    sequences = []
    with open(input_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in selected_ids:
                sequences.append(record)
    
    SeqIO.write(sequences, output_file, 'fasta')

def cleanup_temp_files(temp_files):
    """Remove temporary files."""
    for temp_file in temp_files:
        os.remove(temp_file)

def main():
    args = parse_arguments()
    
    # Define temporary file names
    name_output = os.path.splitext(args.output_file)[0]
    temp_query_file = f"{name_output}.tmp.query.fasta"
    temp_database_file = f"{name_output}.tmp.db"
    blast_results = f"{name_output}.blast.tsv"
    
    # Extract the first sequence as the query
    with open(args.input_file, 'r') as handle:
        sequence_of_interest = next(SeqIO.parse(handle, 'fasta'))
    SeqIO.write(sequence_of_interest, temp_query_file, 'fasta')
    
    # Create BLAST database and run BLAST
    create_blast_database(args.input_file, temp_database_file)
    run_blast(temp_query_file, temp_database_file, blast_results)
    
    # Read BLAST results and select sequences
    selected_ids = read_blast_results(blast_results, args.n_seqs)
    
    # Write selected sequences to the output file
    if len(selected_ids) < 10:
        print("Warning: Fewer than 10 sequences were selected. Writing the first 10 sequences instead.")
        with open(args.input_file, 'r') as handle:
            sequences = [record for _, record in zip(range(10), SeqIO.parse(handle, 'fasta'))]
        SeqIO.write(sequences, args.output_file, 'fasta')
    else:
        write_selected_sequences(args.input_file, args.output_file, selected_ids)
    
    # Clean up temporary files
    cleanup_temp_files(glob.glob(f"{name_output}.tmp*"))

if __name__ == "__main__":
    main()
