import subprocess
import glob
from Bio import SeqIO
import pandas as pd
import argparse
import os

description = """The script gets fasta file with several sequences and select n sequences that are the most similar to the first with blast. The selected sequnces are written to output fasta. The script creates several tmp files in a folder of the output file """

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-i', '--input_file', type=str, help='An input file')
parser.add_argument('-o', '--output_file', type=str, help='An output file')
parser.add_argument('-n', '--n_seqs', type=int, help="Number of sequences that should be selected")

args = parser.parse_args()



if args.n_seqs:
    n_seqs = int(args.n_seqs)
else:
    n_seqs = 50



# Define the BLAST command

file = args.input_file
output_file = args.output_file

name_output = ".".join(output_file.split('.')[:-1])
temp_query_file = f"{name_output}.tmp.query.fasta"
temp_database_file = f"{name_output}.tmp.db"
blast_results = f"{name_output}.blast.tsv"

# Load your sequence of interest
seq_io = SeqIO.parse(file, "fasta")
sequence_of_interest = next(seq_io)
SeqIO.write(sequence_of_interest, temp_query_file, 'fasta')
make_db = f"makeblastdb -in {file} -dbtype prot -out {temp_database_file}"
subprocess.run(make_db, shell=True)

blast_cmd = f"blastp -query {temp_query_file} -db {temp_database_file} -out {blast_results} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -gapopen 11 -gapextend 1 -evalue 100"

# Run the BLAST command using subprocess
subprocess.run(blast_cmd, shell=True)

for tm_file in glob.glob(f"{name_output}.tmp*"):
    os.remove(tm_file)

# Read the tabular results file into a pandas DataFrame

column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

df = pd.read_csv(blast_results, sep='\t', header=None, names=column_names).sort_values(by='bitscore', ascending=False).drop_duplicates(subset='sseqid')

# Display the DataFrame
selected_ids = df.head(n_seqs+1).sseqid.to_list()
its_too_short = True if len(selected_ids) < 10 else False

if its_too_short:
    i=0
    seqs = []
    with open(file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            i+=1
            seqs.append(record)
            if i > 10:
                break
    SeqIO.write(seqs, output_file, 'fasta')
else:
    i=0
    seqs = []
    with open(file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in selected_ids or its_too_short:
                seqs.append(record)
                
    SeqIO.write(seqs, output_file, 'fasta')
