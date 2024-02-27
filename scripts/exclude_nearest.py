import pandas as pd
import argparse

description = """ Thi script analyse phylip matrix of distances and exclude seqeunces that are the nearest to the sequence of interest)"""
        
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--input_seq' help='list of files')
parser.add_argument('-o', '--output', type=str, help='An output fasta file')
parser.add_argument('-m', '--distance_matrix', type=str, help='An phylip distance matrix')
parser.add_argument('-n', '--n_seq', type=int, help='Maximum amout of sequnces that should be included to MSA')
args = parser.parse_args()

input_seq = args.input_seq
out = args.output
distance_matrix = args.distance_matrix
n_seq = args.n_seq




df = pd.read_csv(distance_matrix, skiprows=5, sep='\t', index_col=None).reset_index(level=1)
df = df[df.columns[::-1]]

df = df.set_index(df.iloc[:, 0]).iloc[:, 2:]

df.columns = list(df.index[::-1])


interest = [i.split(' ')[0] for i in df.iloc[0, :].sort_values().index][:n_seq]

seqs= []
with open(input_seq, "r") as file:
    sequences = SeqIO.parse(file, "fasta")

    # Iterate over each sequence in the FASTA file
    for sequence in sequences:
        if sequence.id in interest:
            sequence.seq = sequence.seq.replace('-', '')
            seqs.append(sequence)
        

SeqIO.write(seqs, output, 'fasta')

