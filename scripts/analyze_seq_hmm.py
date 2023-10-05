from Bio import SeqIO
from Bio.Seq import Seq
import os
import pyhmmer
import glob
import sys
import numpy as np
import pandas as pd
import re
import gzip
import argparse

description = """The script gets a fastq or a fasta file, translate it and analyzes sequences with HMM profiles"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-f', '--file', type=str, help='A file or directory with fastq/fasta extenstion')
parser.add_argument('-o', '--output', type=str, help='An output folder')
parser.add_argument('-m', '--hmm', type=str, help='A file or directory with HMM profile')
parser.add_argument('-t', '--threads', type=str, help='Amount of threads to use')
parser.add_argument('-b', '--batch', type=str, help='Batch size in which sequences will be analyzed')
parser.add_argument('-s', '--to_stop', type=str, help='Translate sequences only to stop codon')
parser.add_argument('-r', '--threshold', type=str, help='Minimum amino acid sequence length for analysis')

args = parser.parse_args()

def if_condition(x, message):
    if not x:
        print(message)
        exit()

print(args)
    

if args.file:
    if os.path.isfile(args.file):
        file = args.file
    elif os.path.isdir(args.file):
        file = f'{args.file}/*'
    else:
        if_condition(False, "An input file does not exists")
else:
    if_condition(False, 'Missed -f argument')


if args.hmm:
    if os.path.isfile(args.hmm):
        hmm_path = args.hmm
    elif os.path.isdir(args.hmm):
        hmm_path = f'{args.hmm}/*'
    else:
        if_condition(False, "An HMM file does not exists")
else:
    if_condition(False, 'Missed -m argument')
    

if args.output:
    if os.path.isfile(args.file):
        output = args.output
    elif os.path.isdir(args.file):
        output_dir = args.output
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)    
else:
    if_condition(False, 'Missed -o argument')

    
if args.threads:
    threads = int(args.threads)
else:
    threads = 1
    
if args.to_stop:
    to_stop = args.to_stop
else:
    to_stop = False
    
if args.batch:
    batch = int(args.batch)
else:
    batch = 1000
   

if args.threshold:
    threshold = int(args.threshold)
else:
    threshold = 50

nnn_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA': '*', 'TAG': '*',
        'TGC':'C', 'TGT':'C', 'TGA': '*', 'TGG':'W',
    }

def translate(seqs, nnn_table, threshold=90, to_stop=True):
    """
    This function tranlate a list of sequence and return only translations that are longer than threshold
    seqs
    
    """
    proteins = []
    if to_stop:
        translate_func = _translate_to_stop
    else:
        translate_func = _translate_non_stop
    
    for seq in seqs:
        proteins.append(translate_func(seq, nnn_table, threshold=threshold))
    
    return np.concatenate(proteins)

def _translate_to_stop(seq_array, nnn_table, threshold):
    
    
    proteins= []
    seq, seq_range= seq_array
    count = 0
    seq_len = len(seq)
    skews = {0:seq_len, 1:seq_len, 2:seq_len}
    for match in re.finditer(r'TGA|TAA|TAG', seq):
        count += 1
        skew = match.start() % 3
        skews[skew] = min(skews[skew], match.start())

    skews_a = [key for key, value in skews.items() if value >= threshold]

    for skew in skews_a:

        protein = []
        count = 0
        seq_trunc = seq[skew:]

        end = min(skews[skew], len(seq_trunc))
        


        for i in range(0, (end // 3) * 3, 3):
            aa = nnn_table.get(seq_trunc[i:i + 3], 'X')
            protein.append(aa)
        else:
            proteins.append(["".join(protein)])
    return proteins

def _translate_non_stop(seq_array, nnn_table, threshold):
    proteins= []
    count = 0
    
    seq, seq_info = seq_array
    
    skews_a = [0, 1, 2]

    for skew in skews_a:

        protein = []
        count = 0
        seq_trunc = seq[skew:]
        end = len(seq_trunc)


        for i in range(0, (end // 3) * 3, 3):
            aa = nnn_table.get(seq_trunc[i:i + 3], 'X')
            protein.append(aa)
        else:
            proteins.append(["".join(protein), seq_info+f'{skew}'])
    return np.array(proteins)    

def separate_string(string, max_length, info):
    parts = []
    for i in range(0, len(string), max_length):
        parts.append([string[i:i+max_length], f'{info}.{i}_'])
    return parts

def analyse_seqs(hmms, seqs, scores_global):
    
    global threads
    
    hits = pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=threads)
    for hit in hits:
        for hit2 in hit:
            scores = np.array([(hit.query_name, hit2.name, i.score, i.alignment.target_from, i.alignment.target_to) for i in hit2.domains])

    try:
        scores_global = np.concatenate([scores_global, scores])
    except ValueError:
        pass
    return scores_global

def analyse_file(seq_file, hmms, threshold=90, batch=100000, type_of_file='fastq', to_stop=False):
    """
    Get a fasta/fastq file
    Iterate through seqs
    Split reads to parts less then 100000
    Translate reads. By default, add to futher analysis only translations with distance to stop-codon bigger then threshold
    Runs reads through the HMMer3 with pyhmmerv in batchs size of "batch" parameter
    Return a dataframe with columns: 'Query'(name of HMM), 'Name'(name of read), 'Score'(best domain score)
    
    seq_file - seq_file
    hmms - list of hmms of pyhmmer format
    threshold - minimum len of translated sequence
    batch - size of batch that go to the pyhmmer
    type_of_file - fastq/fasta
    to_stop - translate sequence only to stop-codon(True) or translate entire seq with stop-codon(False)
    
    """
    
    
    fastq_seq = SeqIO.parse(seq_file, type_of_file)
    i=0
    c=0
    scores_global = np.empty([0, 5])
    seqs = []
    for dna_record in fastq_seq:
        i+=1
        dna_seqs = np.concatenate([separate_string(str(dna_record.seq), 290000, f'{dna_record.name}_dir'), separate_string(str(dna_record.seq.reverse_complement()), 290000, f'{dna_record.name}_rev')])
        

        # generate all translation frames
        aa_seqs = translate(dna_seqs, nnn_table, threshold=threshold, to_stop=to_stop)
        aa_seqs = [pyhmmer.easel.TextSequence(sequence=str(i[0]), name=bytes(i[1], 'utf-8')).digitize(hmm.alphabet) for i in aa_seqs]

        seqs += aa_seqs
        c+=len(aa_seqs)

        if i % batch == 0:

            scores_global = analyse_seqs(hmms, seqs, scores_global)
            del seqs
            seqs = []
    scores_global = analyse_seqs(hmms, seqs, scores_global)
    del seqs
    
    result_frame = pd.DataFrame(np.array(scores_global), columns=['Query', 'Name', 'Score', 'From', 'To'])
    result_frame.Score = result_frame.Score.astype('float')
    result_frame.Name = result_frame.Name.astype('str')
    result_frame[['From', 'To']] = result_frame[['From', 'To']].astype('int') * 3
    result_frame.Query = result_frame.Query.astype('str')
    return result_frame

def check_type_of_file(extension):
    if extension in ['fasta', 'fna', 'ffn', 'faa', 'frn', 'fa']:
        return 'fasta'
    elif extension in ['fastq', 'fq']:
        return 'fastq'
    else:
        raise  OSError('Input file must be a fasta or fastq')

hmms = []
list_of_files = glob.glob(hmm_path)
for i in list_of_files:
    with pyhmmer.plan7.HMMFile(i) as hmm_file:
        hmms.append(hmm_file.read())
hmm = hmms[0]


results_dict = {}
for i in glob.glob(file):
    if i.split('.')[-1] == 'gz':
        i_open = gzip.open(i, "rt")
        type_of_file = check_type_of_file(i.split('.')[-2])
        file_name =  '.'.join(i.split("/")[-1].split(".")[:-2])
    else:
        i_open = i
        file_name =  '.'.join(i.split("/")[-1].split(".")[:-1])
        type_of_file = check_type_of_file(i.split('.')[-1])
    print(f'Iterate HMMs through file "{i.split("/")[-1].split(".")[0]}"')
    results_dict[i] = analyse_file(i_open, hmms,threshold=threshold, batch=batch, type_of_file=type_of_file, to_stop=to_stop)
     
    if os.path.isfile(args.file):
        results_dict[i].to_csv(f'{output}')
    elif os.path.isdir(args.file):
        results_dict[i].to_csv(f'{output_dir}/{file_name}.csv')
    