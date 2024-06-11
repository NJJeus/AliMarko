
    
import pandas as pd
from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqIO import SeqRecord
import subprocess
import argparse
import os


"""
# Дизайн кода
1. Вытаскиваем множество пары HMM-контиги, и n сработавших HMM выполняем трансляцию
2. Вытаскииваем nHMM-геномы
> При этом по GENBANK ID мы вытаскиваем название вируса и род. Эт делается с помощью словаря {GENBAMK:VIRUS_GEN}. При создании словаря неуникальные ключи объединяются
3. Складываем геномы и контиги в MSA, контиги по одному. Файл "контиг-HMM". Если она HMM сработала несколько раз, сшиваем последовательности
> Далее надо будет что-то придумать с перекрывающимися контигами от разных семейств  
> Далее надо будет придумать, что делать, если одна HMM сработала на разные контиги
4. Запускаем множественные выравнивания. 
5. Обрезаем множественные выравнивания. Пока что под чистую удаляем гепы
> Потом надо будет придумать нормальную систему тримминга
5. Рисуем деревья "контиг-HMM"
6. Вставляем рисунки в отчёт
"""


def if_condition(x, message):
    if not x:
        print(message)
        exit()

description = """A script"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-g', '--genome_reference', type=str, help='')
parser.add_argument('-c', '--contig_fasta', type=str, help='')
parser.add_argument('-i', '--ictv_report', type=str, help='')
parser.add_argument('-r', '--contig_report', type=str, help='')
parser.add_argument('-o', '--output', type=str, help='Output folder')
parser.add_argument('-t', '--ictv_taxo', type=str, help='')

args = parser.parse_args()




if args.genome_reference:
    if_condition(os.path.isfile(args.genome_reference), "coverage file doesn't exist")
    reference_fasta  = args.genome_reference
else:
    if_condition(False, 'Missed -g argument')
    
if args.contig_fasta:
    if_condition(os.path.isfile(args.contig_fasta), "coverage file doesn't exist")
    contigs_fasta  = args.contig_fasta
else:
    if_condition(False, 'Missed -c argument')

if args.ictv_report:
    if_condition(os.path.isfile(args.ictv_report), "coverage file doesn't exist")
    ictv_report  = args.ictv_report
else:
    if_condition(False, 'Missed -c argument')  
                 
if args.contig_report:
    if_condition(os.path.isfile(args.contig_report), "coverage file doesn't exist")
    contigs_report = args.contig_report
else:
    if_condition(False, 'Missed -c argument')  

if args.ictv_taxo:
    if_condition(os.path.isfile(args.ictv_taxo), "coverage file doesn't exist")
    ictv_taxo  = args.ictv_taxo
else:
    
    if_condition(os.path.isfile('ictv_tables/ictv_taxo.xlsx'), "coverage file doesn't exist")    
    ictv_taxo = 'ictv_tables/ictv_taxo.xlsx'
   
if args.output:
    out_folder = args.output
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)    
else:
    if_condition(False, 'Missed -o argument')                 





def is_intersecting(segment1, segment2):
    segment1[0], segment1[1] = sorted(segment1)
    segment2[0], segment2[1] = sorted(segment2)
    return not (segment1[1]*1.1 <= segment2[0] or segment1[0]*1.1 >= segment2[1])

def remove_intersecting_segments(segments):
    """
    Gets list in format [[index, [1, 2]], [index, [3, 4]]]
    """
    result = {segments[0][0]: segments[0][1]}  # Keep the first segment with index
    saved_indices = [segments[0][0]]
    
    for i in range(1, len(segments)):
        current_index, current_segment = segments[i]
        is_intersects = False
        
        for saved_index in saved_indices:
            if is_intersecting(current_segment, result[saved_index]):
                is_intersects = True
                break
        
        if not is_intersects:
            result[current_index] = current_segment
            saved_indices.append(current_index)
    
    return result, saved_indices


ictv = pd.read_csv(ictv_report, index_col=0)
ictv['Name'] = ictv['Name'].apply(lambda x: x.split('|')[1])

ictv_data = pd.read_excel(ictv_taxo)
ictv_data = ictv_data.dropna(subset='Virus GENBANK accession')
ictv_data['genbank_list'] = ictv_data['Virus GENBANK accession'].apply(lambda i :
                                                                       [el.split(":")[-1] for el in i.replace(' ', '').split(';')])
ictv_data['Isolate_id'] = ictv_data['Species'] +  '{' + ictv_data['Sort'].astype('str')  + '}' + '|' + ictv_data['Genus']
indeces = []
accessions = []
for index, item in ictv_data.iterrows():
    for accession in item['genbank_list']:
        indeces.append(item.Isolate_id)
        accessions.append(accession)

        
taxo_index = pd.DataFrame({'ictv_taxo_index':indeces, 'accessions':accessions}).drop_duplicates(subset='accessions').set_index('accessions')
taxo_index = taxo_index.to_dict(orient='index')

contigs = pd.read_csv(contigs_report, index_col=0)
contigs = contigs.query('Lengh > 300')

contigs_names = contigs['Name'].unique()

name_c_contig = contigs_names[0]
for name_c_contig in contigs_names:

    original_array = contigs.query(f"Name == '{name_c_contig}'").reset_index()[['index', 'From', 'To']].to_numpy()
    reshaped_array = [[i[0], [i[1], i[2]]] for i in original_array]
    filtered_segments, saved_indices = remove_intersecting_segments(reshaped_array)
    try:
        contigs_sample = contigs.iloc[saved_indices]
    except Exception:
        contigs_sample = contigs.head(3)
        

    contig_rows = contigs_sample[['From', 'To', 'Frame', 'Query']].to_dict(orient='index')
    

    hmms = list(contigs_sample.Query.unique())


    for cur_hmm in hmms[:3]:

        ictv_names = ictv.query(f'Query == "{cur_hmm}"')[['Name', 'From', 'To', 'Frame']].reset_index()


        genome_seqs = []

        # Open the FASTA file and parse its content
        with open(reference_fasta, "r") as file:
            sequences = SeqIO.parse(file, "fasta")
            i=0

            # Iterate over each sequence in the FASTA file
            for sequence in sequences:
                seq_id = sequence.id.split('|')[1]
                
                ictv_matches = [ictv_match for i, ictv_match in ictv_names.iterrows() if str(ictv_match['Name']) == seq_id]
                match_count=0
                for ictv_match in ictv_matches:
                    
                    From, To, Frame = ictv_match['From'], ictv_match['To'], ictv_match['Frame']
                    if From > To:
                        To, From = From, To
                    shift_dict = {1:0, 2:1, 3:2, 4:0, 5:1, 6:2}
                    shift = shift_dict[Frame]
                    try:
                        if Frame in [1, 2, 3]:
                             translated_seq = sequence.seq[From+shift:To].translate()
                        else:
                            translated_seq = sequence.seq[From:To+1].reverse_complement()[shift:].translate()
                    except Exception:
                        print(translated_seq[From+shift:To])
                        os._exit(1)
                        
                    output_seq = SeqRecord(translated_seq, id=str(taxo_index[seq_id]['ictv_taxo_index']).replace(' ', '_')+f"#{match_count}")

                    genome_seqs.append(output_seq)
                    match_count+=1

        
        contig_seqs = []

        # Open the FASTA file and parse its contents
        with open(contigs_fasta, "r") as file:
            sequences = SeqIO.parse(file, "fasta")
            i=0

            # Iterate over each sequence in the FASTA file
            for sequence in sequences:
                seq_id = sequence.id
                if seq_id != name_c_contig:
                    continue

                for k, values in contig_rows.items(): 
                    if contig_rows[k]['Query'] != cur_hmm:
                        continue
                    From, To, Frame = contig_rows[k]['From'], contig_rows[k]['To'], contig_rows[k]['Frame']
                    if From > To:
                        To, From = From, To
                    shift_dict = {1:0, 2:1, 3:2, 4:0, 5:1, 6:2}
                    shift = shift_dict[Frame]
                    if Frame in [1, 2, 3]:
                        aa_seq = SeqRecord(sequence.seq[From+shift:To].translate())
                    else:
                        aa_seq = SeqRecord(sequence.seq[From:To].reverse_complement()[shift:].translate())
                    aa_seq.id = sequence.name + "__"  + contig_rows[k]['Query'] + '__' + str(k)
                    contig_seqs.append(aa_seq)
                    print(f"ID: {seq_id}")

        for c_seq in contig_seqs:
            SeqIO.write([c_seq]+genome_seqs, f'{out_folder}/test_{c_seq.id}.fasta', 'fasta')
