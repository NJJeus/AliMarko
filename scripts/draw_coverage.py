import os
import pandas as pd
import argparse

def if_condition(x, message):
    if not x:
        print(message)
        exit()

description = """ Script that concatenate a samtools coverage output and a samtools view quality with an ictv metadata"""
        
parser = argparse.ArgumentParser(description=description)
parser.add_argument('-c', '--coverage', type=str, help='A file with samtools coverage output')
parser.add_argument('-o', '--output', type=str, help='An output file')
parser.add_argument('-t', '--tables_folder', type=str, help='A folder with tables ictv_taxo and genbank accessions')
parser.add_argument('-r', '--threshold', type=str, help='A threshold for drawing virus coverage. 0.2 by default')

args = parser.parse_args()

if args.coverage:
    if_condition(os.path.isfile(args.coverage), "coverage file doesn't exist")
    coverage_file = args.coverage
else:
    if_condition(False, 'Missed -c argument')
    
if args.output:
    output_file = args.output
else:
    if_condition(False, 'Missed -o argument')
if args.tables_folder:
    if_condition(os.path.isdir(args.tables_folder), "tables folder doesn't exist")
    tables_folder = args.tables_folder
else:
    if_condition(False, 'Missed -t argument') 

if args.threshold:
    threshold = args.threshold
else:
    threshold = 0.05

dirname = os.path.dirname(os.path.normpath(output_file))
print(dirname)
if not os.path.exists(dirname):
    os.makedirs(dirname)    
    
ictv_data = pd.read_excel(f'{tables_folder}/ictv_taxo.xlsx')
taxo_index = pd.read_csv(f'{tables_folder}/genbank_accessions.csv', index_col=0)    
    
def get_data_table(coverage_file):
    data = pd.read_table(coverage_file).rename(columns={'#rname': 'rname'}).drop(['numreads', 'covbases',
           'meandepth', 'meanbaseq', 'meanmapq'], axis=1)

    data['rname_short'] = data['rname'].apply(lambda x: x.split('|')[1])
    data = data.set_index('rname_short')

    data = data.copy()
    data = pd.merge(left=data.reset_index(), 
                    right=taxo_index.reset_index(), left_on='rname_short', right_on='index') # Get index of all accession number
    data = pd.merge(left=data, right=ictv_data.reset_index(), 
                    left_on='ictv_taxo_index', right_on='index') # Get ictv metadata by accession number


    # Generalize fragments of one virus
    data[['coverage']] = data[['coverage']]  * 0.01
    data['len'] = data.endpos - data.startpos
    data[['weighted_coverage']] = (data[['coverage']].T * data['len']).T


    agg_all_columns_dict = dict()
    agg_all_columns_dict.update({'len':'sum', 'weighted_coverage':'sum'})
    agg_all_columns_dict.update({'coverage': lambda x: [i for i in list(x)]})
    agg_all_columns_dict.update({'endpos': lambda x: [i for i in list(x)]})
    agg_all_columns_dict.update({'rname': lambda x: [i for i in list(x)]})
    agg_all_columns_dict.update({i:'first' for i in ['Realm', 'Kingdom',
           'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Order', 'Family', 'Genus', 'Species',
                  'Virus name(s)', 'Virus GENBANK accession', 'Host source']})
    data = data.groupby('Species').agg(agg_all_columns_dict)
    data['fragments_coverage'] = data['coverage'].copy()


    data['coverage'] = data.weighted_coverage / data.len
    data = data.sort_values(by='coverage', ascending=False)
    data = data[['endpos', 'rname', 'coverage', 'Species']]
    data = data.query(f'coverage > {threshold}')
    return data

data = get_data_table(coverage_file)
name = os.path.splitext(os.path.basename(coverage_file))[0].replace('_coverage', '')
bam_folder = 'DATA/external_rna_not_paired/unclassified_sorted_bam/'

if not os.path.exists(output_file):
    with open(output_file, 'w') as f:
        pass

for endposes, rnames, spec in zip(data.endpos, data.rname, data.Species):
    for endpos, rname in zip (endposes, rnames):
        with open(output_file, 'a') as handle:
            handle.write(f"{spec.replace(' ', '_')}, {rname}, {endpos}\n")
            