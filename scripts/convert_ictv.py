import argparse
import os
import pandas as pd

def if_condition(x, message):
    if not x:
        print(message)
        exit()

description = """ Script that concatenate a samtools coverage output and a samtools view quality with an ictv metadata"""
        
parser = argparse.ArgumentParser(description=description)

parser.add_argument('-c', '--coverage', type=str, help='A file with samtools coverage output')
parser.add_argument('-q', '--quality', type=str, help='A file with rname and quality columns from samtools view output')
parser.add_argument('-o', '--output', type=str, help='An output file')
parser.add_argument('-t', '--tables_folder', type=str, help='A folder with tables ictv_taxo and genbank accessions')

args = parser.parse_args()

if args.coverage:
    if_condition(os.path.isfile(args.coverage), "coverage file doesn't exist")
    coverage_file = args.coverage
else:
    if_condition(False, 'Missed -c argument')

if args.quality:
    if_condition(os.path.isfile(args.quality), "quality file doesn't exist")
    quality_file = args.quality
else:
    if_condition(False, 'Missed -q argument')    
    
if args.output:
    output_file = args.output
else:
    if_condition(False, 'Missed -o argument')
if args.tables_folder:
    if_condition(os.path.isdir(args.tables_folder), "tables folder doesn't exist")
    tables_folder = args.tables_folder
else:
    if_condition(False, 'Missed -t argument') 

ictv_data = pd.read_excel(f'{tables_folder}/ictv_taxo.xlsx')
taxo_index = pd.read_csv(f'{tables_folder}/genbank_accessions.csv', index_col=0)

# Rad quality file
quality = pd.read_csv(quality_file, header=None, names=['rname', 'quality'], sep='\t')
quality.rname = quality.rname.apply(lambda x: x.split('|')[1])
quality = quality.groupby('rname').mean()

# Read samtools coverage output and drop excess columns
data = pd.read_table(coverage_file).rename(columns={'#rname': 'rname'}).drop(['numreads', 'covbases',
       'meandepth', 'meanbaseq', 'meanmapq'], axis=1)

# rname contains database name, accession number and version. We need only an accession number
data['rname'] = data['rname'].apply(lambda x: x.split('|')[1])
data = data.set_index('rname')
data['quality'] = quality.quality

name = os.path.splitext(os.path.basename(coverage_file))[0]

# Concatenate a coverage data with an ictv metadata 
data = data.copy()
data = pd.merge(left=data.reset_index(), right=taxo_index.reset_index(), left_on='rname', right_on='index') # Get index of all accession number
data = pd.merge(left=data, right=ictv_data.reset_index(), left_on='ictv_taxo_index', right_on='index') # Get ictv metadata by accession number


# Generalize fragments of one virus
data[['coverage', 'quality']] = data[['coverage', 'quality']]  * 0.01
data['len'] = data.endpos - data.startpos
data[['weighted_coverage', 'weighted_quality']] = (data[['coverage', 'quality']].T * data['len']).T


agg_all_columns_dict = dict()
agg_all_columns_dict.update({'len':'sum', 'weighted_coverage':'sum', 'weighted_quality':'sum'})
agg_all_columns_dict.update({'coverage': lambda x: [i for i in list(x)]})
agg_all_columns_dict.update({i:'first' for i in ['rname', 'Realm', 'Kingdom',
       'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Order', 'Family', 'Genus', 'Species',
              'Virus name(s)', 'Virus GENBANK accession', 'Host source']})
data = data.groupby('Species').agg(agg_all_columns_dict)
data['fragments_coverage'] = data['coverage'].copy()


data['coverage'] = data.weighted_coverage / data.len
data['quality'] = data.weighted_quality / data.len
data = data.sort_values(by='coverage', ascending=False)

# Clean and output
data = data[['Virus name(s)', 'Host source', 'coverage', 'quality', 'len', 'fragments_coverage', 'Realm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class',
   'Order', 'Family', 'Genus', 'Species',
   'Virus GENBANK accession']]

data.to_csv(output_file)

