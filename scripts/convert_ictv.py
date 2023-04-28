import argparse
import os
import pandas as pd

def if_condition(x, message):
    if not x:
        print(message)
        exit()

parser = argparse.ArgumentParser(description='Description of your script')

parser.add_argument('-c', '--coverage', type=str, help='A file with samtools coverage output')
parser.add_argument('-o', '--output', type=str, help='An output file')
parser.add_argument('-t', '--tables_folder', type=str, help='A folder with tables ictv_taxo and genbank accessions')

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

ictv_data = pd.read_excel(f'{tables_folder}/ictv_taxo.xlsx')
taxo_index = pd.read_csv(f'{tables_folder}/genbank_accessions.csv', index_col=0)


agg_all_columns_dict = dict()
agg_all_columns_dict.update({'len':'sum', 'weighted_coverage':'sum'})
agg_all_columns_dict.update({i:'first' for i in ['rname', 'Realm', 'Kingdom',
       'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Order', 'Family', 'Genus', 'Species',
              'Virus name(s)', 'Virus GENBANK accession', 'Host source']})    

data = pd.read_table(coverage_file).rename(columns={'#rname': 'rname'}).drop(['startpos', 'numreads', 'covbases',
       'meandepth', 'meanbaseq', 'meanmapq'], axis=1)
data['rname'] = data['rname'].apply(lambda x: x.split('|')[1])

name = os.path.splitext(os.path.basename(coverage_file))[0]

data = data.copy()
data = data.set_index('rname')
data = pd.merge(left=data.reset_index(), right=taxo_index.reset_index(), left_on='rname', right_on='index')
data = pd.merge(left=data, right=ictv_data.reset_index(), left_on='ictv_taxo_index', right_on='index')


data['weighted_coverage'] = data.endpos * 0.01 * data.coverage
data['len'] = data.endpos

data = data.groupby('Species').agg(agg_all_columns_dict)

data['coverage'] = data.weighted_coverage / data.len
data = data.sort_values(by='coverage', ascending=False)

data = data[['Virus name(s)', 'Host source', 'coverage','Realm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class',
   'Order', 'Family', 'Genus', 'Species',
   'Virus GENBANK accession']]

data.to_csv(output_file)

