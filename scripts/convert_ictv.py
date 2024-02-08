import argparse
import os
import pandas as pd
import numpy as np

def if_condition(x, message):
    if not x:
        print(message)
        exit()

description = """ Script that concatenate a samtools coverage output and a samtools view quality with an ictv metadata"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-c', '--coverage', type=str, help='A file with samtools coverage output')
parser.add_argument('-o', '--output', type=str, help='An output file')
parser.add_argument('-m', '--tmp_output', type=str, help='An output file with data for the plot_coverage script')
parser.add_argument('-t', '--tables_folder', type=str, help='A folder with tables ictv_taxo and genbank accessions')
parser.add_argument('-s', '--snps_file', type=str, help='A file with rname, number of snp, number of deep sites as columns and virus samples as rows')

args = parser.parse_args()




if args.coverage:
    if_condition(os.path.isfile(args.coverage), "coverage file doesn't exist")
    coverage_file = args.coverage
else:
    if_condition(False, 'Missed -c argument')

    
if args.tmp_output:
    tmp_output = args.tmp_output
else:
    if_condition(False, 'Missed -m argument')
    
if args.output:
    output_file = args.output
else:
    if_condition(False, 'Missed -o argument')

if args.tables_folder:
    if_condition(os.path.isdir(args.tables_folder), "tables folder doesn't exist")
    tables_folder = args.tables_folder
else:
    if_condition(False, 'Missed -t argument')

if args.snps_file:
    if_condition(os.path.isfile(args.snps_file), "snps file doesn't exist")
    snps_file = args.snps_file
else:
    if_condition(False, 'Missed -s argument')    

    
name = os.path.splitext(os.path.basename(coverage_file))[0]


# Read samtools coverage output and drop excess columns
data = pd.read_table(coverage_file).rename(columns={'#rname': 'rname', 'meanmapq':'meanmapq'}).drop(['numreads', 'covbases', 'meanbaseq'], axis=1)

# read an snps file
snps = pd.read_csv(snps_file)

ictv_data = pd.read_excel(f'{tables_folder}/ictv_taxo.xlsx').dropna(subset='Virus GENBANK accession')
ictv_data['genbank_list'] = ictv_data['Virus GENBANK accession'].apply(lambda i :
                                                                       [el.split(":")[-1] for el in i.replace(' ', '').split(';')])

ictv_data['Isolate_id'] = ictv_data['Species'] +  '{' + ictv_data['Sort'].astype('str')  + '}'
indeces = []
accessions = []

for index, item in ictv_data.iterrows():
    for accession in item['genbank_list']:
        indeces.append(item.Isolate_id)
        accessions.append(accession)
taxo_index = pd.DataFrame({'ictv_taxo_index':indeces}, index=accessions)

data = pd.merge(left=data, right=snps, left_on='rname', right_on='rname', how='left')
data['nucleotide_similarity'] = 1 - data['snps'] / (data['deep_sites'] + 1)

# rname contains database name, accession number and version. We need only an accession number
data['tmp_rname'] = data['rname'].copy()
data['rname'] = data['rname'].apply(lambda x: x.split('|')[1])
data = data.set_index('rname')

# Concatenate a coverage data with an ictv metadata 

data = data.copy()
data = pd.merge(left=data.reset_index(), how='left', right=taxo_index.reset_index(), left_on='rname', right_on='index') # Get index of all accession number
data = pd.merge(left=data, right=ictv_data.reset_index(), left_on='ictv_taxo_index', right_on='Isolate_id') # Get ictv metadata by accession number

# Generalize fragments of one virus
data[['coverage', 'meanmapq']] = data[['coverage', 'meanmapq']]  * 0.01
data['len'] = data.endpos - data.startpos + 1 
data[['weighted_coverage', 'weighted_quality', 'weighted_meandepth']] = (data[['coverage', 'meanmapq', 'meandepth']].T * data['len']).T

per_fragment_columns = ['len', 'coverage', 'meanmapq', 'meandepth', 'deep_sites', 'nucleotide_similarity', 'snps']
data[[f'fragments_{c}' for c in per_fragment_columns]] = data[per_fragment_columns].copy()


# create dict for further grouping
agg_all_columns_dict = dict()
## summing
agg_all_columns_dict.update({i:'sum' for i in ['len', 'snps', 'deep_sites', 'weighted_coverage', 'weighted_quality', 'weighted_meandepth']})
## Generate list with values for all fragment
agg_all_columns_dict.update({i:lambda x: [round(c, 4) for c in list(x)] for i in [f'fragments_{c}' for c in per_fragment_columns]})
agg_all_columns_dict.update({'tmp_rname': lambda x: list(x)})
agg_all_columns_dict.update({i:'first' for i in ['rname', 'Realm', 'Kingdom',
       'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus isolate designation',
              'Virus name(s)', 'Virus GENBANK accession', 'Host source', 'genbank_list']})

data = data.groupby('Isolate_id').agg(agg_all_columns_dict)


data['coverage'] = (data.weighted_coverage / data.len).round(2)
data['meanmapq'] = (data.weighted_quality / data.len * 100).round(2)
data['meandepth'] = (data.weighted_meandepth / data.len).round(2)
data = data.sort_values(by='coverage', ascending=False)
data['nucleotide_similarity'] = 1 - data['snps'].div(data['deep_sites'].replace(0, np.inf))
data['deep_sites'] = data['deep_sites'].astype('int')

indeces = []
tmp_names = []
fragments_len = []

for index, item in data.query('coverage > 0.05 & meanmapq > 19').iterrows():
    for i in range(len(item['tmp_rname'])):
        indeces.append(index)
        tmp_names.append(item.tmp_rname[i])
        fragments_len.append(item.fragments_len[i])
tmp_pic_data = pd.DataFrame(index=indeces, data={'tmp_names':tmp_names, 'fragments_len':fragments_len})
tmp_pic_data.to_csv(tmp_output, header=False)



# Clean and output
data = data[['Virus name(s)', 'Host source', 'len', 'deep_sites', 'coverage', 'meandepth',  'meanmapq', 'snps', 'nucleotide_similarity',
             *[f'fragments_{c}' for c in per_fragment_columns], 'Virus GENBANK accession', 'Realm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class',
   'Order', 'Family', 'Genus', 'Species']]

print(output_file)

data.to_csv(output_file)