import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import numpy as np
import matplotlib.cm as cm
import os

description = """The script gets a fastq or a fasta file, translate it and analyzes sequences with HMM profiles"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-i', '--cov_dir', type=str, help='Directory with coverage files')
parser.add_argument('-c', '--output_coverage_table', type=str, help='Output coverage general file')
parser.add_argument('-p', '--output_coverage_pic', type=str, help='Output coverage plot file')


parser.add_argument('-m', '--hmm_dir', type=str, help='Directory with hmm files')
parser.add_argument('-t', '--output_hmm_table', type=str, help='Output hmm general file')
parser.add_argument('-u', '--output_hmm_pic', type=str, help='Output hmm plot file')

def filter_queries(group, column='Taxon'):
    first_query = group[column].iloc[0]
    return group[group[column] == first_query]

def if_condition(x, message):
    if not x:
        print(message)
        exit()

args = parser.parse_args()
print('hete')

if args.cov_dir:
    if os.path.isdir(args.cov_dir):
        cov_dir = args.cov_dir
    else:
        if_condition(False, "An input folder with coverage files does not exist")
else:
    if_condition(False, 'Missed -i argument')
    
if args.hmm_dir:
    if os.path.isdir(args.hmm_dir):
        hmm_dir = args.hmm_dir
    else:
        if_condition(False, "An input folder with hmm report files does not exist")
else:
    if_condition(False, 'Missed -m argument')


# Step 1: Get a list of all CSV files in a folder using glob
csv_files = glob.glob(os.path.join(cov_dir, '*.csv'))

# Step 2, 3, and 4: Read each CSV file into a DataFrame, set "Sample" column to the filename without extension, and concatenate all the DataFrames
dfs = []
for file in csv_files:
    df = pd.read_csv(file)
    filename = os.path.splitext(os.path.basename(file))[0]
    df['Sample'] = filename
    dfs.append(df.query('coverage > 0.05'))

combined_align = pd.concat(dfs, ignore_index=True)

combined_align = combined_align.sort_values('coverage', ascending=False)[['Isolate_id', 'coverage', 'Sample']]

combined_align = combined_align.pivot(index='Isolate_id', columns='Sample', values='coverage')

combined_align['max_value'] = combined_align.max(axis=1) + combined_align.mean(axis=1) * 0.5
combined_align = combined_align.sort_values(by='max_value', ascending=False)
combined_align.drop(columns='max_value', inplace=True)

combined_align = combined_align.fillna(0)[:200]

combined_align.to_csv(args.output_coverage_table)

plt.figure(figsize=(14, 9))
fig = sns.heatmap(combined_align, cmap=sns.color_palette("mako_r", as_cmap=True))
plt.title('References coverage. General heatmap ')
plt.subplots_adjust(left=0.35, bottom=0.25)
fig.collections[0].colorbar.set_label("Coverage width")
plt.savefig(args.output_coverage_pic, format='png')

# HMM part

# Step 1: Get a list of all CSV files in a folder using glob
csv_files = glob.glob(os.path.join(hmm_dir, '*.csv'))

# Step 2, 3, and 4: Read each CSV file into a DataFrame, set "Sample" column to the filename without extension, and concatenate all the DataFrames
dfs = []
for file in csv_files:
    df = pd.read_csv(file, index_col=0).query('Score_ratio > 1').groupby('Name').apply(filter_queries).reset_index(drop=True)
    filename = os.path.splitext(os.path.basename(file))[0]
    df['Sample'] = filename
    dfs.append(df)

combined_hmm = pd.concat(dfs, ignore_index=True)

combined_hmm.shape

combined_hmm = combined_hmm[['Taxon', 'Sample', 'Score_ratio']].groupby(['Taxon', 'Sample']).sum().reset_index().pivot(index='Taxon', columns='Sample', values='Score_ratio')

combined_hmm['max_value'] = combined_hmm.max(axis=1)
combined_hmm = combined_hmm.sort_values(by='max_value', ascending=False)
combined_hmm.drop(columns='max_value', inplace=True)

combined_hmm = combined_hmm.round(2)

combined_hmm.to_csv(args.output_hmm_table)


max_score = combined_hmm.max().max()
min_score = combined_hmm.min().min()
if max_score != min_score:
    pass  # Set the values for the colorbar
else:
    min_score, max_score = 0, min_score
cmap = cm.ScalarMappable(cmap='jet') 
cmap.set_array([min_score, max_score]) 
cmap.autoscale()


plt.figure(figsize=(16, 11), dpi=250)
plt.subplots_adjust(left=0.1, bottom=0.25, right=1)
fig = sns.heatmap(np.log(combined_hmm).fillna(0)[:35], cmap=sns.color_palette("mako_r", as_cmap=True))
plt.title('HMM general heatmap')
plt.ylabel('Taxon')
plt.yticks(fontsize=8)
fig.collections[0].colorbar.set_label("Sum of normalized score")

plt.savefig(args.output_hmm_pic, format='png')

