
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
sns.set_theme()
import os

import argparse

cmap = cm.ScalarMappable(cmap='jet')  # Choose the colormap


description = """The script gets a report csv hmm file and returns a folder with plots fro all contings"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-i', '--input_file', type=str, help='A file or directory with fastq/fasta extenstion')
parser.add_argument('-o', '--output_dir', type=str, help='An output folder')


args = parser.parse_args()
print(args)

def if_condition(x, message):
    if not x:
        print(message)
        exit()
        
if args.input_file:
    if os.path.isfile(args.input_file):
        input_file = args.input_file
    else:
        if_condition(True, "An input file does not exists")
else:
    if_condition(False, 'Missed -f argument')
    
if args.output_dir:
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)    
else:
    if_condition(False, 'Missed -o argument')


    
data = pd.read_csv(input_file, index_col=0)

for contig in data.Name.unique():
    plt.clf()
    fig = plt.figure(dpi=400)
    ax = fig.add_subplot(1,1,1)
    contig_data = data.query(f'Name == "{contig}"')
    length_contig = contig_data.Length_contig.max()
    max_score = contig_data.Score.max()
    min_score = contig_data.Score.min()
    if max_score != min_score:
        pass  # Set the values for the colorbar
    else:
        min_score, max_score = 0, min_score
    cmap = cm.ScalarMappable(cmap='jet') 
    cmap.set_array([min_score, max_score]) 
    cmap.autoscale()
    i=0
    for index, row in contig_data.iterrows():
        i-=4
        color = cmap.to_rgba(row.Score)

        ax.arrow(row.From, i, row.To-row.From, 0, 
         color=color,
        width=0.2, head_width=0.4, length_includes_head=True, 
                 head_length=max(abs(row.To-row.From)/10, length_contig/100))
        arrowtext = f'{row.Query}:{row["Positive terms"]}, {row.Taxon}'
        ax.text((row.From + row.To)/2, i+1, arrowtext, fontdict={'size':6}, ha='center')
        
        

    ax.set_yticks([])
    cbar = plt.colorbar(cmap, orientation='vertical', ticks=np.linspace(min_score, max_score, num=5), ax=ax)
    cbar.set_label('Score')
        

    
    ax.plot([0, length_contig], [0, 0],linewidth=4.0)
    ax.text(length_contig/2, 1, contig, ha='center')
    
    
    ax.set_ylim([i-2, 0.2])
    ax.set_xlim([0, length_contig]);
    plt.savefig(f'{output_dir}/{contig}.png')
    plt.close()