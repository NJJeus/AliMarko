
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()
import os

import argparse

import colour
red = colour.Color('#2B65EC')
colors = list(red.range_to(colour.Color('red'), 21))


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
    plt.figure(dpi=600)
    contig_data = data.query(f'Name == "{contig}"')
    length_contig = contig_data.Length_contig.max()
    max_score = np.log10(contig_data.Score.max())
    min_score = np.log10(contig_data.Score.min())
    
    i=0
    for index, row in contig_data.iterrows():
        i-=4
        try:
            try:
                color = str(colors[int(((np.log10(row.Score)-min_score)/(max_score-min_score) * 20))])
            except Exception:
                color = str(colors[0])
            
            plt.arrow(row.From, i, row.To-row.From, 0, 
             color=color,
            width=0.2, head_width=0.4, length_includes_head=True, head_length=abs(row.To-row.From)/10)
            
            plt.text((row.From + row.To)/2, i+1, f'{row.Query}:{row["Positive terms"]}', fontdict={'size':6}, ha='center')
        except Exception:
            print(f'error with {contig}, {row.Query}')
            pass
        
        

    plt.yticks([])

    import matplotlib.patches as mpatches
    handles = []
    for col, lab in zip(colors[::-4], range(20, -4, -4)):
        l = int(lab * ((contig_data.Score.max()/20))) // 100 * 100
        handles.append(mpatches.Patch(color=str(col), label=l))
    
    plt.plot([0, length_contig], [0, 0],linewidth=4.0)
    plt.text(length_contig/2, 1, contig, ha='center')
    
    plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), title='Score')
    plt.ylim([i-2, 0.2])
    plt.xlim([0, length_contig]);
    plt.savefig(f'{output_dir}/{contig}.png')