from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import glob
import argparse
import os
import re


def if_condition(x, message):
    if not x:
        print(message)
        exit()
        
    
    
def label_func(x):
    x = re.sub(r'{[^{}]*}', '', str(x)).split('#')[0]
    return str(x)[:30]

def round_bootstrap(clade):
    if clade.confidence:
         clade.confidence = int(clade.confidence * 10) / 10
    for child in clade.clades:
        round_bootstrap(child)
        
def bootstrap_colors(bootstrap):
    return f"hsl({bootstrap}, 100%, 50%)"

description = """A script"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-i', '--input_folder', type=str, help='')
parser.add_argument('-o', '--output_folder', type=str, help='')

args = parser.parse_args()




if args.input_folder:
    if_condition(os.path.exists(args.input_folder), "coverage file doesn't exist")
    input_folder  = args.input_folder
else:
    if_condition(False, 'Missed -i argument')

if args.output_folder:
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)    
else:
    if_condition(False, 'Missed -o argument')     

files = glob.glob(f'{input_folder}/*.treefile')

for file in files:
    name = file.split('/')[-1].replace('_tree.treefile', '')
    print(name)
    try:
        tree = Phylo.read(file, "newick")
    except Exception:
        print(f"Trere is no tree in treefile {file}")
        continue
        
    tree.root_at_midpoint()
    round_bootstrap(tree.root)
    leaf_nodes = tree.get_terminals()
    fams = {}
    # Print the names of the leaf nodes
    for leaf in leaf_nodes:
        fams[leaf.name] = leaf.name.split('|')[-1].split('#')[0]
        if "NODE" in leaf.name:
            fams[leaf.name] = "_".join(leaf.name.split('_')[:2])
    palete_set = set(fams.values())
    legend_dict = {}
    palette = sns.color_palette("husl", n_colors=len(palete_set))
    c=0
    for i in palete_set:
        legend_dict[i] = [i*0.9 for i in palette[c]]
        if i[:4] == 'NODE':
            legend_dict[i] = (0.8901960784313725, 0.10196078431372549, 0.10980392156862745)
        
        c+=1
    labels_dict = {label_func(k):legend_dict[v] for k, v in fams.items()}
    legend_items = list(legend_dict.items())
    patches = [mpatches.Patch(color=color, label=label) for label, color in legend_items]
    
    fig, axes = plt.subplots(1, 1, figsize=(11.5, 5.94), dpi=250)
    # Draw the tree
    axes.legend(handles=patches, bbox_to_anchor=(1, 1), loc='upper left')
    Phylo.draw(tree, axes=axes, label_func=label_func,
               label_colors=labels_dict)

    
    axes.set_title(f"Phylogeny of {name.split('__')[1]}-Matched AA Sequences from {name.split('__')[0].replace('test_', '')} contig")
    x_min, x_max = axes.get_xlim()
    axes.set_xlim([x_min, x_max+(x_max-x_min)*0.5])
    axes.set_ylabel('')
    axes.set_xlabel('Branch Length')
    axes.set_yticklabels([])
    axes.tick_params(axis='y', which='both', length=0)
    fig.subplots_adjust(right=0.8, left=0.05)
    fig.savefig(f"{output_folder}/{name}_tree.jpg", format='jpg')
