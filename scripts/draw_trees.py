from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import glob
import argparse
import os


def if_condition(x, message):
    if not x:
        print(message)
        exit()
def label_func(x):
    return str(x)[:30]

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
    try:
        tree = Phylo.read(file, "newick")
    except Exception:
        continue
    tree.root_at_midpoint()
    leaf_nodes = tree.get_terminals()
    fams = {}
    # Print the names of the leaf nodes
    for leaf in leaf_nodes:
        fams[leaf.name] = leaf.name.split('|')[-1]
    palete_set = set(fams.values())
    legend_dict = {}
    c=0
    for i in palete_set:
        legend_dict[i] = sns.color_palette(n_colors=len(palete_set))[c]
        c+=1
    labels_dict = {k[:30]:legend_dict[v] for k, v in fams.items()}
    legend_items = list(legend_dict.items())
    patches = [mpatches.Patch(color=color, label=label) for label, color in legend_items]
    
    fig, axes = plt.subplots(1, 1, figsize=(16, 8), dpi=100)
    # Draw the tree
    axes.legend(handles=patches, bbox_to_anchor=(1, 1), loc='upper left')
    Phylo.draw(tree, axes=axes, label_func=label_func,
               label_colors=labels_dict)

    axes.set_title(f"{name}")
    axes.set_ylabel('')
    fig.savefig(f"{output_folder}/{name}_tree.svg", format='svg')
