import base64
import pandas as pd
import os
import sys
sys.path.append('scripts')
import styles
import numpy as np
import glob
import argparse
import matplotlib

import matplotlib.cm as cm
import matplotlib.colors as mcolors


def if_condition(x, message):
    if not x:
        print(message)
        exit()
        
def create_palette(min_score, max_score):
    return None


def get_color(palette, value):
    try:
        color = mcolors.to_hex(palette.to_rgba(float(value)),keep_alpha=False)
    except Exception:
        return 'white'
    return color
        
description = """ Script that concatenate an ictv coverage data with drawings of coverage"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-c', '--coverage', type=str, help='A file fith ictv coverage')
parser.add_argument('-o', '--output', type=str, help='An output file')
parser.add_argument('-d', '--drawings', type=str, help='A folder with drawings')
parser.add_argument('-a', '--hmm_drawings', type=str, help='A folder with hmm drawings')
parser.add_argument('-m', '--hmm_report', type=str, help='A hmm report file')
parser.add_argument('-t', '--trees_folder', type=str, help='A folder with trees ')


args = parser.parse_args()

palette = create_palette(0, 10)

if args.coverage:
    if_condition(os.path.isfile(args.coverage), "An input file does not exists")
    ictv_coverage_file = args.coverage
else:
    if_condition(False, 'Missed -c argument')   
if args.output:
    output = args.output
else:
    if_condition(False, 'Missed -o argument')
if args.drawings:
    if_condition(os.path.isdir(args.drawings), "drawings folder doesn't exist")
    coverage_draw_folder_path = args.drawings
else:
    if_condition(False, 'Missed -t argument')

if args.hmm_drawings:
    if_condition(os.path.isdir(args.hmm_drawings), "hmm drawings folder doesn't exist")
    hmm_draw_folder_path = args.hmm_drawings
else:
    if_condition(False, 'Missed -a argument')    
    
if args.hmm_report:
    if_condition(os.path.isfile(args.hmm_report), "A hmm report file does not exists")
    hmm_report = args.hmm_report
else:
    if_condition(False, 'Missed -m argument') 

if args.trees_folder:
    if_condition(os.path.isdir(args.trees_folder), "drawings folder doesn't exist")
    trees_folder = args.trees_folder
else:
    if_condition(False, 'Missed -t argument')


sample_name = os.path.splitext(os.path.basename(ictv_coverage_file))[0]
ictv_coverage = pd.read_csv(ictv_coverage_file)
tree_files = glob.glob(f'{trees_folder}/*')


def create_palette(min_score, max_score):
    if max_score != min_score:
        pass  # Set the values for the colorbar
    else:
        min_score, max_score = 0, min_score
    cmap = cm.ScalarMappable(cmap='Blues') 
    cmap.set_array([min_score, max_score*1.5]) 
    cmap.autoscale()
    return cmap


def get_color(palette, value):
    try:
        color = mcolors.to_hex(palette.to_rgba(float(value)),keep_alpha=False)
        print(f'{color} {value}')
    except Exception:
        print(f'error {value}')
        return 'white'
    return color

def add_string(x, string_to_add):
    if pd.notna(x):
        return  string_to_add + "Sequences of this virus have association with " + str(x) + "." 
    else:
        return x


ictv_coverage['Feature'] = ictv_coverage['Feature'].apply(add_string, string_to_add='CONTAMINATION_INFO:')
ictv_coverage['Species'] = ictv_coverage['Species'] + ictv_coverage['Feature'].fillna('')

introduction_frame = ictv_coverage[['Species', 'Host source',
               'coverage', 'meandepth', 'Genus', 'Family', 'Realm']].query('coverage != 0')



introduction_header = ['Species', 'Host source', 'Coverage width', 'Mean depth', 'Genus', 'Family', 'Realm']
introduction_table = introduction_frame.to_numpy()


hmm_frame = pd.read_csv(hmm_report)[['Query', 'Taxon', 'Name', 'Positive terms', 'Score', 'Threshold', 'From', 'To', 'Score_ratio']]
hmm_frame[['Score', 'Threshold', 'Score_ratio']] = hmm_frame[['Score', 'Threshold', 'Score_ratio']].astype('int')
hmm_frame = hmm_frame.rename(columns={'Positive terms': 'Putative Protein', 'Query': 'HMM'})
hmm_header = hmm_frame.drop(['Threshold', 'From', 'To', 'Score_ratio'], axis=1).columns
hmm_table = hmm_frame.drop(['Threshold', 'From', 'To', 'Score_ratio'], axis=1).to_numpy()




ictv_drawings = ictv_coverage.query('coverage > 0.05').set_index('Host source')

ictv_drawings['genbank_list'] = ictv_drawings['Virus GENBANK accession'].apply(lambda i :
                                                                       [el.split(":")[-1] for el in i.replace(' ', '').split(';')])


ictv_drawings[['fragments_len', 'fragments_meandepth', 'fragments_coverage', 'fragments_nucleotide_similarity', 'fragments_meanmapq', 'fragments_snps']] = ictv_drawings[['fragments_len',
       'fragments_meandepth', 'fragments_coverage', 'fragments_nucleotide_similarity', 'fragments_meanmapq', 'fragments_snps']].applymap(lambda x: x.strip('][').split(', '))

ictv_drawings = ictv_drawings.groupby(level=0).agg({i:lambda x: list(x) for i in ictv_drawings.columns})




host_dict = {}

# Coverage images
for host, row in ictv_drawings.iterrows():
    viruses = {}
    for virus_loc in range(len(row.Isolate_id)):
        virus_tech_name = row['Isolate_id'][virus_loc]
        virus_name = row['Species'][virus_loc]
        if 'CONTAMINATION_INFO:' in virus_name:
            virus_name = f'<h3 class="tooltip" style="color:#C80000">{virus_name.split("CONTAMINATION_INFO:")[0]} <span class="tooltip-text">{virus_name.split("CONTAMINATION_INFO:")[1]}</span></h3>'
        else:
            virus_name = f'<h3>{virus_name}</h3>'
        
        virus_list = np.array([row.genbank_list[virus_loc], row.fragments_len[virus_loc], row.fragments_coverage[virus_loc], 
                               row.fragments_nucleotide_similarity[virus_loc], row.fragments_meandepth[virus_loc], row.fragments_meanmapq[virus_loc], row.fragments_snps[virus_loc]]).T


        images = []
        for i in row.genbank_list[virus_loc]:
            try:
                picture_path = glob.glob(f"{coverage_draw_folder_path}/{virus_tech_name}/*{i}*")[0]
            except Exception:
                continue
            with open(picture_path, "rb") as image_file:
                encoded_image = base64.b64encode(image_file.read()).decode()
                images.append(f'<div style="overflow: hidden; max-height:700px;"><img alt="" src="data:image/jpeg;base64,{encoded_image}" alt="Ooops! This should have been a picture" style="width: 60%; border: 2px solid #959494; min-width: 700px;"/></div>')
        images = '\n'.join(images)

        viruses.update({f'{virus_name}':styles.Table(virus_list, ['Fragment', 'Len', 'Coverage width', 'Nucleotide similarity', 'Mean Depth', 'MeanMAPQ',' SNP Count'], palette, get_color).make_table() + images})

    host_dict.update({f'<h2>Host: {host}</h2>':styles.Details(viruses).make_details()})


order = [f'<h2>Host: {host}</h2>' for host in ['vertebrates', 'marine (S)', 'invertebrates',  'invertebrates, vertebrates', 'plants, invertebrates', 'fungi', 'plants', 'protists (S)', 'algae', 'protists', 'bacteria', 'archaea', 'sewage (S)']]
order = list(set(order + [f'<h2>{i}</h2>' for i in host_dict.keys()]))
def key_func(x):
    return order.index(x)
host_dict = sorted(host_dict.items(), key=lambda x: key_func(x[0]))
host_dict = {f"{key}":value for key, value in host_dict}


positive_contigs = hmm_frame.query('Score_ratio > 2').groupby('Name').count().query('Score > 0').reset_index().Name
positive_contigs = positive_contigs[positive_contigs.isin(positive_contigs)]

list_of_contigs = hmm_frame['Name'][hmm_frame['Name'].isin(positive_contigs)].unique()
contigs_dict = {}
images_hmm = []
table_contigs = []
for contig in list_of_contigs:
    
    table_contig = hmm_frame.query(f'Name == "{contig}"').drop(['Threshold', 'From', 'To', 'Score_ratio'], axis=1)
    models = hmm_frame.query(f'Name == "{contig}"').Taxon.unique().tolist()
    models = ",".join(models) if len(models) < 3 else f"({len(models)} taxa)"
    
    try:
        picture_path = glob.glob(f"{hmm_draw_folder_path}/*{contig}*")[0]
    except Exception:
        continue
    
    trees_pictures = ''
    for tree in [i for i in tree_files if contig in i]:
        print(tree)
        with open(tree, "rb") as tree_image_file:
            tree_encoded_image= base64.b64encode(tree_image_file.read()).decode()
            tree_HMM = tree.split('__')[1]
            protein_HMM = table_contig.query(f'HMM == "{tree_HMM}"')['Putative Protein'].unique()[0]
            tree_name = f"<h3>Phylogenetic Tree of {tree_HMM}-Matched Amino Acid Sequences</h3><h3>Putative protein: {protein_HMM}</h3>"
            html_tree = f'<div style="overflow: hidden; max-height:700px;"><img alt="" src="data:image/jpeg;base64,{tree_encoded_image}" alt="Ooops! This should have been a picture" style="width: 60%; border: 2px solid #959494; min-width: 700px;"/></div><p>\n</p>'
            trees_pictures = trees_pictures + '\n' + tree_name + html_tree       
        
        
    with open(picture_path, "rb") as image_file:
        encoded_image = base64.b64encode(image_file.read()).decode()
        html_image = f'<div style="overflow: hidden; max-height:700px;"><img alt="" src="data:image/svg+xml;charset=utf-8;base64,{encoded_image}" alt="Ooops! This should have been a picture" style="width: 60%; border: 2px solid #959494; min-width: 700px;"/></div><p>\n</p>'
        details_image = styles.Details({f'<h3>{contig}:{models}<h3>':styles.Table(table_contig.to_numpy(), table_contig.columns, palette, get_color).make_table()+'\n'+ html_image + '\n' + trees_pictures}).make_details()
        images_hmm.append(details_image)
images_hmm = '\n'.join(images_hmm)
            
images_hmm = styles.Details({'<h2>Contigs</h2>\n': images_hmm}).make_details()

    
style = styles.set_style

greatings = f"""
<p>\n</p>
<h1> <center>Sample: {sample_name} </h1>

""".format(name=os.path.splitext(os.path.basename(ictv_coverage_file))[0])

table_html =  styles.Table(introduction_table, introduction_header, palette, get_color).make_table()

hmm_greetings = f"<p>\n</p><h2> <center> HMM Hit Summary </h2><p>\n</p>"

table_hmm = styles.Table(hmm_table, hmm_header, palette, get_color).make_table()

details = styles.Details(host_dict).make_details()

out = style.replace('Sample Name', sample_name) + greatings + "<h2>Mapping Summary</h2>" + table_html + hmm_greetings +   table_hmm + '<p>\n</p><hr>' + "<h2>Mapping Details</h2>" + details + '<p>\n</p><p>\n</p><hr>' + "<h2>HMM Module Results by Contig</h2>" + images_hmm + '<p>\n</p>'

with open(output, 'w') as f:
    f.write(out)
