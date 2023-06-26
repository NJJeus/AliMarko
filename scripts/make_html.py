import base64
import pandas as pd
import os
import sys
sys.path.append('scripts')
import styles
import numpy as np
import glob
import argparse


def if_condition(x, message):
    if not x:
        print(message)
        exit()
        
description = """ Script that concatenate an ictv coverage data with drawings of coverage"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-c', '--coverage', type=str, help='A file fith ictv coverage')
parser.add_argument('-o', '--output', type=str, help='An output file')
parser.add_argument('-d', '--drawings', type=str, help='A folder with drawings')
args = parser.parse_args()



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
    drawings_folder_path = args.drawings
else:
    if_condition(False, 'Missed -t argument')

sample_name = os.path.splitext(os.path.basename(ictv_coverage_file))[0]
ictv_coverage = pd.read_csv(ictv_coverage_file)

introduction_frame = ictv_coverage[['Virus name(s)', 'Host source',
               'coverage', 'meandepth', 'Genus', 'Family', 'Realm']][:50]
introduction_header = ['Virus name(s)', 'Host source', 'Coverage width', 'Mean depth', 'Genus', 'Family', 'Realm']
introduction_table = introduction_frame.to_numpy()

ictv_drawings = ictv_coverage.query('coverage > 0.1').set_index('Host source')

ictv_drawings['genbank_list'] = ictv_drawings['Virus GENBANK accession'].apply(lambda i :
                                                                       [el.split(":")[-1] for el in i.replace(' ', '').split(';')])
ictv_drawings[['fragments_len', 'fragments_meandepth', 'fragments_coverage', 'fragments_nucleotide_similarity', 'fragments_meanmapq', 'fragments_snps']] = ictv_drawings[['fragments_len',
       'fragments_meandepth', 'fragments_coverage', 'fragments_nucleotide_similarity', 'fragments_meanmapq', 'fragments_snps']].applymap(lambda x: x.strip('][').split(', '))

ictv_drawings = ictv_drawings.groupby(level=0).agg({i:lambda x: list(x) for i in ictv_drawings.columns})

host_dict = {}
for host, row in ictv_drawings.iterrows():
    viruses = {}
    for virus_loc in range(len(row.Isolate_id)):
        virus_name = row.Isolate_id[virus_loc]
        virus_list = np.array([row.genbank_list[virus_loc], row.fragments_len[virus_loc], row.fragments_coverage[virus_loc], 
                               row.fragments_nucleotide_similarity[virus_loc], row.fragments_meandepth[virus_loc], row.fragments_meanmapq[virus_loc], row.fragments_snps[virus_loc]]).T


        images = []
        for i in row.genbank_list[virus_loc]:
            try:
                picture_path = glob.glob(f"{drawings_folder_path}/{virus_name}/*{i}*")[0]
            except Exception:
                continue
            with open(picture_path, "rb") as image_file:
                encoded_image = base64.b64encode(image_file.read()).decode()
                images.append(f'<div style="overflow: hidden; max-height:700px;"><img alt="" src="data:image/jpeg;base64,{encoded_image}" alt="Ooops! This should have been a picture" style="width: 60%; border: 2px solid #959494; min-width: 700px;"/></div>')
        images = '\n'.join(images)

        viruses.update({f'<h3>{virus_name}</h3>':styles.Table(virus_list, ['Fragment', 'Len', 'Coverage width', 'Nucleotide similarity', 'Mean Depth', 'MeanMAPQ',' SNP Count']).make_table() + images})

    host_dict.update({f'<h2>{host}</h2>':styles.Details(viruses).make_details()})


order = [f'<h2>{host}</h2>' for host in ['vertebrates', 'invertebrates', 'fungi', 'plants', 'algae', 'protists', 'bacteria', 'archaea']]
def key_func(x):
    return order.index(x)
host_dict = sorted(host_dict.items(), key=lambda x: key_func(x[0]))
host_dict = {key:value for key, value in host_dict}

    
style = styles.set_style

greatings = f"""
<p>\n</p>
<h1> <center> {sample_name} </h1>

""".format(name=os.path.splitext(os.path.basename(ictv_coverage_file))[0])

table_html =  styles.Table(introduction_table, introduction_header).make_table()

details = styles.Details(host_dict).make_details()

out = style.replace('Sample Name', sample_name) + greatings + table_html + details

with open(output, 'w') as f:
    f.write(out)