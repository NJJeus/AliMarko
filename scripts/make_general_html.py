import base64
import pandas as pd
import os
import sys
sys.path.append('scripts')
import styles
import numpy as np
import glob
import argparse
import matplotlib.cm as cm
import matplotlib.colors as mcolors


def if_condition(x, message):
    if not x:
        print(message)
        exit()
        
description = """ Script that concatenate an ictv coverage data with drawings of coverage"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-c', '--coverage', type=str, help='A file fith ictv coverage')
parser.add_argument('-e', '--coverage_heatmap', type=str, help='A file fith ictv heatmap coverage')
parser.add_argument('-m', '--hmm_report', type=str, help='A hmm report file')
parser.add_argument('-p', '--hmm_heatmap', type=str, help='A file fith hmm heatmap coverage')
parser.add_argument('-o', '--output', type=str, help='An output file')



args = parser.parse_args()

if args.coverage:
    if_condition(os.path.isfile(args.coverage), "An coverage_report input file does not exists")
    ictv_coverage_file = args.coverage
else:
    if_condition(False, 'Missed -c argument')   

if args.coverage_heatmap:
    if_condition(os.path.isfile(args.coverage_heatmap), "An coverage_report input file does not exists")
    coverage_heatmap = args.coverage_heatmap
else:
    if_condition(False, 'Missed -e argument')   

if args.hmm_heatmap:
    if_condition(os.path.isfile(args.hmm_heatmap), "An coverage_report input file does not exists")
    hmm_heatmap = args.hmm_heatmap
else:
    if_condition(False, 'Missed -p argument')   
    
if args.hmm_report:
    if_condition(os.path.isfile(args.hmm_report), "An hmm_report input file does not exists")
    hmm_report = args.hmm_report
else:
    if_condition(False, 'Missed -m argument')   
    
if args.output:
    output = args.output
else:
    if_condition(False, 'Missed -o argument')


    
sample_name = os.path.splitext(os.path.basename(ictv_coverage_file))[0]
    
ictv_coverage = pd.read_csv(ictv_coverage_file)
hmm_report = pd.read_csv(hmm_report)




style = styles.set_style

introduction_header = list(ictv_coverage.columns)
introduction_table = ictv_coverage.fillna(0).to_numpy()
max_coverage, min_coverage = introduction_table[:, 1:].astype('float').max().max(), introduction_table[:, 1:].astype('float').min().min()


hmm_header = list(hmm_report.columns)
hmm_table = hmm_report.fillna(0).to_numpy()
max_hmm, min_hmm = hmm_table[:, 1:].astype('float').astype('float').max().max(), hmm_table[:, 1:].astype('float').min().min()



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





style = styles.set_style
palette = create_palette(min_coverage, max_coverage)
table_html =  styles.Table(introduction_table, introduction_header, palette, get_color).make_table()

palette = create_palette(min_hmm, max_hmm)
table_hmm =  styles.Table(hmm_table, hmm_header, palette, get_color).make_table()



with open(coverage_heatmap, "rb") as image_file:
    encoded_image = base64.b64encode(image_file.read()).decode()
    coverage_image = f'<div style="overflow: hidden"><img alt="" src="data:image/jpeg;base64,{encoded_image}" alt="Ooops! This should have been a picture" style="width: 60%; border: 2px solid #959494; min-width: 700px;"/></div>'
    
    
with open(hmm_heatmap, "rb") as image_file:
    encoded_image = base64.b64encode(image_file.read()).decode()
    hmm_image = f'<div style="overflow: hidden"><img alt="" src="data:image/jpeg;base64,{encoded_image}" alt="Ooops! This should have been a picture" style="width: 60%; border: 2px solid #959494; min-width: 700px;"/></div>'

greatings = f"""
<p>\n</p>
<h1> <center>AliMarko multisample report</h1>

""".format(name=os.path.splitext(os.path.basename(ictv_coverage_file))[0])

out = style.replace('Sample Name', sample_name) + greatings + '<p>\n</p>'  + coverage_image + '<p>\n</p>' + table_html + '<p>\n</p>' + hmm_image + '<p>\n</p>' + table_hmm + '<p>\n</p>' 

with open(output, 'w') as f:
    f.write(out)