import argparse
import os
import pandas as pd
import glob
import numpy as np

def per_fragments_sum_coverage(x, y):
    return [[round(i_x + i_y, 2) for i_x, i_y in zip(el_x, el_y)] for el_x, el_y in zip(x, y)]

def mean_per_fragments_coverage(x, n):
    return [int(round(i/n, 2) * 100) for i in x]

description = """ Script that generalize many ictv coverage tables(output of convert_ictv.py)"""
        
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--files', nargs='+', help='list of files')
parser.add_argument('-o', '--output', type=str, help='An output tsv file')
args = parser.parse_args()

files = args.files
out = args.output

common_file = pd.read_csv(files[0], converters={'fragments_coverage': pd.eval}).sort_values(by='Species')
max_coverage = np.array(common_file.coverage)
min_coverage = np.array(common_file.coverage)
for i in files[1:]:
    file = pd.read_csv(i, converters={'fragments_coverage': pd.eval}).sort_values(by='Species')
    common_file.coverage = np.array(common_file.coverage) + np.array(file.coverage)
    common_file.fragments_coverage = per_fragments_sum_coverage(common_file.fragments_coverage.values.tolist(),
                                                           file.fragments_coverage.values.tolist())
    max_coverage = np.maximum(max_coverage, np.array(file.coverage))
    min_coverage = np.minimum(min_coverage, np.array(file.coverage))
common_file.insert(value=max_coverage, column='max_coverage',  loc=4)
common_file.insert(value=min_coverage, column='min_coverage',  loc=5)
common_file.max_coverage = (common_file.max_coverage.round(2) * 100).astype('int')
common_file.min_coverage = (common_file.min_coverage.round(2) * 100).astype('int')
common_file.coverage = ((common_file.coverage/len(files)).apply(lambda x: round(x, 2)) * 100).astype('int')
common_file.fragments_coverage = common_file.fragments_coverage.apply(lambda x: mean_per_fragments_coverage(x, len(files)))

common_file = common_file.sort_values(by='coverage', ascending=False)

common_file.reset_index(drop=True).to_csv(out, sep='\t')