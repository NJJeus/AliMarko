#!/usr/bin/env python3
"""
Script for analyzing sequence coverage and HMM profiles from FASTQ/FASTA files.
Processes coverage data and HMM search results to generate summary tables and heatmap visualizations.
"""

import argparse
import glob
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm


def parse_arguments():
    """Parse and validate command line arguments."""
    description = """
    The script processes coverage and HMM profile data from sequence files.
    Generates summary tables and heatmap visualizations for comparative analysis.
    """
    
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument('-i', '--cov_dir', type=str, 
                       help='Directory with coverage files')
    parser.add_argument('-c', '--output_coverage_table', type=str, 
                       help='Output coverage summary table file')
    parser.add_argument('-p', '--output_coverage_pic', type=str, 
                       help='Output coverage heatmap file')
    parser.add_argument('-m', '--hmm_dir', type=str, 
                       help='Directory with HMM profile files')
    parser.add_argument('-t', '--output_hmm_table', type=str, 
                       help='Output HMM summary table file') 
    parser.add_argument('-u', '--output_hmm_pic', type=str, 
                       help='Output HMM heatmap file')
    
    args = parser.parse_args()
    
    # Validate input directories
    if args.cov_dir:
        if not os.path.isdir(args.cov_dir):
            raise ValueError("Input folder with coverage files does not exist")
    else:
        raise ValueError('Missing required argument: -i/--cov_dir')
        
    if args.hmm_dir:
        if not os.path.isdir(args.hmm_dir):
            raise ValueError("Input folder with HMM report files does not exist")
    else:
        raise ValueError('Missing required argument: -m/--hmm_dir')
    
    return args


def filter_queries(group, column='Taxon'):
    """Filter a group to keep only rows matching the first value in specified column."""
    first_query = group[column].iloc[0]
    return group[group[column] == first_query]


def add_contamination_info(x, string_to_add):
    """Add contamination information string if value is not NA."""
    if pd.notna(x):
        return string_to_add + "Sequences of this virus have association with " + str(x) + "."
    return x


def process_coverage_data(cov_dir, output_table_path):
    """
    Process coverage data from CSV files in the input directory.
    Returns processed DataFrame and indexes of contaminated samples.
    """
    # Read and combine all coverage CSV files
    csv_files = glob.glob(os.path.join(cov_dir, '*.csv'))
    coverage_dfs = []
    
    for file in csv_files:
        df = pd.read_csv(file)
        filename = os.path.splitext(os.path.basename(file))[0]
        df['Sample'] = filename
        coverage_dfs.append(df.query('coverage > 0.05'))
    
    combined_coverage = pd.concat(coverage_dfs, ignore_index=True)
    
    # Process contamination information
    combined_coverage['Feature'] = combined_coverage['Feature'].apply(
        add_contamination_info, 
        string_to_add='CONTAMINATION_INFO:'
    )
    combined_coverage['Isolate_id'] = (
        combined_coverage['Isolate_id'] + combined_coverage['Feature'].fillna('')
    )
    combined_coverage['Isolate_id'] = (
        combined_coverage['Isolate_id']
        .str.replace('{', '[')
        .str.replace('}', ']')
    )
    
    # Create pivot table
    combined_coverage = combined_coverage.sort_values('coverage', ascending=False)
    coverage_pivot = combined_coverage.pivot(
        index='Isolate_id', 
        columns='Sample', 
        values='coverage'
    )
    
    # Sort by max coverage
    coverage_pivot['max_value'] = (
        coverage_pivot.max(axis=1) + coverage_pivot.mean(axis=1) * 0.5
    )
    coverage_pivot = coverage_pivot.sort_values(by='max_value', ascending=False)
    coverage_pivot.drop(columns='max_value', inplace=True)
    
    # Final processing and save
    coverage_pivot = coverage_pivot.fillna(0).head(25)
    coverage_pivot.index = pd.Series(coverage_pivot.index).apply(lambda x: x.split('[')[0])
    coverage_pivot.to_csv(output_table_path)
    
    # Get indexes of contaminated samples
    red_indexes = coverage_pivot.reset_index()[
        coverage_pivot.index.str.contains('CONTAMINATION_INFO:')
    ].index
    coverage_pivot.index = pd.Series(coverage_pivot.index).apply(
        lambda x: x.split('CONTAMINATION_INFO:')[0]
    )
    
    return coverage_pivot, red_indexes


def create_coverage_heatmap(coverage_data, red_indexes, output_path):
    """Create and save a heatmap visualization of coverage data."""
    plt.figure(figsize=(17, 11), dpi=300)
    heatmap = sns.heatmap(
        coverage_data, 
        cmap=sns.color_palette("mako_r", as_cmap=True), 
        fmt=".2f"
    )
    
    # Set labels and ticks
    heatmap.set_yticks(np.array(list(range(coverage_data.index.shape[0]))) + 0.5)
    heatmap.set_yticklabels(list(coverage_data.index), fontsize=14)
    plt.xticks(fontsize=14, rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.xlabel('\n Sample', fontsize=18)
    plt.ylabel('Taxon', fontsize=18)
    
    # Adjust layout and title
    plt.subplots_adjust(bottom=0.25, right=1.01, left=0.4, top=0.92)
    plt.title('Coverage Width Multisample Heatmap', fontsize=24)
    
    # Highlight contaminated samples
    y_labels = plt.gca().get_yticklabels()
    for index in red_indexes:
        y_labels[index].set_color('#C80000')
    
    # Add colorbar label
    heatmap.collections[0].colorbar.set_label("\n Coverage Width", fontsize=14)
    plt.savefig(output_path, format='jpg')
    plt.close()


def process_hmm_data(hmm_dir, output_table_path):
    """Process HMM data from CSV files in the input directory."""
    csv_files = glob.glob(os.path.join(hmm_dir, '*.csv'))
    hmm_dfs = []
    
    for file in csv_files:
        df = pd.read_csv(file, index_col=0).query('Score_ratio > 0.5')
        df = df.groupby('Name').apply(filter_queries).reset_index(drop=True)
        filename = os.path.splitext(os.path.basename(file))[0]
        df['Sample'] = filename
        hmm_dfs.append(df)
    
    combined_hmm = pd.concat(hmm_dfs, ignore_index=True)
    
    # Create pivot table
    hmm_pivot = combined_hmm[['Taxon', 'Sample', 'Score_ratio']].groupby(
        ['Taxon', 'Sample']
    ).sum().reset_index().pivot(
        index='Taxon', 
        columns='Sample', 
        values='Score_ratio'
    )
    
    # Sort by max score
    hmm_pivot['max_value'] = hmm_pivot.max(axis=1)
    hmm_pivot = hmm_pivot.sort_values(by='max_value', ascending=False)
    hmm_pivot.drop(columns='max_value', inplace=True)
    hmm_pivot = hmm_pivot.round(2)
    
    # Filter columns and save
    cols = [col for col in hmm_pivot.columns 
            if not (col.endswith('_') and 'MERGED' not in col)]
    hmm_pivot = hmm_pivot[cols].head(25)
    hmm_pivot.to_csv(output_table_path)
    
    return hmm_pivot


def create_hmm_heatmap(hmm_data, output_path):
    """Create and save a heatmap visualization of HMM data."""
    max_score = hmm_data.max().max() 
    min_score = hmm_data.min().min()
    
    if max_score == min_score:
        min_score, max_score = 0, min_score
    
    plt.figure(figsize=(17, 11), dpi=300)
    plt.subplots_adjust(left=0.2, bottom=0.22, right=1.01, top=0.95)
    
    heatmap = sns.heatmap(
        np.log(hmm_data).fillna(0)[:35], 
        cmap=sns.color_palette("mako_r", as_cmap=True)
    )
    
    # Set labels and title
    plt.title('HMM Multisample Heatmap', fontsize=18)
    plt.xlabel('\n Sample', fontsize=18)
    plt.ylabel('Taxon \n', fontsize=18)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=13, rotation=45, ha='right')
    
    # Add colorbar label
    heatmap.collections[0].colorbar.set_label("\n Sum of normalized scores", fontsize=14)
    plt.savefig(output_path, format='jpg')
    plt.close()


def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Process coverage data
    coverage_data, red_indexes = process_coverage_data(
        args.cov_dir, 
        args.output_coverage_table
    )
    create_coverage_heatmap(
        coverage_data, 
        red_indexes, 
        args.output_coverage_pic
    )
    
    # Process HMM data
    hmm_data = process_hmm_data(
        args.hmm_dir, 
        args.output_hmm_table
    )
    create_hmm_heatmap(
        hmm_data, 
        args.output_hmm_pic
    )


if __name__ == '__main__':
    main()
