import argparse
import os
import pandas as pd
import numpy as np

def check_condition(condition, message):
    """Check a condition and exit with a message if it fails."""
    if not condition:
        print(message)
        exit()

def parse_arguments():
    """Parse command-line arguments."""
    description = "Script that concatenates a samtools coverage output and a samtools view quality with an ICTV metadata."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-c', '--coverage', type=str, required=True, help='A file with samtools coverage output')
    parser.add_argument('-o', '--output', type=str, required=True, help='An output file')
    parser.add_argument('-m', '--tmp_output', type=str, required=True, help='An output file with data for the plot_coverage script')
    parser.add_argument('-t', '--tables_folder', type=str, required=True, help='A folder with tables ictv_taxo and genbank accessions')
    parser.add_argument('-s', '--snps_file', type=str, required=True, help='A file with rname, number of snp, number of deep sites as columns and virus samples as rows')

    args = parser.parse_args()

    # Validate file and directory paths
    check_condition(os.path.isfile(args.coverage), "Coverage file doesn't exist")
    check_condition(os.path.isdir(args.tables_folder), "Tables folder doesn't exist")
    check_condition(os.path.isfile(args.snps_file), "SNPs file doesn't exist")

    return args

def load_and_preprocess_data(coverage_file, snps_file, tables_folder):
    """Load and preprocess data from input files."""
    # Read samtools coverage output and drop excess columns
    coverage_data = pd.read_table(coverage_file).rename(columns={'#rname': 'rname', 'meanmapq': 'meanmapq'})
    coverage_data = coverage_data.drop(['numreads', 'covbases', 'meanbaseq'], axis=1)

    # Read SNPs file
    snps_data = pd.read_csv(snps_file)

    # Load ICTV metadata and preprocess
    ictv_data = pd.read_excel(f'{tables_folder}/ictv_taxo.xlsx').dropna(subset='Virus GENBANK accession')
    ictv_data['genbank_list'] = ictv_data['Virus GENBANK accession'].apply(
        lambda x: [el.split(":")[-1] for el in x.replace(' ', '').split(';')]
    )
    ictv_data['Isolate_id'] = ictv_data['Species'] + '{' + ictv_data['Isolate Sort'].astype(str) + '}'

    # Create a mapping between accession numbers and Isolate_id
    accessions = []
    indeces = []
    for _, row in ictv_data.iterrows():
        for accession in row['genbank_list']:
            accessions.append(accession)
            indeces.append(row['Isolate_id'])
    taxo_index = pd.DataFrame({'accession': accessions, 'ictv_taxo_index': indeces})

    return coverage_data, snps_data, ictv_data, taxo_index

def merge_and_process_data(coverage_data, snps_data, taxo_index, ictv_data):
    """Merge and process data to generate final output."""
    # Merge coverage data with SNPs data
    merged_data = pd.merge(coverage_data, snps_data, on='rname', how='left')
    merged_data['nucleotide_similarity'] = 1 - merged_data['snps'] / (merged_data['deep_sites'] + 1)

    # Extract accession number from rname
    merged_data['tmp_rname'] = merged_data['rname']
    merged_data['rname'] = merged_data['rname'].apply(lambda x: x.split('.')[0])
    merged_data = merged_data.set_index('rname')

    # Merge with taxo_index and ICTV metadata
    merged_data = pd.merge(merged_data.reset_index(), taxo_index, left_on='rname', right_on='accession', how='left')
    merged_data = pd.merge(merged_data, ictv_data, left_on='ictv_taxo_index', right_on='Isolate_id', how='left')

    # Generalize fragments of one virus
    merged_data[['coverage', 'meanmapq']] = merged_data[['coverage', 'meanmapq']] * 0.01
    merged_data['len'] = merged_data['endpos'] - merged_data['startpos'] + 1
    merged_data[['weighted_coverage', 'weighted_quality', 'weighted_meandepth']] = (
        merged_data[['coverage', 'meanmapq', 'meandepth']].T * merged_data['len']).T

    # Create per-fragment columns
    per_fragment_columns = ['len', 'coverage', 'meanmapq', 'meandepth', 'deep_sites', 'nucleotide_similarity', 'snps']
    for col in per_fragment_columns:
        merged_data[f'fragments_{col}'] = merged_data[col]

    # Aggregate data
    agg_dict = {
        **{col: 'sum' for col in ['len', 'snps', 'deep_sites', 'weighted_coverage', 'weighted_quality', 'weighted_meandepth']},
        **{col: lambda x: [round(c, 4) for c in x] for col in [f'fragments_{c}' for c in per_fragment_columns]},
        'tmp_rname': lambda x: list(x),
        **{col: 'first' for col in ['rname', 'Realm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus isolate designation', 'Virus name(s)', 'Virus GENBANK accession', 'Host source', 'genbank_list', 'Feature']}
    }
    aggregated_data = merged_data.groupby('Isolate_id').agg(agg_dict)

    # Calculate final metrics
    aggregated_data['coverage'] = (aggregated_data['weighted_coverage'] / aggregated_data['len']).round(2)
    aggregated_data['meanmapq'] = (aggregated_data['weighted_quality'] / aggregated_data['len'] * 100).round(2)
    aggregated_data['meandepth'] = (aggregated_data['weighted_meandepth'] / aggregated_data['len']).round(2)
    aggregated_data['nucleotide_similarity'] = 1 - aggregated_data['snps'].div(aggregated_data['deep_sites'].replace(0, np.inf))
    aggregated_data['deep_sites'] = aggregated_data['deep_sites'].astype(int)
    aggregated_data = aggregated_data.sort_values(by='coverage', ascending=False)

    return aggregated_data

def save_output(data, output_file, tmp_output):
    """Save the final output and temporary data for plotting."""
    # Save temporary data for plotting
    plot_data = data.query('coverage > 0.05 & meanmapq > 19')
    tmp_data = pd.DataFrame({
        'tmp_names': [name for sublist in plot_data['tmp_rname'] for name in sublist],
        'fragments_len': [length for sublist in plot_data['fragments_len'] for length in sublist]
    }, index=[index for index, _ in plot_data.iterrows() for _ in range(len(plot_data.loc[index, 'tmp_rname']))])
    tmp_data.to_csv(tmp_output, header=False)

    # Save final output
    output_columns = [
        'Virus name(s)', 'Host source', 'len', 'deep_sites', 'coverage', 'meandepth', 'meanmapq', 'snps', 'nucleotide_similarity',
        *[f'fragments_{col}' for col in ['len', 'coverage', 'meanmapq', 'meandepth', 'deep_sites', 'nucleotide_similarity', 'snps']],
        'Virus GENBANK accession', 'Realm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Feature'
    ]
    data[output_columns].to_csv(output_file)

def main():
    args = parse_arguments()
    coverage_data, snps_data, ictv_data, taxo_index = load_and_preprocess_data(args.coverage, args.snps_file, args.tables_folder)
    processed_data = merge_and_process_data(coverage_data, snps_data, taxo_index, ictv_data)
    save_output(processed_data, args.output, args.tmp_output)

if __name__ == "__main__":
    main()
