import base64 


import pandas as pd
import os
import sys
import numpy as np
import glob
import argparse
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import urllib.parse
from Bio import SeqIO # New Import for reading FASTA files

# Add custom scripts to the Python path
sys.path.append('scripts')
# Assuming styles.py is available and contains the necessary classes
import styles 
from styles import Header


def check_condition(condition, message, mode='fail'):
    if condition:
        return True
    elif not condition and mode=='fail':
        print(message)
        sys.exit(1)
    elif not condition and mode=='ignore':
        return False


def create_palette(min_score, max_score):
    """Create a color palette for visualization."""
    if max_score == min_score:
        min_score, max_score = 0, min_score
    cmap = cm.ScalarMappable(cmap='Blues')
    cmap.set_array([min_score, max_score * 1.5])
    cmap.autoscale()
    return cmap


def get_color(palette, value):
    """Get a color from the palette based on the value."""
    try:
        color = mcolors.to_hex(palette.to_rgba(float(value)), keep_alpha=False)
    except Exception:
        return 'white'
    return color


def add_string(x, string_to_add):
    """Add a string to a value if it is not NaN."""
    if pd.notna(x):
        return string_to_add + "Sequences of this virus have association with " + str(x) + "."
    return x


def load_blast_results(blast_results_path):
    """Load and process BLAST results."""
    try:
        blast_results = pd.read_table(blast_results_path)
        blast_results.columns = ['qseqid', 'sseqid', 'stitle', 'salltitles', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    except Exception:
        blast_results = pd.DataFrame({i:[] for i in ['qseqid', 'sseqid', 'stitle', 'salltitles', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']})
    blast_results = blast_results.sort_values('bitscore', ascending=False).drop_duplicates('qseqid')
    blast_res_dict = blast_results[['qseqid', 'stitle', 'pident', 'length']].set_index('qseqid').to_dict('index')
    return {k: f"The best blastn hit: {v['stitle']}. \n Alignment length: {v['length']} bp, identity: {v['pident']}%" 
            for k, v in blast_res_dict.items()}


def validate_inputs(args):
    """Validate all input files and directories."""
    if not args.coverage is None:
        check_condition(os.path.isfile(args.coverage), "An input file does not exist")
        check_condition(os.path.isdir(args.drawings), "Drawings folder doesn't exist")
    check_condition(os.path.isdir(args.hmm_drawings), "HMM drawings folder doesn't exist")
    check_condition(os.path.isfile(args.hmm_report), "A HMM report file does not exist")
    if not check_condition(os.path.isdir(args.trees_folder), "Trees folder doesn't exist", mode='ignore'):
        args.trees_folder = False
    check_condition(os.path.isfile(args.blast_results), "An input file with BLAST results does not exist")
    check_condition(os.path.isfile(args.fasta_file), "The input FASTA file does not exist") # New validation


def load_fasta_sequences(fasta_file_path):
    """
    Load sequences from a FASTA file into a dictionary mapping ID to sequence string.
    """
    sequence_dict = {}
    try:
        for record in SeqIO.parse(fasta_file_path, "fasta"):
            sequence_dict[record.id] = str(record.seq)
        return sequence_dict
    except Exception as e:
        print(f"Error reading FASTA file {fasta_file_path}: {e}")
        sys.exit(1)


def process_ictv_coverage(coverage_file):
    """Process ICTV coverage data."""
    sample_name = os.path.splitext(os.path.basename(coverage_file))[0]
    ictv_coverage = pd.read_csv(coverage_file)
    try:
        # Add contamination info if present
        ictv_coverage['Feature'] = ictv_coverage['Feature'].apply(add_string, string_to_add='CONTAMINATION_INFO:')
        ictv_coverage['Species'] = ictv_coverage['Species'] + ictv_coverage['Feature'].fillna('')
    except Exception:
        None
    
    return sample_name, ictv_coverage


def create_introduction_table(ictv_coverage):
    """Create introduction table from ICTV coverage data."""
    introduction_frame = ictv_coverage[['Species', 'Host source', 'coverage', 'meandepth', 'Genus', 'Family', 'Realm']].query('coverage != 0')
    introduction_header = ['Species', 'Host source', 'Coverage width', 'Mean depth', 'Genus', 'Family', 'Realm']
    return introduction_frame.to_numpy(), introduction_header


def process_hmm_report(hmm_report_path):
    """Process HMM report data."""
    hmm_frame = pd.read_csv(hmm_report_path)[['Query', 'Taxon', 'Name', 'Positive terms', 'Score', 'Threshold', 'From', 'To', 'Score_ratio']]
    hmm_frame[['Score', 'Threshold', 'Score_ratio']] = hmm_frame[['Score', 'Threshold', 'Score_ratio']].astype('int')
    hmm_frame = hmm_frame.rename(columns={'Positive terms': 'Putative Protein', 'Query': 'HMM'})
    hmm_header = hmm_frame.drop(['Threshold', 'From', 'To', 'Score_ratio'], axis=1).columns
    hmm_table = hmm_frame.drop(['Threshold', 'From', 'To', 'Score_ratio'], axis=1).to_numpy()
    return hmm_frame, hmm_header, hmm_table


def process_ictv_drawings(ictv_coverage, drawings_folder):
    """Process ICTV drawings data."""
    ictv_drawings = ictv_coverage.query('coverage > 0.05').set_index('Host source')
    ictv_drawings['genbank_list'] = ictv_drawings['Virus GENBANK accession'].apply(
        lambda i: [el.split(":")[-1] for el in i.replace(' ', '').split(';')])
    fragment_cols = ['fragments_len', 'fragments_meandepth', 'fragments_coverage', 
                    'fragments_nucleotide_similarity', 'fragments_meanmapq', 'fragments_snps']
    ictv_drawings[fragment_cols] = ictv_drawings[fragment_cols].applymap(
        lambda x: x.strip('][').split(', '))
    ictv_drawings = ictv_drawings.groupby(level=0).agg({i: lambda x: list(x) for i in ictv_drawings.columns})
    
    return ictv_drawings


def create_host_dict(ictv_drawings, drawings_folder, palette):
    """Create host dictionary with coverage images."""
    host_dict = {}
    for host, row in ictv_drawings.iterrows():
        viruses = {}
        for virus_loc in range(len(row.Isolate_id)):
            virus_tech_name = row['Isolate_id'][virus_loc]
            virus_name = row['Species'][virus_loc]
            
            # contamination_info is empty if there is no CONTAMINATION_INFO:
            virus_name, _, contamination_info = virus_name.partition('CONTAMINATION_INFO:')
            virus_name = Header(virus_name, l=3, tooltip_text=contamination_info)

            # Prepare virus data for table
            virus_list = np.array([
                row.genbank_list[virus_loc],
                row.fragments_len[virus_loc],
                row.fragments_coverage[virus_loc],
                row.fragments_nucleotide_similarity[virus_loc],
                row.fragments_meandepth[virus_loc],
                row.fragments_meanmapq[virus_loc],
                row.fragments_snps[virus_loc]
            ]).T

            # Load images for this virus
            images = load_virus_images(drawings_folder, virus_tech_name, row.genbank_list[virus_loc])
            images_html = '\n'.join(images)

            # Create table and add to viruses dict
            table_html = styles.HTMLTable(
                virus_list,
                ['Fragment', 'Len', 'Coverage width', 'Nucleotide similarity', 'Mean Depth', 'MeanMAPQ', 'SNP Count'],
                palette,
                get_color
            ).render() + images_html
            
            viruses.update({virus_name: table_html})

        host_dict.update({Header(f'Host: {host}', l=2): styles.HTMLDetails(viruses).render()})

    return host_dict


def load_virus_images(drawings_folder, virus_tech_name, genbank_ids):
    """Load and encode images for a virus."""
    images = []
    for genbank_id in genbank_ids:
        try:
            picture_path = glob.glob(f"{drawings_folder}/{virus_tech_name}/*{genbank_id}*")[0]
            with open(picture_path, "rb") as image_file:
                encoded_image = base64.b64encode(image_file.read()).decode()
                images.append(styles.HTMLImage(encoded_image, type='jpg'))
        except Exception as e:
            print(e)
            continue
    return images


def order_host_dict(host_dict):
    """Order hosts in the host dictionary according to predefined order."""
    predefined_order = [
        'vertebrates', 'marine (S)', 'invertebrates', 'invertebrates, vertebrates',
        'invertebrates (S)', 'invertebrates, plants', 'plants, invertebrates',
        'fungi', 'plants', 'protists (S)', 'algae', 'protists', 'bacteria',
        'archaea', 'sewage (S)', 'soil (S)', 'freshwater (S)'
    ]
    order = [Header(f'Host: {host}', l=2) for host in predefined_order]
    missing_hosts = [key for key in host_dict.keys() if key not in order]
    order += missing_hosts
    return {key: host_dict[key] for key in sorted(host_dict.keys(), key=lambda x: order.index(x))}




def process_hmm_results(hmm_frame, hmm_drawings_folder, trees_folder, blast_res_dict, palette, sequences_dict):
    """Process HMM results and generate HTML content."""
    positive_contigs = hmm_frame.query('Score_ratio > 0.5').groupby('Name').count().query('Score > 0').reset_index().Name
    list_of_contigs = hmm_frame['Name'][hmm_frame['Name'].isin(positive_contigs)].unique()
    images_hmm = []

    for contig in list_of_contigs:
        # 2. Select the sequence record where ID matches 'contig'
        contig_sequence = sequences_dict.get(contig, "") 

        # 3. Add sequence utility HTML
        sequence_utility_html = styles.create_sequence_utility_html(contig, contig_sequence)

        table_contig = hmm_frame.query(f'Name == "{contig}"').drop(['Threshold', 'From', 'To', 'Score_ratio'], axis=1)
        models = hmm_frame.query(f'Name == "{contig}"').Taxon.unique().tolist()
        models = ",".join(models) if len(models) < 3 else f"({len(models)} taxa)"

        blastn_info = Header(blast_res_dict.get(contig, 'No blast hits for this contig'), l=4)

        try:
            picture_path = glob.glob(f"{hmm_drawings_folder}/*{contig}*")[0]
            with open(picture_path, "rb") as image_file:
                encoded_image = base64.b64encode(image_file.read()).decode()
                html_image = styles.HTMLImage(encoded_image, 'svg')
                
                trees_pictures = process_trees_for_contig(trees_folder, contig, table_contig)
                
                details_image = styles.HTMLDetails({
                    Header(f'{contig}:{models}', l=3): (
                        # Inject sequence utility here
                        sequence_utility_html + 
                        styles.HTMLTable(table_contig.to_numpy(), table_contig.columns, palette, get_color).render() + 
                        '\n' + blastn_info + '\n' + html_image + '\n' + trees_pictures
                    )
                }).render()
                images_hmm.append(details_image)
        except Exception:
            # print(f"Warning: Could not process drawings/images for contig {contig}")
            continue

    return '\n'.join(images_hmm)


def process_trees_for_contig(trees_folder, contig, table_contig):
    """Process phylogenetic trees for a contig."""
    trees_pictures = ''
    if not trees_folder:
        return '\n'
    for tree in [i for i in glob.glob(f"{trees_folder}/*") if contig in i]:
        with open(tree, "rb") as tree_image_file:
            tree_encoded_image = base64.b64encode(tree_image_file.read()).decode()
            tree_HMM = tree.split('__')[1]
            try:
                protein_HMM = table_contig.query(f'HMM == "{tree_HMM}"')['Putative Protein'].unique()[0]
            except Exception:
                continue
            tree_name = Header(f"Phylogenetic Tree of {tree_HMM}-Matched Amino Acid Sequences Putative protein: {protein_HMM}", l=3)
            html_tree = styles.HTMLImage(tree_encoded_image, 'jpg')
            trees_pictures = trees_pictures + '\n' + tree_name + html_tree
    return trees_pictures



def generate_html_output(sample_name, style, introduction_table, introduction_header, 
                        hmm_table, hmm_header, host_dict, images_hmm):
    """Generate final HTML output, appending the necessary JavaScript."""
    greetings = Header(f"Sample: {sample_name}", l=1)
    table_html = styles.HTMLTable(introduction_table, introduction_header, None, get_color).render()
    hmm_greetings = Header("HMM Hit Summary", l=2)
    table_hmm = styles.HTMLTable(hmm_table, hmm_header, None, get_color).render()
    details = styles.HTMLDetails(host_dict).render()
    
    # 1. Generate the main report body
    report_body = (
        greetings + 
        Header("Mapping Summary", l=2) + table_html + 
        hmm_greetings + table_hmm + 
        Header("Mapping Details", l=2) + details + 
        Header("HMM Module Results by Contig", l=2) + 
        styles.HTMLDetails({Header('Contigs', l=2): images_hmm}).render() +
        Header('', l=2)
    )
    
    # 2. Append the utility script and custom message box to the body end
    utility_script = styles.get_utility_script()
    
    # 3. Insert the report body and script into the HTML template (assuming styles.get_style() provides a full HTML structure)
    # This assumes styles.get_style() contains the opening <html>, <head>, and <body> tags, and the closing tags are appended later.
    final_output = style.replace('Sample Name', sample_name) + report_body + utility_script
    
    return final_output


def main():
    # Argument parsing
    description = "Script that concatenates ICTV coverage data with drawings of coverage."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-c', '--coverage', type=str, required=False, default=None, help='A file with ICTV coverage')
    parser.add_argument('-o', '--output', type=str, required=False, default=None, help='An output file')
    parser.add_argument('-d', '--drawings', type=str, required=False, default=None, help='A folder with drawings')
    parser.add_argument('-a', '--hmm_drawings', type=str, required=True,  help='A folder with HMM drawings')
    parser.add_argument('-m', '--hmm_report', type=str, required=True, help='A HMM report file')
    parser.add_argument('-t', '--trees_folder', type=str, required=True, help='A folder with trees')
    parser.add_argument('-b', '--blast_results', type=str, required=True, help='BLAST tsv table')
    # 1. Add FASTA file argument
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='Path to the input FASTA file with all sequences.')
    args = parser.parse_args()

    # Validate inputs
    validate_inputs(args)

    # Load BLAST results
    try:
        blast_res_dict = load_blast_results(args.blast_results)
    except Exception:
        blast_res_dict = {}

    # Load FASTA sequences (New Step)
    sequences_dict = load_fasta_sequences(args.fasta_file)

    sample_name = os.path.splitext(os.path.basename(args.hmm_report))[0]

    if not args.coverage is None:
        # Process ICTV coverage data
        sample_name, ictv_coverage = process_ictv_coverage(args.coverage)

        introduction_table, introduction_header = create_introduction_table(ictv_coverage)

        # Process ICTV drawings
        ictv_drawings = process_ictv_drawings(ictv_coverage, args.drawings)
        host_dict = create_host_dict(ictv_drawings, args.drawings, None)
        host_dict = order_host_dict(host_dict)
    else:
        ictv_coverage, introduction_header, introduction_table, host_dict = '', '', '', dict()

    # Process HMM report
    hmm_frame, hmm_header, hmm_table = process_hmm_report(args.hmm_report)


    # Process HMM results (Pass sequences_dict)
    images_hmm = process_hmm_results(hmm_frame, args.hmm_drawings, args.trees_folder, blast_res_dict, None, sequences_dict)

    # Generate and write final output
    out = generate_html_output(
        sample_name,
        styles.get_style(),
        introduction_table,
        introduction_header,
        hmm_table,
        hmm_header,
        host_dict,
        images_hmm
    )

    with open(args.output, 'w') as f:
        f.write(out)


if __name__ == "__main__":
    main()

