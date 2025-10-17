import base64
import pandas as pd
import os
import sys
import numpy as np
import argparse
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Add custom scripts to the Python path
sys.path.append('scripts')
from styles import HTMLTable, HTMLImage, Header, get_style


def check_condition(condition, message):
    """Check a condition and exit with a message if it fails."""
    if not condition:
        print(message)
        sys.exit(1)


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
        return None
    return color


def validate_inputs(args):
    """Validate all input files."""

    if not args.coverage is None:
        check_condition(os.path.isfile(args.coverage), "Input file does not exist: coverage")
        check_condition(os.path.isfile(args.coverage_heatmap), "Input file does not exist: coverage_heatmap")
    check_condition(os.path.isfile(args.hmm_report), "Input file does not exist: hmm_report")
    check_condition(os.path.isfile(args.hmm_heatmap), "Input file does not exist: hmm_heatmap")


def process_coverage_data(coverage_file):
    """Process coverage data."""
    coverage_df = pd.read_csv(coverage_file)
    sample_name = os.path.splitext(os.path.basename(coverage_file))[0]
    coverage_table = coverage_df.fillna(0).to_numpy()
    max_cov, min_cov = coverage_table[:, 1:].astype('float').max().max(), coverage_table[:, 1:].astype('float').min().min()
    
    return {
        'header': list(coverage_df.columns),
        'table': coverage_table,
        'min': min_cov,
        'max': max_cov
    }, sample_name


def process_hmm_data(hmm_file):
    """Process HMM data."""
    hmm_df = pd.read_csv(hmm_file)
    hmm_table = hmm_df.fillna(0).to_numpy()
    max_hmm, min_hmm = hmm_table[:, 1:].astype('float').max().max(), hmm_table[:, 1:].astype('float').min().min()
    
    return {
        'header': list(hmm_df.columns),
        'table': hmm_table,
        'min': min_hmm,
        'max': max_hmm
    }


def load_image_as_html(image_path):
    """Load an image and return as HTML."""
    with open(image_path, "rb") as f:
        return HTMLImage(base64.b64encode(f.read()).decode(), 'jpg')


def generate_report_tables(coverage_data, hmm_data):
    """Generate HTML tables for coverage and HMM data."""
    
    if coverage_data:
        cov_palette = create_palette(coverage_data['min'], coverage_data['max'])
        
        coverage_table = HTMLTable(
            coverage_data['table'],
            coverage_data['header'],
            cov_palette,
            get_color
        ).render()
    else:
        coverage_table = ''

    hmm_palette = create_palette(hmm_data['min'], hmm_data['max'])
    hmm_table = HTMLTable(
        hmm_data['table'],
        hmm_data['header'],
        hmm_palette,
        get_color
    ).render()
    
    return coverage_table, hmm_table


def compose_html_report(sample_name, coverage_table, hmm_table, coverage_image, hmm_image):
    """Compose the complete HTML report."""
    html_parts = [
        get_style().replace('Sample Name', sample_name),
        Header('AliMarko multisample report', l=1),
        Header('Mapping Multisample Results', l=2),
        Header('Mapping Coverage Width Heatmap', l=3),
        coverage_image,
        Header('Mapping Coverage Width', l=3),
        coverage_table,
        Header('HMM Multisample Results', l=2),
        Header('HMM Score Heatmap', l=3),
        hmm_image,
        Header('HMM Score', l=3),
        hmm_table,
        Header('', l=3)
    ]
    
    return ''.join(html_parts)


def main():
    """Main function to execute the script."""
    description = "Script that concatenates an ICTV coverage data with drawings of coverage"
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-c', '--coverage', type=str, default=None,
                       help='A file with ICTV coverage')
    parser.add_argument('-e', '--coverage_heatmap', type=str, default=None,
                       help='A file with ICTV heatmap coverage')
    parser.add_argument('-m', '--hmm_report', type=str, required=True,
                       help='A HMM report file')
    parser.add_argument('-p', '--hmm_heatmap', type=str, required=True,
                       help='A file with HMM heatmap coverage')
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='An output file')
    
    args = parser.parse_args()
    validate_inputs(args)
    
    if not args.coverage is None:
        coverage_data, sample_name = process_coverage_data(args.coverage)
        coverage_image = load_image_as_html(args.coverage_heatmap)
    else:
        coverage_data = False
        sample_name = ''
        coverage_image = ''
    hmm_data = process_hmm_data(args.hmm_report)
    
    coverage_table, hmm_table = generate_report_tables(coverage_data, hmm_data)
    
    hmm_image = load_image_as_html(args.hmm_heatmap)
    
    html_report = compose_html_report(
        sample_name,
        coverage_table,
        hmm_table,
        coverage_image,
        hmm_image
    )
    
    with open(args.output, 'w') as f:
        f.write(html_report)


if __name__ == "__main__":
    main()
