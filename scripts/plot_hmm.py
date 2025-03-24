import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import os
import argparse

# Set Seaborn theme
sns.set_theme()

def parse_arguments():
    """Parse command-line arguments."""
    description = """The script gets a report CSV file and returns a folder with plots for all contigs."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Input CSV file')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Output directory for plots')
    return parser.parse_args()

def validate_inputs(input_file, output_dir):
    """Validate input file and output directory."""
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file does not exist: {input_file}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def correct_text_position(ax, x, y, text, fontsize, max_iterations=130):
    """
    Adjust text position to avoid overlap with plot boundaries.
    
    Args:
        ax: Matplotlib axis object.
        x, y: Initial text position.
        text: Text to display.
        fontsize: Initial font size.
        max_iterations: Maximum number of iterations to adjust position.
    
    Returns:
        str: Adjusted text.
    """
    for _ in range(max_iterations):
        text_obj = ax.text(x, y, text, fontdict={"fontsize": fontsize}, ha='center')
        bbox = ax.get_window_extent()
        xmin, xmax = bbox.x0 + 250, bbox.x1 - 250
        tbox = text_obj.get_window_extent()
        x0, x1 = tbox.x0, tbox.x1
        xlim = ax.get_xlim()
        step = (xlim[1] - xlim[0]) / 100

        if xmin < x0 and x1 <= xmax:
            return text
        if xmin > x0 and x1 <= xmax:
            x += step
        elif x0 >= xmin and x1 > xmax:
            x -= step
        elif x0 < xmin and x1 > xmax:
            fontsize -= 1
        text_obj.remove()
    return text

def setup_plot(ax, contig_name, length_contig):
    """
    Set up the plot with title and limits.
    
    Args:
        ax: Matplotlib axis object.
        contig_name: Name of the contig.
        length_contig: Length of the contig.
    """
    ax.text(length_contig / 2, 0.7, f"HMM hits to {contig_name}", ha='center')
    ax.set_xlim([0, length_contig])
    ax.set_yticks([])

def add_colorbar(ax, cmap, min_score, max_score):
    """
    Add a colorbar to the plot.
    
    Args:
        ax: Matplotlib axis object.
        cmap: Colormap object.
        min_score: Minimum score for colorbar.
        max_score: Maximum score for colorbar.
    """
    cbar = plt.colorbar(cmap, orientation='vertical', ticks=np.linspace(min_score, max_score, num=5), ax=ax)
    cbar.set_label('Score')

def plot_arrows(ax, contig_data, cmap, length_contig):
    """
    Plot arrows for HMM hits.
    
    Args:
        ax: Matplotlib axis object.
        contig_data: DataFrame containing contig data.
        cmap: Colormap object.
        length_contig: Length of the contig.
    """
    for i, (index, row) in enumerate(contig_data.iterrows(), start=-4):
        color = cmap.to_rgba(row.Score)
        ax.arrow(row.From, i, row.To - row.From, 0,
                 color=color, width=0.1, head_width=0.2,
                 length_includes_head=True,
                 head_length=max(abs(row.To - row.From) / 10, length_contig / 100))
        arrow_text = f'{row.Query}:{row["Positive terms"]}, {row.Taxon}'
        correct_text_position(ax, (row.From + row.To) / 2, i + 0.3, arrow_text, 9)

def plot_contig(ax, contig_data, contig_name):
    """
    Plot HMM hits for a single contig.
    
    Args:
        ax: Matplotlib axis object.
        contig_data: DataFrame containing contig data.
        contig_name: Name of the contig.
    """
    length_contig = contig_data.Length_contig.max()
    setup_plot(ax, contig_name, length_contig)

    max_score = contig_data.Score.max()
    min_score = contig_data.Score.min()
    if max_score == min_score:
        min_score, max_score = 0, min_score

    cmap = cm.ScalarMappable(cmap='jet')
    cmap.set_array([min_score, max_score])
    cmap.autoscale()

    plot_arrows(ax, contig_data, cmap, length_contig)
    add_colorbar(ax, cmap, min_score, max_score)

    ax.plot([0, length_contig], [0, 0], linewidth=4.0)
    ax.set_ylim([-contig_data.shape[0] - 2, 0.2])

def main():
    """Main function to generate plots for contigs."""
    args = parse_arguments()
    validate_inputs(args.input_file, args.output_dir)

    data = pd.read_csv(args.input_file, index_col=0)

    for contig in data.Name.unique():
        plt.clf()
        fig = plt.figure(dpi=400, figsize=(8, 4.5))
        ax = fig.add_subplot(1, 1, 1)

        contig_data = data.query(f'Name == "{contig}"').head(4)
        plot_contig(ax, contig_data, contig)

        plt.savefig(f'{args.output_dir}/{contig}.svg')
        plt.close()

if __name__ == "__main__":
    main()
