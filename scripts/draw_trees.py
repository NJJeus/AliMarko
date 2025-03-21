from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import glob
import argparse
import os
import re


def check_condition(condition, message):
    """Check a condition and exit with a message if the condition is not met."""
    if not condition:
        print(message)
        exit()


def format_label(label):
    """Format the label by removing unwanted patterns and truncating it."""
    label = re.sub(r'{[^{}]*}', '', str(label)).split('#')[0].split("|")[0]
    return str(label)[:28]


def round_bootstrap_values(clade):
    """Round bootstrap values to one decimal place recursively."""
    if clade.confidence:
        clade.confidence = int(clade.confidence * 10) / 10
    for child in clade.clades:
        round_bootstrap_values(child)


def get_bootstrap_color(bootstrap):
    """Generate a color based on the bootstrap value."""
    return f"hsl({bootstrap}, 100%, 50%)"


def parse_arguments():
    """Parse command-line arguments."""
    description = "A script to generate phylogenetic trees from tree files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input_folder', type=str, required=True, help='Input folder containing tree files')
    parser.add_argument('-o', '--output_folder', type=str, required=True, help='Output folder to save generated trees')
    return parser.parse_args()


def create_output_folder(output_folder):
    """Create the output folder if it does not exist."""
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


def extract_family_info(leaf_nodes):
    """Extract family information from leaf nodes."""
    families = {}
    for leaf in leaf_nodes:
        family = leaf.name.split('|')[-1].split('#')[0]
        if "CONTIG" in leaf.name:
            family = "_".join(leaf.name.split('_')[:2])
        families[leaf.name] = family
    return families


def generate_color_palette(unique_families):
    """Generate a color palette for the families."""
    palette = sns.color_palette("husl", n_colors=len(unique_families))
    legend_colors = {family: color for family, color in zip(unique_families, palette)}

    # Assign specific color for CONTIG families
    for family in unique_families:
        if family.startswith('CONTIG'):
            legend_colors[family] = (0.890, 0.102, 0.110)
    return legend_colors


def create_legend_patches(legend_colors):
    """Create legend patches for the plot."""
    return [mpatches.Patch(color=color, label=label) for label, color in legend_colors.items()]


def plot_tree(tree, label_colors, legend_patches, name, output_folder):
    """Plot the phylogenetic tree and save it to a file."""
    fig, ax = plt.subplots(figsize=(11.5, 5.94), dpi=250)

    # Draw the tree
    Phylo.draw(tree, axes=ax, label_func=format_label, label_colors=label_colors)

    # Add the legend
    ax.legend(handles=legend_patches, bbox_to_anchor=(1, 1), loc='upper left')

    # Set plot title and labels
    ax.set_title(f"Phylogeny of {name.split('__')[1]}-Matched AA Sequences from {name.split('__')[0].replace('test_', '')} contig")
    ax.set_xlabel('Branch Length')
    ax.set_ylabel('')
    ax.set_yticklabels([])
    ax.tick_params(axis='y', which='both', length=0)

    # Adjust layout to make space for the legend
    fig.subplots_adjust(right=0.8, left=0.05)
    x_min, x_max = ax.get_xlim()
    ax.set_xlim([x_min, x_max+(x_max-x_min)*0.5])

    # Save the figure
    output_path = os.path.join(output_folder, f"{name}_tree.jpg")
    fig.savefig(output_path, format='jpg')
    plt.close(fig)

def process_tree_file(file, output_folder):
    """Process a single tree file and generate a phylogenetic tree plot."""
    name = os.path.basename(file).replace('_tree.treefile', '')
    print(f"Processing: {name}")

    try:
        tree = Phylo.read(file, "newick")
    except Exception as e:
        print(f"No tree found in {file}: {e}")
        return

    try:
        tree.root_at_midpoint()
    except Exception as e:
        print(f"Failed to root tree at midpoint: {e}")
        return

    round_bootstrap_values(tree.root)
    leaf_nodes = tree.get_terminals()

    # Extract family information
    families = extract_family_info(leaf_nodes)
    unique_families = set(families.values())

    # Generate color palette and legend
    legend_colors = generate_color_palette(unique_families)
    label_colors = {format_label(leaf): legend_colors[family] for leaf, family in families.items()}
    legend_patches = create_legend_patches(legend_colors)

    # Plot and save the tree
    plot_tree(tree, label_colors, legend_patches, name, output_folder)


def main():
    """Main function to execute the script."""
    args = parse_arguments()

    check_condition(os.path.exists(args.input_folder), "Input folder does not exist")
    create_output_folder(args.output_folder)

    tree_files = glob.glob(os.path.join(args.input_folder, '*.treefile'))
    for file in tree_files:
        process_tree_file(file, args.output_folder)


if __name__ == "__main__":
    main()
