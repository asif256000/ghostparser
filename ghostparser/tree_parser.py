"""Tree parsing and standardization module using BioPython."""

import argparse
from io import StringIO
from pathlib import Path

from Bio import Phylo


def read_tree_file(filepath):
    """Read Newick format trees from a file using BioPython.

    Args:
        filepath: Path to the Newick format tree file.

    Returns:
        A list of Bio.Phylo tree objects.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file contains invalid Newick format.
    """
    try:
        trees = list(Phylo.parse(filepath, "newick"))
        if not trees:
            raise ValueError(f"No valid trees found in {filepath}")
        return trees
    except FileNotFoundError:
        raise FileNotFoundError(f"Tree file not found: {filepath}")
    except Exception as e:
        raise ValueError(f"Error parsing tree file {filepath}: {e}")


def calculate_average_support(tree):
    """Calculate the average support value for a tree.

    Args:
        tree: A Bio.Phylo tree object.

    Returns:
        The average support value across all internal nodes, or None if no support values found.
    """
    support_values = []
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.confidence is not None:
            support_values.append(clade.confidence)

    if not support_values:
        return None

    return sum(support_values) / len(support_values)


def remove_support_values(tree):
    """Remove support values (bootstrap/posterior probabilities) from internal nodes.

    This function removes confidence values that appear after internal nodes
    in Newick format trees like: (A,B)0.95:0.5 -> (A,B):0.5

    Args:
        tree: A Bio.Phylo tree object.

    Returns:
        The tree with support values removed from all internal nodes.
    """
    # Traverse all internal nodes (non-terminals) and remove confidence values
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.confidence = None
    return tree


def standardize_tree(tree):
    """Standardize a BioPython tree object.

    This function takes a BioPython tree object and standardizes it by:
    - Removing support values (bootstrap/posterior probabilities) from internal nodes
    - Ensuring consistent Newick format output
    - Preserving branch lengths and taxa names

    Args:
        tree: A Bio.Phylo tree object.

    Returns:
        A Bio.Phylo tree object (standardized).
    """
    tree = remove_support_values(tree)
    return tree


def format_newick_with_precision(tree, decimal_places=10):
    """Format a tree to Newick string with custom decimal precision for branch lengths.

    Trailing zeros after decimal places are removed to keep the output clean while
    preserving necessary precision.

    Args:
        tree: A Bio.Phylo tree object.
        decimal_places: Number of decimal places for branch lengths. Default: 10.

    Returns:
        A Newick format string with the specified decimal precision and trailing zeros removed.
    """

    def format_branch_length(branch_length):
        """Format branch length with precision and remove trailing zeros."""
        # Format with specified decimal places
        formatted = f"{branch_length:.{decimal_places}f}"
        # Remove trailing zeros after the decimal point
        if "." in formatted:
            formatted = formatted.rstrip("0").rstrip(".")
        return formatted

    def format_clade(clade):
        """Recursively format a clade to Newick string."""
        if clade.is_terminal():
            # Leaf node: just the name and branch length
            result = clade.name or ""
        else:
            # Internal node: format all children recursively
            children = [format_clade(c) for c in clade.clades]
            result = "(" + ",".join(children) + ")"

        # Add branch length if present
        if clade.branch_length is not None:
            # Format branch length with custom precision and remove trailing zeros
            result += f":{format_branch_length(clade.branch_length)}"

        return result

    # Format the entire tree
    newick_str = format_clade(tree.root) + ";"
    return newick_str


def write_clean_trees(trees, output_filepath, decimal_places=15):
    """Write standardized trees to a file using BioPython with custom precision.

    Args:
        trees: A list of Bio.Phylo tree objects.
        output_filepath: Path where the cleaned trees will be written.
        decimal_places: Number of decimal places for branch lengths. Default: 10.
    """
    with open(output_filepath, "w") as f:
        for tree in trees:
            # Write Newick format with custom decimal precision
            newick_str = format_newick_with_precision(tree, decimal_places)
            f.write(newick_str + "\n")


def clean_and_save_trees(input_filepath, output_filepath, min_avg_support=0.5):
    """Read, standardize, and save trees from a file using BioPython.

    Trees with average support values below the threshold are filtered out.

    Args:
        input_filepath: Path to the input Newick format tree file.
        output_filepath: Path where the cleaned trees will be written.
        min_avg_support: Minimum average support value threshold. Default: 0.5.

    Returns:
        A tuple of (cleaned_trees_list, dropped_trees_info_dict).
        dropped_trees_info_dict maps tree index to average support value.
    """
    # Read the trees using BioPython
    trees = read_tree_file(input_filepath)

    # Track dropped trees
    dropped_trees = {}
    cleaned_trees = []

    # Process each tree
    for idx, tree in enumerate(trees, start=1):
        # Calculate average support before standardization
        avg_support = calculate_average_support(tree)

        # Check if tree meets minimum support threshold
        if avg_support is not None and avg_support < min_avg_support:
            dropped_trees[idx] = avg_support
            continue

        # Standardize the tree (removes support values)
        standardized = standardize_tree(tree)
        cleaned_trees.append(standardized)

    # Write to output file using BioPython
    write_clean_trees(cleaned_trees, output_filepath)

    return cleaned_trees, dropped_trees


def get_clean_filename(filepath):
    """Generate a clean filename with _clean prefix.

    Args:
        filepath: The original file path.

    Returns:
        The path with _clean inserted before the file extension.
    """
    path = Path(filepath)
    return str(path.parent / f"{path.stem}_clean{path.suffix}")


def main():
    """Main entry point for the tree parser module."""
    parser = argparse.ArgumentParser(
        description="Ghost parser for identifying ghost introgressions in phylogenetic trees."
    )

    parser.add_argument("-st", "--species_tree", required=True, help="Path to the species tree file in Newick format")
    parser.add_argument("-gt", "--gene_trees", required=True, help="Path to the gene trees file in Newick format")
    parser.add_argument("-og", "--outgroup", required=True, help="Outgroup species identifier")

    args = parser.parse_args()

    # Convert to Path objects and resolve relative to current directory
    species_tree_path = Path(args.species_tree)
    gene_trees_path = Path(args.gene_trees)

    # Validate input files exist
    if not species_tree_path.exists():
        print(f"Error: Species tree file not found: {args.species_tree}")
        return

    if not gene_trees_path.exists():
        print(f"Error: Gene trees file not found: {args.gene_trees}")
        return

    # Get output filenames (keep them relative)
    species_tree_clean = get_clean_filename(str(species_tree_path))
    gene_trees_clean = get_clean_filename(str(gene_trees_path))

    print(f"Processing species tree: {args.species_tree}")
    print(f"Processing gene trees: {args.gene_trees}")
    print(f"Outgroup: {args.outgroup}")
    print("Support threshold: 0.5")
    print()

    # Clean and save species tree
    try:
        species_trees, dropped_species = clean_and_save_trees(str(species_tree_path), species_tree_clean)
        print(f"✓ Species tree cleaned and saved to: {species_tree_clean}")
        print(f"  Processed {len(species_trees)} tree(s)")
        if dropped_species:
            print(f"  ⚠ Dropped {len(dropped_species)} tree(s) with avg support < 0.5:")
            for idx, avg_support in dropped_species.items():
                print(f"    - Index {idx} from {args.species_tree} (avg support: {avg_support:.4f})")
    except Exception as e:
        print(f"✗ Error processing species tree: {e}")
        return

    # Clean and save gene trees
    try:
        gene_trees, dropped_genes = clean_and_save_trees(str(gene_trees_path), gene_trees_clean)
        print(f"✓ Gene trees cleaned and saved to: {gene_trees_clean}")
        print(f"  Processed {len(gene_trees)} tree(s)")
        if dropped_genes:
            print(f"  ⚠ Dropped {len(dropped_genes)} tree(s) with avg support < 0.5:")
            for idx, avg_support in dropped_genes.items():
                print(f"    - Index {idx} from {args.gene_trees} (avg support: {avg_support:.4f})")
    except Exception as e:
        print(f"✗ Error processing gene trees: {e}")
        return

    print("\nProcessing complete!")


if __name__ == "__main__":
    main()
