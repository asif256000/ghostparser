"""Tree parsing and standardization module using BioPython."""

import argparse
from io import StringIO
from itertools import combinations
from pathlib import Path

from Bio import Phylo


def read_tree_file(filepath):
    """Read Newick trees from a file.

    Args:
        filepath: Path to the Newick tree file.

    Returns:
        A list of Bio.Phylo tree objects.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file contains invalid Newick format.
    """
    try:
        trees = list(Phylo.parse(filepath, "newick"))
        if not trees:
            raise ValueError(f"Invalid Newick format in {filepath}")
        # Validate that each tree has at least one terminal (leaf node)
        for idx, tree in enumerate(trees, start=1):
            if not tree.get_terminals():
                raise ValueError(f"Invalid Newick format in {filepath}: Tree {idx} has no terminal nodes")
        return trees
    except FileNotFoundError:
        raise FileNotFoundError(f"Tree file not found: {filepath}")
    except ValueError:
        # Re-raise ValueError as-is (our own validation errors)
        raise
    except Exception as e:
        # Catch any parsing errors from BioPython
        raise ValueError(f"Invalid Newick format in {filepath}: {e}")


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


def get_taxa_from_tree(tree):
    """Extract all terminal taxa names from a phylogenetic tree.

    Args:
        tree: A Bio.Phylo tree object.

    Returns:
        A sorted list of terminal taxa names.
    """
    taxa = [terminal.name for terminal in tree.get_terminals()]
    return sorted(taxa)


def generate_triplets(taxa_list, outgroup):
    """Generate all unique triplet combinations from taxa excluding the outgroup.

    Args:
        taxa_list: List of all taxa names.
        outgroup: The outgroup taxon to exclude.

    Returns:
        A list of tuples, where each tuple contains 3 taxa names.
    """
    # Remove outgroup from taxa list
    ingroup_taxa = [taxon for taxon in taxa_list if taxon != outgroup]

    # Generate all possible triplet combinations
    triplets = list(combinations(ingroup_taxa, 3))

    return triplets


def write_triplets_to_file(triplets, output_filepath):
    """Write triplets to a file, one triplet per line.

    Args:
        triplets: List of triplet tuples.
        output_filepath: Path where the triplets will be written.
    """
    with open(output_filepath, "w") as f:
        for triplet in triplets:
            f.write(",".join(triplet) + "\n")


def extract_triplet_subtree(tree, triplet_taxa):
    """Extract a subtree containing only the specified triplet taxa.

    Args:
        tree: A Bio.Phylo tree object.
        triplet_taxa: List or tuple of 3 taxa names to extract.

    Returns:
        A new Bio.Phylo tree containing only the specified taxa, or None if not all taxa are present.
    """
    # Get all terminal names in the tree
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}

    # Check if all triplet taxa are present in the tree
    if not all(taxon in tree_taxa for taxon in triplet_taxa):
        return None

    # Create a copy of the tree to avoid modifying the original
    import copy

    subtree = copy.deepcopy(tree)

    # Get all terminals to remove (all except our triplet)
    terminals_to_remove = [t.name for t in subtree.get_terminals() if t.name not in triplet_taxa]

    # Prune away unwanted terminals
    for terminal_name in terminals_to_remove:
        # Find the terminal in the tree
        for terminal in subtree.get_terminals():
            if terminal.name == terminal_name:
                subtree.prune(terminal)
                break

    return subtree


def process_gene_trees_for_triplets(gene_trees, triplets):
    """Process gene trees to extract subtrees for each triplet.

    Args:
        gene_trees: List of Bio.Phylo tree objects.
        triplets: List of triplet tuples.

    Returns:
        A dictionary mapping triplets to lists of subtree Newick strings.
    """
    triplet_gene_trees = {triplet: [] for triplet in triplets}

    for gene_tree in gene_trees:
        for triplet in triplets:
            subtree = extract_triplet_subtree(gene_tree, triplet)
            if subtree:
                # Convert to Newick string using our custom formatting
                newick_str = format_newick_with_precision(subtree)
                triplet_gene_trees[triplet].append(newick_str)

    return triplet_gene_trees


def write_triplet_gene_trees(triplet_gene_trees, output_filepath):
    """Write triplet gene trees to a file in the specified format.

    Format:
    TaxonA,TaxonB,TaxonC<tab><count>
    <blank line>
    newick_tree_1
    newick_tree_2
    ...
    <blank line>
    <blank line>
    NextTriplet...

    Args:
        triplet_gene_trees: Dictionary mapping triplets to lists of Newick strings.
        output_filepath: Path where the results will be written.
    """
    with open(output_filepath, "w") as f:
        for i, (triplet, newick_trees) in enumerate(triplet_gene_trees.items()):
            # Write header line: triplet and count
            triplet_str = ",".join(triplet)
            count = len(newick_trees)
            f.write(f"{triplet_str}\t{count}\n")

            # Write blank line
            f.write("\n")

            # Write all Newick trees for this triplet
            for newick in newick_trees:
                f.write(f"{newick}\n")

            # Write two blank lines between triplets (except after the last one)
            if i < len(triplet_gene_trees) - 1:
                f.write("\n\n")


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

    # Get output filenames
    species_tree_clean = get_clean_filename(str(species_tree_path))
    gene_trees_clean = get_clean_filename(str(gene_trees_path))

    print(f"Processing species tree: {args.species_tree}")
    print(f"Processing gene trees: {args.gene_trees}")
    print(f"Outgroup: {args.outgroup}")
    print("Support threshold: 0.5")
    print()

    # Clean and save species tree, then re-read from *_clean file for downstream steps
    try:
        species_trees, dropped_species = clean_and_save_trees(str(species_tree_path), species_tree_clean)
        print(f"✓ Species tree cleaned and saved to: {species_tree_clean}")
        print(f"  Processed {len(species_trees)} tree(s)")
        if dropped_species:
            print(f"  ⚠ Dropped {len(dropped_species)} tree(s) with avg support < 0.5:")
            for idx, avg_support in dropped_species.items():
                print(f"    - Index {idx} from {args.species_tree} (avg support: {avg_support:.4f})")

        # Always use the cleaned file for subsequent processing
        species_trees = read_tree_file(species_tree_clean)
    except Exception as e:
        print(f"✗ Error processing species tree: {e}")
        return

    # Generate triplets from species tree
    try:
        # Get taxa from the first species tree
        if species_trees:
            taxa = get_taxa_from_tree(species_trees[0])
            print(f"\n✓ Found {len(taxa)} taxa in species tree")

            # Check if outgroup exists in taxa
            if args.outgroup not in taxa:
                print(f"⚠ Warning: Outgroup '{args.outgroup}' not found in species tree taxa")
                print(f"  Available taxa: {', '.join(taxa)}")
            else:
                # Generate triplets excluding outgroup
                triplets = generate_triplets(taxa, args.outgroup)

                print(f"✓ Generated {len(triplets)} unique triplets")
                print(f"  ({len(taxa) - 1} taxa choose 3 = {len(triplets)} combinations)")
    except Exception as e:
        print(f"✗ Error generating triplets: {e}")
        return

    # Clean and save gene trees, then re-read from *_clean file for downstream steps
    try:
        gene_trees, dropped_genes = clean_and_save_trees(str(gene_trees_path), gene_trees_clean)
        print(f"\n✓ Gene trees cleaned and saved to: {gene_trees_clean}")
        print(f"  Processed {len(gene_trees)} tree(s)")
        if dropped_genes:
            print(f"  ⚠ Dropped {len(dropped_genes)} tree(s) with avg support < 0.5:")
            for idx, avg_support in dropped_genes.items():
                print(f"    - Index {idx} from {args.gene_trees} (avg support: {avg_support:.4f})")

        # Always use the cleaned file for subsequent processing
        gene_trees = read_tree_file(gene_trees_clean)
    except Exception as e:
        print(f"✗ Error processing gene trees: {e}")
        return

    # Process gene trees for triplets
    try:
        print("\n✓ Processing gene trees for triplets...")
        triplet_gene_trees = process_gene_trees_for_triplets(gene_trees, triplets)

        # Generate output filename for triplet gene trees
        triplet_output_path = str(species_tree_path.parent / f"{species_tree_path.stem}_unique_triplet_gene_trees.txt")

        # Write results to file
        write_triplet_gene_trees(triplet_gene_trees, triplet_output_path)

        # Print statistics
        total_subtrees = sum(len(trees) for trees in triplet_gene_trees.values())
        triplets_with_trees = sum(1 for trees in triplet_gene_trees.values() if trees)

        print(f"✓ Triplet gene trees saved to: {triplet_output_path}")
        print(f"  Total triplets: {len(triplets)}")
        print(f"  Triplets with gene trees: {triplets_with_trees}")
        print(f"  Total subtrees extracted: {total_subtrees}")
        if triplets_with_trees > 0:
            avg_trees_per_triplet = total_subtrees / triplets_with_trees
            print(f"  Average trees per triplet: {avg_trees_per_triplet:.2f}")
    except Exception as e:
        print(f"✗ Error processing triplet gene trees: {e}")
        return

    print("\nProcessing complete!")


if __name__ == "__main__":
    main()
