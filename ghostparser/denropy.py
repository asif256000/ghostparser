"""DendroPy-based tree parsing and standardization module."""

import argparse
from itertools import combinations
from multiprocessing import cpu_count
import multiprocessing as mp
from pathlib import Path
import time

import dendropy


def read_tree_file(filepath):
    """Read Newick trees from a file using DendroPy.

    Args:
        filepath: Path to the Newick tree file.

    Returns:
        A list of dendropy.Tree objects.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file contains invalid Newick format.
    """
    try:
        trees = dendropy.TreeList.get(path=filepath, schema="newick", preserve_underscores=True)
        if not trees:
            raise ValueError(f"Invalid Newick format in {filepath}")
        for idx, tree in enumerate(trees, start=1):
            if not tree.leaf_nodes():
                raise ValueError(f"Invalid Newick format in {filepath}: Tree {idx} has no terminal nodes")
        return list(trees)
    except FileNotFoundError:
        raise FileNotFoundError(f"Tree file not found: {filepath}")
    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Invalid Newick format in {filepath}: {e}")


def calculate_average_support(tree):
    """Calculate the average support value for a tree.

    Args:
        tree: A dendropy.Tree object.

    Returns:
        The average support value across all internal nodes, or None if no support values found.
    """
    support_values = []
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        if node.label is not None:
            try:
                support_values.append(float(node.label))
            except ValueError:
                continue

    if not support_values:
        return None

    return sum(support_values) / len(support_values)


def remove_support_values(tree):
    """Remove support values (internal node labels) from a tree."""
    for node in tree.preorder_node_iter():
        if not node.is_leaf():
            node.label = None
    return tree


def standardize_tree(tree):
    """Standardize a DendroPy tree by removing support values."""
    return remove_support_values(tree)


def format_newick_with_precision(tree, decimal_places=10):
    """Format a tree to Newick string with custom decimal precision for branch lengths."""

    def format_branch_length(branch_length):
        formatted = f"{branch_length:.{decimal_places}f}"
        if "." in formatted:
            formatted = formatted.rstrip("0").rstrip(".")
        return formatted

    def format_node(node):
        if node.is_leaf():
            result = node.taxon.label if node.taxon else ""
        else:
            children = [format_node(child) for child in node.child_node_iter()]
            result = "(" + ",".join(children) + ")"

        if node.edge_length is not None:
            result += f":{format_branch_length(node.edge_length)}"

        return result

    return format_node(tree.seed_node) + ";"


def write_clean_trees(trees, output_filepath, decimal_places=15):
    """Write standardized trees to a file using custom precision."""
    with open(output_filepath, "w") as f:
        for tree in trees:
            newick_str = format_newick_with_precision(tree, decimal_places)
            f.write(newick_str + "\n")


def clean_and_save_trees(input_filepath, output_filepath, min_avg_support=0.5):
    """Read, standardize, and save trees from a file using DendroPy.

    Returns:
        A tuple of (cleaned_trees_list, dropped_trees_info_dict).
    """
    trees = read_tree_file(input_filepath)

    dropped_trees = {}
    cleaned_trees = []

    for idx, tree in enumerate(trees, start=1):
        avg_support = calculate_average_support(tree)
        if avg_support is not None and avg_support < min_avg_support:
            dropped_trees[idx] = avg_support
            continue

        standardized = standardize_tree(tree)
        cleaned_trees.append(standardized)

    write_clean_trees(cleaned_trees, output_filepath)

    return cleaned_trees, dropped_trees


def get_clean_filename(filepath):
    """Generate a clean filename with _clean prefix."""
    path = Path(filepath)
    return str(path.parent / f"{path.stem}_clean{path.suffix}")


def get_taxa_from_tree(tree):
    """Extract all terminal taxa names from a phylogenetic tree."""
    taxa = [leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon]
    return sorted(taxa)


def generate_triplets(taxa_list, outgroup):
    """Generate all unique triplet combinations from taxa excluding the outgroup."""
    ingroup_taxa = [taxon for taxon in taxa_list if taxon != outgroup]
    return list(combinations(ingroup_taxa, 3))


def write_triplets_to_file(triplets, output_filepath):
    """Write triplets to a file, one triplet per line."""
    with open(output_filepath, "w") as f:
        for triplet in triplets:
            f.write(",".join(triplet) + "\n")


def extract_triplet_subtree(tree, triplet_taxa):
    """Extract a subtree containing only the specified triplet taxa.

    Returns:
        A new dendropy.Tree containing only the specified taxa, or None if not all taxa are present.
    """
    tree_taxa = {leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon}
    if not set(triplet_taxa).issubset(tree_taxa):
        return None

    subtree = tree.extract_tree_with_taxa_labels(triplet_taxa)
    return subtree


def process_gene_trees_for_triplets(gene_trees, triplets):
    """Process gene trees to extract subtrees for each triplet."""
    triplet_gene_trees = {triplet: [] for triplet in triplets}

    for gene_tree in gene_trees:
        for triplet in triplets:
            subtree = extract_triplet_subtree(gene_tree, triplet)
            if subtree:
                newick_str = format_newick_with_precision(subtree)
                triplet_gene_trees[triplet].append(newick_str)

    return triplet_gene_trees


def _chunk_list(lst, chunk_size):
    """Split a list into chunks of specified size."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i : i + chunk_size]


def _get_mp_context():
    """Get a multiprocessing context that avoids fork in multi-threaded processes."""
    if hasattr(mp, "get_context"):
        methods = mp.get_all_start_methods()
        if "forkserver" in methods:
            return mp.get_context("forkserver")
        if "spawn" in methods:
            return mp.get_context("spawn")
    return mp


_TRIPLET_SET = None


def _init_triplet_worker(triplet_set):
    """Initializer for multiprocessing workers."""
    global _TRIPLET_SET
    _TRIPLET_SET = triplet_set


def _process_single_gene_tree(newick_str):
    """Process a single gene tree (as Newick string) for all relevant triplets."""
    tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
    tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
    taxa_sorted = sorted(tree_taxa)
    results = []

    for triplet in combinations(taxa_sorted, 3):
        if triplet not in _TRIPLET_SET:
            continue
        subtree = extract_triplet_subtree(tree, triplet)
        if subtree:
            newick = format_newick_with_precision(subtree)
            results.append((triplet, newick))

    return results


def _calculate_worker_count(total_items, use_multiprocessing, processes):
    """Calculate worker count based on multiprocessing flag and available items."""
    if not use_multiprocessing:
        return 1
    available_cpus = cpu_count()
    worker_count = processes or available_cpus
    return max(1, min(worker_count, total_items))


def write_triplet_gene_trees_multiprocess(
    triplets,
    gene_tree_newicks,
    output_filepath,
    use_multiprocessing=True,
    processes=None,
    chunksize=None,
):
    """Write triplet gene trees using multiprocessing over gene trees.

    Parallelizes across gene trees to avoid being bounded by the number of triplets.

    Args:
        triplets: List of triplet tuples to process.
        gene_tree_newicks: List of Newick format tree strings.
        output_filepath: Output file path for results.
        use_multiprocessing: Whether to use multiprocessing (default: True).
                           If False, processes sequentially on single worker.
        processes: Number of worker processes when multiprocessing is enabled.
                   Defaults to cpu_count(). Only used if use_multiprocessing=True.
        chunksize: Chunk size for gene tree distribution to workers.
                   Only used if use_multiprocessing=True.

    Returns:
        Tuple of (total_subtrees, triplets_with_trees, worker_count).
    """
    # Truncate output file before writing results
    with open(output_filepath, "w") as f:
        f.write("")

    if not triplets or not gene_tree_newicks:
        return 0, 0, 0

    worker_count = _calculate_worker_count(
        len(gene_tree_newicks),
        use_multiprocessing=use_multiprocessing,
        processes=processes,
    )

    if chunksize is None:
        chunksize = max(1, len(gene_tree_newicks) // (worker_count * 4))

    triplet_set = set(triplets)
    triplet_gene_trees = {triplet: [] for triplet in triplets}

    if use_multiprocessing and worker_count > 1:
        ctx = _get_mp_context()
        with ctx.Pool(
            processes=worker_count,
            initializer=_init_triplet_worker,
            initargs=(triplet_set,),
        ) as pool:
            for result in pool.imap(_process_single_gene_tree, gene_tree_newicks, chunksize=chunksize):
                for triplet, newick in result:
                    triplet_gene_trees[triplet].append(newick)
    else:
        _init_triplet_worker(triplet_set)
        for newick_str in gene_tree_newicks:
            for triplet, newick in _process_single_gene_tree(newick_str):
                triplet_gene_trees[triplet].append(newick)

    write_triplet_gene_trees(triplet_gene_trees, output_filepath)

    total_subtrees = sum(len(trees) for trees in triplet_gene_trees.values())
    triplets_with_trees = sum(1 for trees in triplet_gene_trees.values() if trees)

    return total_subtrees, triplets_with_trees, worker_count


def write_triplet_gene_trees(triplet_gene_trees, output_filepath):
    """Write triplet gene trees to a file in the specified format."""
    with open(output_filepath, "w") as f:
        for i, (triplet, newick_trees) in enumerate(triplet_gene_trees.items()):
            triplet_str = ",".join(triplet)
            count = len(newick_trees)
            f.write(f"{triplet_str}\t{count}\n")
            f.write("\n")

            for newick in newick_trees:
                f.write(f"{newick}\n")

            if i < len(triplet_gene_trees) - 1:
                f.write("\n" + "=" * 60 + "\n")


class MetricsLogger:
    """Logger that writes to both stdout and a metrics file."""

    def __init__(self, filepath):
        """Initialize the metrics logger.

        Args:
            filepath: Path to the metrics file.
        """
        self.filepath = filepath
        self.metrics_file = None
        self.lines = []

    def __enter__(self):
        """Context manager entry: open the metrics file."""
        try:
            self.metrics_file = open(self.filepath, "w")
        except Exception as e:
            print(f"Warning: Could not open metrics file {self.filepath}: {e}")
            self.metrics_file = None
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit: close the metrics file."""
        if self.metrics_file:
            self.metrics_file.close()
        return False

    def log(self, message):
        """Print message to both stdout and metrics file."""
        print(message)
        if self.metrics_file:
            self.metrics_file.write(message + "\n")
            self.metrics_file.flush()
        self.lines.append(message)


def main():
    """Main entry point for the DendroPy-based tree parser module."""
    parser = argparse.ArgumentParser(
        description="DendroPy-based ghost parser for identifying ghost introgressions in phylogenetic trees."
    )

    parser.add_argument("-st", "--species_tree", required=True, help="Path to the species tree file in Newick format")
    parser.add_argument("-gt", "--gene_trees", required=True, help="Path to the gene trees file in Newick format")
    parser.add_argument("-og", "--outgroup", required=True, help="Outgroup species identifier")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output folder relative to input data folder (default: same folder as input data)",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=None,
        help="Number of worker processes for triplet extraction (default: cpu_count)",
    )
    parser.add_argument(
        "--no-multiprocessing",
        action="store_true",
        help="Disable multiprocessing for triplet extraction",
    )

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

    # Determine output directory
    if args.output:
        output_dir = species_tree_path.parent / args.output
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = species_tree_path.parent

    # Get output filenames
    species_tree_clean = str(output_dir / f"{species_tree_path.stem}_clean{species_tree_path.suffix}")
    gene_trees_clean = str(output_dir / f"{gene_trees_path.stem}_clean{gene_trees_path.suffix}")
    metrics_filepath = str(output_dir / f"{species_tree_path.stem}_dendropy_metrics.txt")

    # Initialize metrics logger as context manager
    with MetricsLogger(metrics_filepath) as metrics:
        metrics.log(f"Processing species tree: {args.species_tree}")
        metrics.log(f"Processing gene trees: {args.gene_trees}")
        metrics.log(f"Outgroup: {args.outgroup}")
        metrics.log("Support threshold: 0.5")
        metrics.log("")

        # Clean and save species tree
        try:
            species_start = time.time()
            species_trees, dropped_species = clean_and_save_trees(str(species_tree_path), species_tree_clean)
            metrics.log(f"✓ Species tree cleaned and saved to: {species_tree_clean}")
            metrics.log(f"  Processed {len(species_trees)} tree(s)")
            if dropped_species:
                metrics.log(f"  ⚠ Dropped {len(dropped_species)} tree(s) with avg support < 0.5:")
                for idx, avg_support in dropped_species.items():
                    metrics.log(f"    - Index {idx} from {args.species_tree} (avg support: {avg_support:.4f})")

            # Re-read the cleaned file for subsequent processing
            species_trees = read_tree_file(species_tree_clean)
            metrics.log(f"  Time taken: {time.time() - species_start:.2f}s")
        except Exception as e:
            metrics.log(f"✗ Error processing species tree: {e}")
            return

        # Generate triplets from species tree
        try:
            # Get taxa from the first species tree
            if species_trees:
                taxa = get_taxa_from_tree(species_trees[0])
                metrics.log(f"\n✓ Found {len(taxa)} taxa in species tree")

                # Check if outgroup exists in taxa
                if args.outgroup not in taxa:
                    metrics.log(f"⚠ Warning: Outgroup '{args.outgroup}' not found in species tree taxa")
                    metrics.log(f"  Available taxa: {', '.join(taxa)}")
                else:
                    # Generate triplets excluding outgroup
                    triplets = generate_triplets(taxa, args.outgroup)

                    metrics.log(f"✓ Generated {len(triplets)} unique triplets")
                    metrics.log(f"  ({len(taxa) - 1} taxa choose 3 = {len(triplets)} combinations)")
        except Exception as e:
            metrics.log(f"✗ Error generating triplets: {e}")
            return

        # Clean and save gene trees
        try:
            genes_start = time.time()
            gene_trees, dropped_genes = clean_and_save_trees(str(gene_trees_path), gene_trees_clean)
            metrics.log(f"\n✓ Gene trees cleaned and saved to: {gene_trees_clean}")
            metrics.log(f"  Processed {len(gene_trees)} tree(s)")
            if dropped_genes:
                metrics.log(f"  ⚠ Dropped {len(dropped_genes)} tree(s) with avg support < 0.5:")
                for idx, avg_support in dropped_genes.items():
                    metrics.log(f"    - Index {idx} from {args.gene_trees} (avg support: {avg_support:.4f})")

            # Re-read cleaned gene trees for downstream processing
            gene_trees = read_tree_file(gene_trees_clean)
            metrics.log(f"  Time taken: {time.time() - genes_start:.2f}s")
        except Exception as e:
            metrics.log(f"✗ Error processing gene trees: {e}")
            return

        # Process gene trees for triplets
        try:
            metrics.log("\n✓ Processing gene trees for triplets...")
            triplet_start = time.time()

            # Generate output filename for triplet gene trees
            triplet_output_path = str(output_dir / f"{gene_trees_path.stem}_unique_triplet_gene_trees.txt")

            # Convert gene trees to Newick strings for multiprocessing
            gene_tree_newicks = [format_newick_with_precision(tree, decimal_places=15) for tree in gene_trees]

            worker_count = _calculate_worker_count(
                len(gene_tree_newicks),
                use_multiprocessing=not args.no_multiprocessing,
                processes=args.processes,
            )
            metrics.log(f"  Workers used: {worker_count}")

            total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
                triplets,
                gene_tree_newicks,
                triplet_output_path,
                use_multiprocessing=not args.no_multiprocessing,
                processes=args.processes,
            )

            # Print statistics
            metrics.log(f"✓ Triplet gene trees saved to: {triplet_output_path}")
            metrics.log(f"  Total triplets: {len(triplets)}")
            metrics.log(f"  Triplets with gene trees: {triplets_with_trees}")
            metrics.log(f"  Total subtrees extracted: {total_subtrees}")
            if triplets_with_trees > 0:
                avg_trees_per_triplet = total_subtrees / triplets_with_trees
                metrics.log(f"  Average trees per triplet: {avg_trees_per_triplet:.2f}")
            metrics.log(f"  Time taken: {time.time() - triplet_start:.2f}s")
        except Exception as e:
            metrics.log(f"✗ Error processing triplet gene trees: {e}")
            return

        metrics.log("\nProcessing complete!")
        metrics.log(f"\nMetrics saved to: {metrics_filepath}")


if __name__ == "__main__":
    main()
