"""Tree parsing and standardization module using BioPython."""

import argparse
from io import StringIO
from itertools import combinations
from multiprocessing import cpu_count
import multiprocessing as mp
from pathlib import Path
import time

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


def read_newick_lines(filepath):
    """Read Newick strings from a file (one tree per line)."""
    with open(filepath, "r") as f:
        return [line.strip() for line in f if line.strip()]


def _chunk_list(items, chunk_size):
    """Yield successive chunks from a list."""
    for i in range(0, len(items), chunk_size):
        yield items[i : i + chunk_size]


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
    return _extract_triplet_subtree(tree, triplet_taxa)


def _extract_triplet_subtree(tree, triplet_taxa, tree_taxa=None):
    """Internal helper to extract a subtree for a triplet with optional taxa cache."""
    # Get all terminal names in the tree
    if tree_taxa is None:
        tree_taxa = {terminal.name for terminal in tree.get_terminals()}

    # Check if all triplet taxa are present in the tree
    if not set(triplet_taxa).issubset(tree_taxa):
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


_TRIPLET_SET = None
_GENE_TREE_NEWICKS = None
_OUTPUT_PATH = None
_OUTPUT_LOCK = None


def _init_triplet_worker(triplet_set):
    """Initializer for multiprocessing workers."""
    global _TRIPLET_SET
    _TRIPLET_SET = triplet_set


def _init_triplet_writer(gene_tree_newicks, output_path, lock):
    """Initializer for multiprocessing writers."""
    global _GENE_TREE_NEWICKS, _OUTPUT_PATH, _OUTPUT_LOCK
    _GENE_TREE_NEWICKS = gene_tree_newicks
    _OUTPUT_PATH = output_path
    _OUTPUT_LOCK = lock


def _process_single_gene_tree(newick_str):
    """Process a single gene tree (as Newick string) for all relevant triplets."""
    tree = Phylo.read(StringIO(newick_str), "newick")
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}
    taxa_sorted = sorted(tree_taxa)
    results = []

    for triplet in combinations(taxa_sorted, 3):
        if triplet not in _TRIPLET_SET:
            continue
        subtree = _extract_triplet_subtree(tree, triplet, tree_taxa=tree_taxa)
        if subtree:
            newick_str = format_newick_with_precision(subtree)
            results.append((triplet, newick_str))

    return results


def _process_triplet_chunk(triplet_chunk):
    """Process a chunk of triplets against all gene trees and append results to file."""
    triplet_results = {triplet: [] for triplet in triplet_chunk}

    for newick_str in _GENE_TREE_NEWICKS:
        tree = Phylo.read(StringIO(newick_str), "newick")
        tree_taxa = {terminal.name for terminal in tree.get_terminals()}
        for triplet in triplet_chunk:
            if not set(triplet).issubset(tree_taxa):
                continue
            subtree = _extract_triplet_subtree(tree, triplet, tree_taxa=tree_taxa)
            if subtree:
                newick = format_newick_with_precision(subtree)
                triplet_results[triplet].append(newick)

    total_subtrees = 0
    triplets_with_trees = 0
    output_lines = []
    for i, triplet in enumerate(triplet_chunk):
        newick_trees = triplet_results[triplet]
        count = len(newick_trees)
        total_subtrees += count
        if count > 0:
            triplets_with_trees += 1

        output_lines.append(",".join(triplet) + f"\t{count}")
        output_lines.append("")
        for newick in newick_trees:
            output_lines.append(newick)
        if i < len(triplet_chunk) - 1:
            output_lines.append("")
            output_lines.append("=" * 60)
            output_lines.append("")

    output_block = "\n".join(output_lines) + "\n"

    if _OUTPUT_LOCK is not None:
        with _OUTPUT_LOCK:
            with open(_OUTPUT_PATH, "a") as f:
                f.write(output_block)
    else:
        with open(_OUTPUT_PATH, "a") as f:
            f.write(output_block)

    return total_subtrees, triplets_with_trees


def process_gene_trees_for_triplets(gene_trees, triplets, use_multiprocessing=True, processes=None, chunksize=None):
    """Process gene trees to extract subtrees for each triplet.

    Processes each gene tree against all relevant triplets.
    Uses multiprocessing to parallelize gene tree processing if enabled.

    Args:
        gene_trees: List of Bio.Phylo tree objects or Newick strings.
        triplets: List of triplet tuples to extract.
        use_multiprocessing: Whether to use multiprocessing (default: True).
                           If False, processes sequentially on single worker.
        processes: Number of worker processes when multiprocessing is enabled.
                   Defaults to cpu_count(). Only used if use_multiprocessing=True.
        chunksize: Chunk size for multiprocessing task distribution.
                   Only used if use_multiprocessing=True.

    Returns:
        A dictionary mapping triplets to lists of subtree Newick strings.
    """
    triplet_gene_trees = {triplet: [] for triplet in triplets}

    if not gene_trees or not triplets:
        return triplet_gene_trees

    triplet_set = set(triplets)

    # Calculate worker count based on multiprocessing flag
    if use_multiprocessing:
        available_cpus = cpu_count()
        worker_count = processes or available_cpus
        worker_count = max(1, min(worker_count, len(gene_trees)))
    else:
        worker_count = 1

    if use_multiprocessing and worker_count > 1:
        # Serialize trees to Newick to ensure compatibility across platforms
        if isinstance(gene_trees[0], str):
            gene_tree_newicks = gene_trees
        else:
            gene_tree_newicks = [format_newick_with_precision(tree, decimal_places=15) for tree in gene_trees]
        if chunksize is None:
            chunksize = max(1, len(gene_tree_newicks) // (worker_count * 4))

        ctx = mp.get_context("fork") if hasattr(mp, "get_context") else mp
        with ctx.Pool(processes=worker_count, initializer=_init_triplet_worker, initargs=(triplet_set,)) as pool:
            for result in pool.imap(_process_single_gene_tree, gene_tree_newicks, chunksize=chunksize):
                for triplet, newick in result:
                    triplet_gene_trees[triplet].append(newick)
    else:
        if isinstance(gene_trees[0], str):
            for newick_str in gene_trees:
                for triplet, newick in _process_single_gene_tree(newick_str):
                    triplet_gene_trees[triplet].append(newick)
        else:
            for gene_tree in gene_trees:
                tree_taxa = {terminal.name for terminal in gene_tree.get_terminals()}
                taxa_sorted = sorted(tree_taxa)
                for triplet in combinations(taxa_sorted, 3):
                    if triplet not in triplet_set:
                        continue
                    subtree = _extract_triplet_subtree(gene_tree, triplet, tree_taxa=tree_taxa)
                    if subtree:
                        # Convert to Newick string using our custom formatting
                        newick_str = format_newick_with_precision(subtree)
                        triplet_gene_trees[triplet].append(newick_str)

    return triplet_gene_trees


def write_triplet_gene_trees_multiprocess(
    triplets,
    gene_tree_newicks,
    output_filepath,
    use_multiprocessing=True,
    processes=None,
    chunksize=None,
):
    """Write triplet gene trees using multiprocessing without storing full results in memory.

    Distributes triplet chunks across workers for parallel processing.
    Each worker writes results directly to output file using a lock.

    Args:
        triplets: List of triplet tuples to process.
        gene_tree_newicks: List of Newick format tree strings.
        output_filepath: Output file path for results.
        use_multiprocessing: Whether to use multiprocessing (default: True).
                           If False, processes sequentially on single worker.
        processes: Number of worker processes when multiprocessing is enabled.
                   Defaults to cpu_count(). Only used if use_multiprocessing=True.
        chunksize: Chunk size for triplet distribution to workers.
                   Only used if use_multiprocessing=True.

    Returns:
        Tuple of (total_subtrees, triplets_with_trees, worker_count).
    """
    # Truncate output file before appending from workers
    with open(output_filepath, "w") as f:
        f.write("")

    if not triplets:
        return 0, 0, 0

    # Calculate worker count based on multiprocessing flag
    if use_multiprocessing:
        available_cpus = cpu_count()
        worker_count = processes or available_cpus
        worker_count = max(1, min(worker_count, len(triplets)))
    else:
        worker_count = 1

    if chunksize is None:
        chunksize = max(1, len(triplets) // (worker_count * 4))

    triplet_chunks = list(_chunk_list(triplets, chunksize))

    if use_multiprocessing and worker_count > 1:
        ctx = mp.get_context("fork") if hasattr(mp, "get_context") else mp
        manager = ctx.Manager()
        lock = manager.Lock()

        with ctx.Pool(
            processes=worker_count,
            initializer=_init_triplet_writer,
            initargs=(gene_tree_newicks, output_filepath, lock),
        ) as pool:
            totals = pool.map(_process_triplet_chunk, triplet_chunks)
    else:
        _init_triplet_writer(gene_tree_newicks, output_filepath, None)
        totals = [_process_triplet_chunk(chunk) for chunk in triplet_chunks]

    total_subtrees = sum(t[0] for t in totals)
    triplets_with_trees = sum(t[1] for t in totals)

    return total_subtrees, triplets_with_trees, worker_count


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

            # Write separator between triplets (except after the last one)
            if i < len(triplet_gene_trees) - 1:
                f.write("\n" + "=" * 60 + "\n")


class MetricsLogger:
    """Logger that writes to both stdout and a metrics file."""

    def __init__(self, metrics_filepath):
        self.metrics_filepath = metrics_filepath
        self.metrics_file = None
        self.lines = []

    def __enter__(self):
        """Context manager entry: open the metrics file."""
        self.metrics_file = open(self.metrics_filepath, "w")
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
    """Main entry point for the tree parser module."""
    parser = argparse.ArgumentParser(
        description="Ghost parser for identifying ghost introgressions in phylogenetic trees."
    )

    parser.add_argument("-st", "--species_tree", required=True, help="Path to the species tree file in Newick format")
    parser.add_argument("-gt", "--gene_trees", required=True, help="Path to the gene trees file in Newick format")
    parser.add_argument("-og", "--outgroup", required=True, help="Outgroup species identifier")
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

    # Get output filenames
    species_tree_clean = get_clean_filename(str(species_tree_path))
    gene_trees_clean = get_clean_filename(str(gene_trees_path))
    metrics_filepath = str(species_tree_path.parent / f"{species_tree_path.stem}_metrics.txt")

    # Initialize metrics logger as context manager
    with MetricsLogger(metrics_filepath) as metrics:
        metrics.log(f"Processing species tree: {args.species_tree}")
        metrics.log(f"Processing gene trees: {args.gene_trees}")
        metrics.log(f"Outgroup: {args.outgroup}")
        metrics.log("Support threshold: 0.5")
        metrics.log("")

        # Clean and save species tree, then re-read from *_clean file for downstream steps
        try:
            species_start = time.time()
            species_trees, dropped_species = clean_and_save_trees(str(species_tree_path), species_tree_clean)
            metrics.log(f"✓ Species tree cleaned and saved to: {species_tree_clean}")
            metrics.log(f"  Processed {len(species_trees)} tree(s)")
            if dropped_species:
                metrics.log(f"  ⚠ Dropped {len(dropped_species)} tree(s) with avg support < 0.5:")
                for idx, avg_support in dropped_species.items():
                    metrics.log(f"    - Index {idx} from {args.species_tree} (avg support: {avg_support:.4f})")

            # Always use the cleaned file for subsequent processing
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

        # Clean and save gene trees, then re-read from *_clean file for downstream steps
        try:
            genes_start = time.time()
            gene_trees, dropped_genes = clean_and_save_trees(str(gene_trees_path), gene_trees_clean)
            metrics.log(f"\n✓ Gene trees cleaned and saved to: {gene_trees_clean}")
            metrics.log(f"  Processed {len(gene_trees)} tree(s)")
            if dropped_genes:
                metrics.log(f"  ⚠ Dropped {len(dropped_genes)} tree(s) with avg support < 0.5:")
                for idx, avg_support in dropped_genes.items():
                    metrics.log(f"    - Index {idx} from {args.gene_trees} (avg support: {avg_support:.4f})")

            # Read cleaned Newick strings for downstream processing
            gene_trees = read_newick_lines(gene_trees_clean)
            metrics.log(f"  Time taken: {time.time() - genes_start:.2f}s")
        except Exception as e:
            metrics.log(f"✗ Error processing gene trees: {e}")
            return

        # Process gene trees for triplets
        try:
            metrics.log("\n✓ Processing gene trees for triplets...")
            triplet_start = time.time()

            # Generate output filename for triplet gene trees
            triplet_output_path = str(gene_trees_path.parent / f"{gene_trees_path.stem}_unique_triplet_gene_trees.txt")

            total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
                triplets,
                gene_trees,
                triplet_output_path,
                use_multiprocessing=not args.no_multiprocessing,
                processes=args.processes,
            )

            # Print statistics
            metrics.log(f"✓ Triplet gene trees saved to: {triplet_output_path}")
            metrics.log(f"  Workers used: {worker_count}")
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
