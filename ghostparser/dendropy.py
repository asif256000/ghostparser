"""Tree parsing and standardization module."""

import argparse
from itertools import combinations
from multiprocessing import cpu_count
import multiprocessing as mp
from pathlib import Path
import time

import dendropy


def read_tree_file(filepath):
    """Read Newick trees from a file.

    Args:
        filepath: Path to the Newick tree file.

    Returns:
        A list of tree objects.

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


def read_triplet_filter_file(filepath):
    """Read triplet filter file with comma-separated taxa per line.

    Returns:
        Tuple of (triplets, invalid_lines) where invalid_lines is a list of
        (line_number, line_text) for lines that do not parse into 3 taxa.
    """
    triplets = []
    invalid_lines = []
    with open(filepath, "r") as f:
        for line_number, line in enumerate(f, start=1):
            raw = line.strip()
            if not raw:
                continue
            parts = [part.strip() for part in raw.split(",")]
            parts = [part for part in parts if part]
            if len(parts) != 3:
                invalid_lines.append((line_number, raw))
                continue
            triplets.append(tuple(parts))
    return triplets, invalid_lines


def filter_triplets_by_taxa(triplets, taxa_set):
    """Filter triplets to those fully contained in taxa_set.

    Returns:
        Tuple of (kept_triplets, skipped_triplets) where skipped_triplets is a
        list of (triplet, missing_taxa).
    """
    kept_triplets = []
    skipped_triplets = []
    for triplet in triplets:
        missing = sorted(set(triplet) - taxa_set)
        if missing:
            skipped_triplets.append((triplet, missing))
            continue
        kept_triplets.append(triplet)
    return kept_triplets, skipped_triplets


def calculate_average_support(tree):
    """Calculate the average support value for a tree.

    Args:
        tree: A tree object.

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
    """Standardize a tree by removing support values."""
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
    """Read, standardize, and save trees from a file.

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
    """Generate all unique triplet combinations from taxa excluding the outgroup(s)."""
    if isinstance(outgroup, str):
        outgroup_taxa = set(_parse_outgroup_arg(outgroup))
    else:
        outgroup_taxa = set(outgroup)
    ingroup_taxa = [taxon for taxon in taxa_list if taxon not in outgroup_taxa]
    return list(combinations(ingroup_taxa, 3))


def write_triplets_to_file(triplets, output_filepath):
    """Write triplets to a file, one triplet per line."""
    with open(output_filepath, "w") as f:
        for triplet in triplets:
            f.write(",".join(triplet) + "\n")


def extract_triplet_subtree(tree, triplet_taxa):
    """Extract a subtree containing only the specified triplet taxa.

    Returns:
        A new tree containing only the specified taxa, or None if not all taxa are present.
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


def _parse_outgroup_arg(outgroup_arg):
    """Parse comma-separated outgroup taxa into a list of names."""
    if isinstance(outgroup_arg, str):
        parts = [p.strip() for p in outgroup_arg.split(",")]
        return [p for p in parts if p]
    return [str(taxon).strip() for taxon in outgroup_arg if str(taxon).strip()]


def _root_tree_on_outgroup(tree, outgroup_taxa):
    """Root a tree on the MRCA of outgroup taxa and prune the outgroup clade.

    Returns:
        Tuple of (pruned_tree, excluded_taxa_set, missing_taxa_set, ingroup_taxa_set).
    """
    tree_taxa = {leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon}
    outgroup_set = set(outgroup_taxa)
    missing = outgroup_set - tree_taxa
    present = [taxon for taxon in outgroup_taxa if taxon in tree_taxa]

    if not present:
        return None, set(), missing, set()

    tree.is_rooted = True
    mrca = tree.mrca(taxon_labels=present)
    if mrca is None:
        return None, set(), missing, set()

    excluded_taxa = {leaf.taxon.label for leaf in mrca.leaf_nodes() if leaf.taxon}
    tree.reroot_at_node(mrca)
    ingroup_taxa = set(tree_taxa) - excluded_taxa
    if not ingroup_taxa:
        return None, excluded_taxa, missing, set()

    # Prune outgroup taxa from the tree
    tree.prune_taxa_with_labels(list(excluded_taxa))

    return tree, excluded_taxa, missing, ingroup_taxa


_GENE_TREES_PATH = None
_CHUNK_DIR = None


def _init_triplet_chunk_worker(gene_trees_path, chunk_dir):
    """Initializer for multiprocessing triplet chunk workers."""
    global _GENE_TREES_PATH, _CHUNK_DIR
    _GENE_TREES_PATH = gene_trees_path
    _CHUNK_DIR = chunk_dir


def _process_triplet_chunk_stream(args):
    """Process a chunk of triplets against all gene trees and write to a chunk file."""
    chunk_index, triplet_chunk = args
    triplet_results = {triplet: [] for triplet in triplet_chunk}

    with open(_GENE_TREES_PATH, "r") as gene_f:
        for line in gene_f:
            newick_str = line.strip()
            if not newick_str:
                continue
            tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
            tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
            for triplet in triplet_chunk:
                if not set(triplet).issubset(tree_taxa):
                    continue
                subtree = extract_triplet_subtree(tree, triplet)
                if subtree:
                    triplet_results[triplet].append(format_newick_with_precision(subtree))

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

    chunk_path = Path(_CHUNK_DIR) / f"chunk_{chunk_index:06d}.txt"
    with open(chunk_path, "w") as out_f:
        out_f.write("\n".join(output_lines))
        if output_lines:
            out_f.write("\n")

    return total_subtrees, triplets_with_trees


def _merge_files_with_separators(input_paths, output_path):
    """Merge input files into output file with separators between file contents."""
    wrote_any = False
    with open(output_path, "w") as out_f:
        for path in input_paths:
            content = Path(path).read_text()
            if not content:
                continue
            if wrote_any:
                out_f.write("\n" + "=" * 60 + "\n")
            out_f.write(content)
            if not content.endswith("\n"):
                out_f.write("\n")
            wrote_any = True


def _merge_chunk_files(chunk_dir, output_filepath, batch_size=1000):
    """Merge chunk files in batches to limit number of temp files."""
    chunk_paths = sorted(Path(chunk_dir).glob("chunk_*.txt"))
    aggregates = []
    remaining = []

    for batch_index in range(0, len(chunk_paths), batch_size):
        batch = chunk_paths[batch_index : batch_index + batch_size]
        if len(batch) == batch_size:
            agg_path = Path(chunk_dir) / f"aggregate_{batch_index // batch_size:06d}.txt"
            _merge_files_with_separators(batch, agg_path)
            for path in batch:
                path.unlink()
            aggregates.append(agg_path)
        else:
            remaining = batch

    final_parts = aggregates + remaining
    _merge_files_with_separators(final_parts, output_filepath)

    for path in remaining:
        path.unlink()
    for path in aggregates:
        path.unlink()


def _calculate_worker_count(total_items, use_multiprocessing, processes):
    """Calculate worker count based on multiprocessing flag and available items."""
    if not use_multiprocessing:
        return 1
    available_cpus = cpu_count()
    worker_count = processes or available_cpus
    return max(1, min(worker_count, total_items))


def write_triplet_gene_trees_multiprocess(
    triplets,
    gene_trees_filepath,
    output_filepath,
    use_multiprocessing=True,
    processes=None,
    chunksize=None,
):
    """Write triplet gene trees using multiprocessing over triplets.

    Each worker processes a chunk of triplets and streams over the gene trees file.
    Results are written to per-chunk files and merged in order to avoid high memory use.

    Args:
        triplets: List of triplet tuples to process.
        gene_trees_filepath: Path to the cleaned gene trees file.
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
    # Truncate output file before writing results
    with open(output_filepath, "w") as f:
        f.write("")

    if not triplets:
        return 0, 0, 0

    worker_count = _calculate_worker_count(
        len(triplets),
        use_multiprocessing=use_multiprocessing,
        processes=processes,
    )

    if chunksize is None:
        chunksize = max(1, len(triplets) // (worker_count * 4))

    triplet_chunks = list(_chunk_list(triplets, chunksize))
    chunk_dir = Path(output_filepath).with_suffix("").parent / f".{Path(output_filepath).stem}_chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    temp_gene_file = None
    if not isinstance(gene_trees_filepath, (str, Path)):
        temp_gene_file = chunk_dir / "gene_trees_tmp.txt"
        temp_gene_file.write_text("\n".join(gene_trees_filepath))
        gene_trees_filepath = str(temp_gene_file)

    args = [(idx, chunk) for idx, chunk in enumerate(triplet_chunks)]

    if use_multiprocessing and worker_count > 1:
        ctx = _get_mp_context()
        with ctx.Pool(
            processes=worker_count,
            initializer=_init_triplet_chunk_worker,
            initargs=(gene_trees_filepath, str(chunk_dir)),
        ) as pool:
            totals = pool.map(_process_triplet_chunk_stream, args)
    else:
        _init_triplet_chunk_worker(gene_trees_filepath, str(chunk_dir))
        totals = [_process_triplet_chunk_stream(item) for item in args]

    # Merge chunk files in batches to avoid too many temp files
    _merge_chunk_files(chunk_dir, output_filepath, batch_size=1000)
    if temp_gene_file and temp_gene_file.exists():
        temp_gene_file.unlink()
    chunk_dir.rmdir()

    total_subtrees = sum(t[0] for t in totals)
    triplets_with_trees = sum(t[1] for t in totals)

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


def write_triplet_gene_trees_streaming(triplets, gene_trees_filepath, output_filepath):
    """Write triplet gene trees by streaming one triplet at a time.

    This avoids keeping all triplet results in memory and only stores
    the current triplet's subtrees.

    Returns:
        Tuple of (total_subtrees, triplets_with_trees).
    """
    total_subtrees = 0
    triplets_with_trees = 0

    with open(output_filepath, "w") as out_f:
        for idx, triplet in enumerate(triplets):
            newick_trees = []
            triplet_set = set(triplet)

            with open(gene_trees_filepath, "r") as gene_f:
                for line in gene_f:
                    newick_str = line.strip()
                    if not newick_str:
                        continue
                    tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
                    tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
                    if not triplet_set.issubset(tree_taxa):
                        continue
                    subtree = extract_triplet_subtree(tree, triplet)
                    if subtree:
                        newick_trees.append(format_newick_with_precision(subtree))

            count = len(newick_trees)
            total_subtrees += count
            if count > 0:
                triplets_with_trees += 1

            out_f.write(",".join(triplet) + f"\t{count}\n")
            out_f.write("\n")
            for newick in newick_trees:
                out_f.write(f"{newick}\n")

            if idx < len(triplets) - 1:
                out_f.write("\n" + "=" * 60 + "\n")

    return total_subtrees, triplets_with_trees


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
    """Main entry point for the tree parser module."""
    parser = argparse.ArgumentParser(
        description="Ghost parser for identifying ghost introgressions in phylogenetic trees."
    )

    parser.add_argument("-st", "--species_tree", required=True, help="Path to the species tree file in Newick format")
    parser.add_argument("-gt", "--gene_trees", required=True, help="Path to the gene trees file in Newick format")
    parser.add_argument("-og", "--outgroup", required=True, help="Outgroup species identifier")
    parser.add_argument(
        "-tf",
        "--triplet-filter",
        type=str,
        default=None,
        help="Path to triplet filter file (comma-separated taxa per line)",
    )
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
    metrics_filepath = str(output_dir / f"{species_tree_path.stem}_metrics.txt")

    # Initialize metrics logger as context manager
    with MetricsLogger(metrics_filepath) as metrics:
        metrics.log(f"Processing species tree: {args.species_tree}")
        metrics.log(f"Processing gene trees: {args.gene_trees}")
        outgroup_taxa = _parse_outgroup_arg(args.outgroup)
        metrics.log(f"Outgroup: {', '.join(outgroup_taxa)}")
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

                pruned_tree, excluded_taxa, missing_taxa, ingroup_taxa = _root_tree_on_outgroup(
                    species_trees[0], outgroup_taxa
                )

                if missing_taxa:
                    metrics.log(
                        f"⚠ Warning: Outgroup taxa not found in species tree: {', '.join(sorted(missing_taxa))}"
                    )

                if pruned_tree is None or not ingroup_taxa:
                    metrics.log("⚠ Warning: Unable to root and prune species tree on outgroups")
                    metrics.log(f"  Available taxa: {', '.join(taxa)}")
                    return

                extra_excluded = sorted(excluded_taxa - set(outgroup_taxa))
                metrics.log(f"  Outgroup taxa excluded: {len(excluded_taxa)}")
                if extra_excluded:
                    metrics.log(
                        "⚠ Warning: Additional taxa excluded while rooting on outgroups: "
                        + ", ".join(extra_excluded)
                    )
                    metrics.log("  Excluded taxa: " + ", ".join(sorted(excluded_taxa)))

                # Persist pruned species tree for downstream use
                species_trees = [pruned_tree]
                write_clean_trees(species_trees, species_tree_clean)

                # Generate triplets from ingroup taxa only (or filter if provided)
                if args.triplet_filter:
                    filter_path = Path(args.triplet_filter)
                    if not filter_path.exists():
                        metrics.log(f"✗ Error: Triplet filter file not found: {args.triplet_filter}")
                        return

                    raw_triplets, invalid_lines = read_triplet_filter_file(str(filter_path))
                    for line_number, raw in invalid_lines:
                        metrics.log(
                            f"⚠ Warning: Skipping invalid triplet line {line_number} in {filter_path}: {raw}"
                        )

                    filtered_triplets, skipped_triplets = filter_triplets_by_taxa(raw_triplets, set(ingroup_taxa))
                    for triplet, missing in skipped_triplets:
                        metrics.log(
                            "⚠ Warning: Skipping triplet with missing taxa: "
                            f"{','.join(triplet)} (missing: {', '.join(missing)})"
                        )

                    seen = set()
                    triplets = []
                    for triplet in filtered_triplets:
                        if triplet not in seen:
                            seen.add(triplet)
                            triplets.append(triplet)

                    metrics.log(f"✓ Using {len(triplets)} filtered triplets from: {filter_path}")
                else:
                    triplets = generate_triplets(sorted(ingroup_taxa), [])
                    metrics.log(f"✓ Generated {len(triplets)} unique triplets")
                    metrics.log(f"  ({len(ingroup_taxa)} taxa choose 3 = {len(triplets)} combinations)")

                triplet_list_path = str(output_dir / f"{species_tree_path.stem}_triplets.txt")
                write_triplets_to_file(triplets, triplet_list_path)
                metrics.log(f"✓ Triplets saved to: {triplet_list_path}")
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

            worker_count = _calculate_worker_count(
                len(triplets),
                use_multiprocessing=not args.no_multiprocessing,
                processes=args.processes,
            )
            metrics.log(f"  Workers used: {worker_count}")

            total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
                triplets,
                gene_trees_clean,
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
