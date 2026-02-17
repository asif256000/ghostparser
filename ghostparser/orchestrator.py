"""End-to-end orchestration for GhostParser.

This module orchestrates both existing stages:
1) tree preprocessing / triplet extraction (from ``tree_parser``)
2) per-triplet introgression inference (from ``triplet_processor``)

It accepts the same CLI options as ``tree_parser`` and runs per-triplet
processing in parallel. Per-triplet statistics dictionaries are accumulated and
written to a final TSV report.
"""

from __future__ import annotations

import argparse
import csv
from multiprocessing import cpu_count
from pathlib import Path
import time

import dendropy

from .tree_parser import (
    _build_species_triplet_metadata,
    _calculate_worker_count,
    _get_mp_context,
    _parse_outgroup_arg,
    _root_tree_on_outgroup,
    clean_and_save_gene_trees,
    clean_and_save_trees,
    extract_triplet_subtree,
    filter_triplets_by_taxa,
    format_newick_with_precision,
    generate_triplets,
    get_taxa_from_tree,
    MetricsLogger,
    read_triplet_filter_file,
    read_tree_file,
    write_clean_trees,
)
from .triplet_processor import (
    _species_topology_from_newick,
    run_triplet_pipeline,
)


_ORCH_GENE_TREE_NEWICKS = None


def _resolve_processes(processes):
    """Resolve process count where 0 means all available CPU cores."""
    if processes == 0:
        return cpu_count()
    return processes


def _init_triplet_worker(gene_tree_newicks):
    """Initializer for per-triplet worker processes."""
    global _ORCH_GENE_TREE_NEWICKS
    _ORCH_GENE_TREE_NEWICKS = gene_tree_newicks


def _extract_triplet_gene_tree_newicks(triplet):
    """Extract rooted gene-tree subtrees for one triplet from in-memory gene trees."""
    triplet_set = set(triplet)
    extracted = []

    for newick_str in _ORCH_GENE_TREE_NEWICKS:
        tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
        tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
        if not triplet_set.issubset(tree_taxa):
            continue

        subtree = extract_triplet_subtree(tree, triplet)
        if subtree:
            extracted.append(format_newick_with_precision(subtree))

    return extracted


def _process_triplet(args):
    """Worker entry for one triplet.

    Returns:
        Tuple of (row_dict, gene_tree_count, triplet, species_triplet_newick, triplet_gene_trees_or_none)
    """
    triplet, species_triplet_newick, alpha_dct, alpha_ks, include_gene_trees = args

    triplet_gene_trees = _extract_triplet_gene_tree_newicks(triplet)
    species_topology = _species_topology_from_newick(species_triplet_newick, triplet)

    result = run_triplet_pipeline(
        triplet,
        triplet_gene_trees,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        species_topology=species_topology,
    )
    row = result.to_dict()

    if include_gene_trees:
        return row, len(triplet_gene_trees), triplet, species_triplet_newick, triplet_gene_trees

    return row, len(triplet_gene_trees), triplet, species_triplet_newick, None


def _append_triplet_gene_trees_log(out_f, triplet, species_triplet_newick, triplet_gene_trees, is_first):
    """Append one triplet section to a log file in tree_parser-compatible format."""
    if not is_first:
        out_f.write("\n" + "=" * 60 + "\n")

    out_f.write(f"{','.join(triplet)}\t{len(triplet_gene_trees)}\t{species_triplet_newick}\n\n")
    for newick in triplet_gene_trees:
        out_f.write(newick + "\n")


def analyze_triplets_parallel(
    triplets,
    species_triplet_trees,
    gene_tree_newicks,
    use_multiprocessing=True,
    processes=None,
    alpha_dct=0.01,
    alpha_ks=0.05,
    log_triplet_gene_trees=False,
    triplet_log_filepath=None,
):
    """Analyze triplets in parallel and return per-triplet dictionaries.

    Args:
        triplets: List of normalized ABC triplets.
        species_triplet_trees: Dict triplet -> species-triplet subtree Newick.
        gene_tree_newicks: List of cleaned rooted gene tree Newick strings.
        use_multiprocessing: Enable multiprocessing.
        processes: Process count (0 means all cores via caller resolution).
        alpha_dct: DCT significance threshold.
        alpha_ks: KS significance threshold.
        log_triplet_gene_trees: Whether to log per-triplet extracted gene trees.
        triplet_log_filepath: Path for optional triplet gene-tree log file.

    Returns:
        Tuple of (rows, worker_count, total_subtrees, triplets_with_trees).
    """
    if not triplets:
        return [], 0, 0, 0

    resolved_processes = _resolve_processes(processes)
    worker_count = _calculate_worker_count(
        len(triplets),
        use_multiprocessing=use_multiprocessing,
        processes=resolved_processes,
    )

    args = [
        (
            triplet,
            species_triplet_trees[triplet],
            alpha_dct,
            alpha_ks,
            log_triplet_gene_trees,
        )
        for triplet in triplets
    ]

    rows = []
    total_subtrees = 0
    triplets_with_trees = 0

    log_handle = None
    if log_triplet_gene_trees and triplet_log_filepath:
        log_handle = open(triplet_log_filepath, "w")

    try:
        if use_multiprocessing and worker_count > 1:
            ctx = _get_mp_context()
            with ctx.Pool(
                processes=worker_count,
                initializer=_init_triplet_worker,
                initargs=(gene_tree_newicks,),
            ) as pool:
                iterator = pool.imap(_process_triplet, args)
                for index, payload in enumerate(iterator):
                    row, count, triplet, species_triplet_newick, triplet_gene_trees = payload
                    rows.append(row)
                    total_subtrees += count
                    if count > 0:
                        triplets_with_trees += 1

                    if log_handle is not None:
                        _append_triplet_gene_trees_log(
                            log_handle,
                            triplet,
                            species_triplet_newick,
                            triplet_gene_trees or [],
                            is_first=(index == 0),
                        )
        else:
            _init_triplet_worker(gene_tree_newicks)
            for index, arg in enumerate(args):
                row, count, triplet, species_triplet_newick, triplet_gene_trees = _process_triplet(arg)
                rows.append(row)
                total_subtrees += count
                if count > 0:
                    triplets_with_trees += 1

                if log_handle is not None:
                    _append_triplet_gene_trees_log(
                        log_handle,
                        triplet,
                        species_triplet_newick,
                        triplet_gene_trees or [],
                        is_first=(index == 0),
                    )
    finally:
        if log_handle is not None:
            log_handle.close()

    return rows, worker_count, total_subtrees, triplets_with_trees


def write_orchestrated_results_tsv(rows, output_filepath):
    """Write orchestrated per-triplet dictionary rows to TSV."""
    fieldnames = [
        "triplet",
        "abc_mapping",
        "species_topology",
        "con_topology",
        "con_topology_display",
        "dis1_topology",
        "dis1_topology_display",
        "dis2_topology",
        "dis2_topology_display",
        "top1_topology",
        "top2_topology",
        "top3_topology",
        "highest_freq_topologies",
        "concordant_diff",
        "n_topology_ab",
        "n_topology_bc",
        "n_topology_ac",
        "n_con",
        "n_dis1",
        "n_dis2",
        "dct_p_value",
        "dct_significant",
        "ks_statistic",
        "ks_p_value",
        "ks_significant",
        "median_con",
        "median_dis",
        "classification",
        "skipped_trees",
    ]

    with open(output_filepath, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            topology_counts = row.get("topology_counts", {})
            ranking = row.get("topology_frequency_ranking") or [None, None, None]
            normalized = {
                "triplet": ",".join(row.get("triplet", ())),
                "abc_mapping": row.get("abc_mapping"),
                "species_topology": row.get("species_topology"),
                "con_topology": row.get("con_topology"),
                "con_topology_display": row.get("con_topology_display"),
                "dis1_topology": row.get("dis1_topology"),
                "dis1_topology_display": row.get("dis1_topology_display"),
                "dis2_topology": row.get("dis2_topology"),
                "dis2_topology_display": row.get("dis2_topology_display"),
                "top1_topology": ranking[0],
                "top2_topology": ranking[1],
                "top3_topology": ranking[2],
                "highest_freq_topologies": ",".join(row.get("highest_freq_topologies", [])),
                "concordant_diff": row.get("concordant_diff"),
                "n_topology_ab": topology_counts.get("((A,B),C)"),
                "n_topology_bc": topology_counts.get("((B,C),A)"),
                "n_topology_ac": topology_counts.get("((A,C),B)"),
                "n_con": row.get("n_con"),
                "n_dis1": row.get("n_dis1"),
                "n_dis2": row.get("n_dis2"),
                "dct_p_value": row.get("dct_p_value"),
                "dct_significant": row.get("dct_significant"),
                "ks_statistic": row.get("ks_statistic"),
                "ks_p_value": row.get("ks_p_value"),
                "ks_significant": row.get("ks_significant"),
                "median_con": row.get("median_con"),
                "median_dis": row.get("median_dis"),
                "classification": row.get("classification"),
                "skipped_trees": row.get("skipped_trees"),
            }
            writer.writerow(normalized)


def main():
    """CLI entry point for orchestrated end-to-end run."""
    parser = argparse.ArgumentParser(
        description="GhostParser orchestrator: run tree_parser and triplet_processor end-to-end."
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
        help="Number of worker processes for triplet extraction/processing (0 = all cores)",
    )
    parser.add_argument(
        "--no-multiprocessing",
        action="store_true",
        help="Disable multiprocessing for triplet extraction/processing",
    )
    parser.add_argument(
        "--log-triplet-gene-trees",
        action="store_true",
        help="Log generated triplets and their extracted gene trees (debug mode)",
    )

    args = parser.parse_args()

    species_tree_path = Path(args.species_tree)
    gene_trees_path = Path(args.gene_trees)

    if not species_tree_path.exists():
        print(f"Error: Species tree file not found: {args.species_tree}")
        return

    if not gene_trees_path.exists():
        print(f"Error: Gene trees file not found: {args.gene_trees}")
        return

    if args.output:
        output_dir = species_tree_path.parent / args.output
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = species_tree_path.parent

    species_tree_clean = str(output_dir / f"processed_{species_tree_path.name}")
    gene_trees_clean = str(output_dir / f"processed_{gene_trees_path.name}")
    metrics_filepath = str(output_dir / "metrics.txt")

    with MetricsLogger(metrics_filepath) as metrics:
        metrics.log(f"Processing species tree: {args.species_tree}")
        metrics.log(f"Processing gene trees: {args.gene_trees}")
        outgroup_taxa = _parse_outgroup_arg(args.outgroup)
        metrics.log(f"Outgroup: {', '.join(outgroup_taxa)}")
        metrics.log("Support threshold: 0.5")
        metrics.log("")

        try:
            species_start = time.time()
            species_trees, dropped_species = clean_and_save_trees(str(species_tree_path), species_tree_clean)
            metrics.log(f"✓ Species tree cleaned and saved to: {species_tree_clean}")
            metrics.log(f"  Processed {len(species_trees)} tree(s)")
            if dropped_species:
                metrics.log(f"  ⚠ Dropped {len(dropped_species)} tree(s) with avg support < 0.5")

            species_trees = read_tree_file(species_tree_clean)
            metrics.log(f"  Time taken: {time.time() - species_start:.2f}s")
        except Exception as exc:
            metrics.log(f"✗ Error processing species tree: {exc}")
            return

        try:
            if species_trees:
                taxa = get_taxa_from_tree(species_trees[0])
                metrics.log(f"\n✓ Found {len(taxa)} taxa in species tree")

                pruned_tree, _excluded_taxa, missing_taxa, ingroup_taxa = _root_tree_on_outgroup(
                    species_trees[0], outgroup_taxa
                )

                if missing_taxa:
                    metrics.log(
                        f"⚠ Warning: Outgroup taxa not found in species tree: {', '.join(sorted(missing_taxa))}"
                    )

                if pruned_tree is None or not ingroup_taxa:
                    metrics.log("⚠ Warning: Unable to root and prune species tree on outgroups")
                    return

                species_trees = [pruned_tree]
                write_clean_trees(species_trees, species_tree_clean)

                if args.triplet_filter:
                    filter_path = Path(args.triplet_filter)
                    if not filter_path.exists():
                        metrics.log(f"✗ Error: Triplet filter file not found: {args.triplet_filter}")
                        return

                    raw_triplets, invalid_lines = read_triplet_filter_file(str(filter_path))
                    for line_number, raw in invalid_lines:
                        metrics.log(f"⚠ Warning: Skipping invalid triplet line {line_number} in {filter_path}: {raw}")

                    triplets, skipped_triplets = filter_triplets_by_taxa(raw_triplets, set(ingroup_taxa))
                    for triplet, missing in skipped_triplets:
                        metrics.log(
                            "⚠ Warning: Skipping triplet with missing taxa: "
                            f"{','.join(triplet)} (missing: {', '.join(missing)})"
                        )
                else:
                    triplets = generate_triplets(sorted(ingroup_taxa), [])
                    metrics.log(f"✓ Generated {len(triplets)} unique triplets")

                species_tree_newick = format_newick_with_precision(pruned_tree)
                species_dendro_tree = dendropy.Tree.get(
                    data=species_tree_newick,
                    schema="newick",
                    preserve_underscores=True,
                )

                triplets, species_triplet_trees, skipped_species_triplets = _build_species_triplet_metadata(
                    species_dendro_tree,
                    triplets,
                )

                if skipped_species_triplets:
                    metrics.log(
                        "⚠ Warning: Skipping triplets that could not be mapped on species tree: "
                        + "; ".join(",".join(triplet) for triplet in skipped_species_triplets)
                    )

                metrics.log(f"✓ Normalized {len(triplets)} triplets to A,B,C (A and B are sisters)")
        except Exception as exc:
            metrics.log(f"✗ Error generating triplets: {exc}")
            return

        try:
            genes_start = time.time()
            gene_trees, dropped_genes, rooted_count, missing_root_indices = clean_and_save_gene_trees(
                str(gene_trees_path),
                gene_trees_clean,
                outgroup_taxa,
            )
            metrics.log(f"\n✓ Gene trees cleaned and saved to: {gene_trees_clean}")
            metrics.log(f"  Processed {len(gene_trees)} tree(s)")
            metrics.log(f"  Rooted {rooted_count} tree(s) on outgroup taxa")
            if missing_root_indices:
                metrics.log(f"  Discarded {len(missing_root_indices)} gene tree(s) without outgroup taxa")
            if dropped_genes:
                metrics.log(f"  ⚠ Dropped {len(dropped_genes)} tree(s) with avg support < 0.5")
            metrics.log(f"  Time taken: {time.time() - genes_start:.2f}s")
        except Exception as exc:
            metrics.log(f"✗ Error processing gene trees: {exc}")
            return

        use_multiprocessing = not args.no_multiprocessing
        processes = _resolve_processes(args.processes)

        try:
            metrics.log("\n✓ Running triplet extraction + per-triplet introgression inference...")
            infer_start = time.time()

            gene_tree_newicks = [format_newick_with_precision(tree) for tree in gene_trees]
            triplet_log_path = str(output_dir / "unique_triplets_gene_trees.txt") if args.log_triplet_gene_trees else None

            rows, infer_workers, total_subtrees, triplets_with_trees = analyze_triplets_parallel(
                triplets,
                species_triplet_trees,
                gene_tree_newicks,
                use_multiprocessing=use_multiprocessing,
                processes=processes,
                log_triplet_gene_trees=args.log_triplet_gene_trees,
                triplet_log_filepath=triplet_log_path,
            )

            final_tsv = str(output_dir / "orchestrator_triplet_results.tsv")
            write_orchestrated_results_tsv(rows, final_tsv)

            if args.log_triplet_gene_trees:
                metrics.log(f"✓ Triplet gene trees log saved to: {triplet_log_path}")

            metrics.log(f"✓ Final orchestrated TSV saved to: {final_tsv}")
            metrics.log(f"  Workers used for inference: {infer_workers}")
            metrics.log(f"  Triplets analyzed: {len(rows)}")
            metrics.log(f"  Triplets with gene trees: {triplets_with_trees}")
            metrics.log(f"  Total subtrees extracted: {total_subtrees}")
            metrics.log(f"  Time taken: {time.time() - infer_start:.2f}s")
        except Exception as exc:
            metrics.log(f"✗ Error in triplet inference stage: {exc}")
            return

        metrics.log("\nProcessing complete!")
        metrics.log(f"\nMetrics saved to: {metrics_filepath}")


if __name__ == "__main__":
    main()
