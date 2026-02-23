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
from multiprocessing import cpu_count
from pathlib import Path
import time

import dendropy

from .tree_parser import (
    _build_species_triplet_metadata,
    _parse_outgroup_arg,
    _root_tree_on_outgroup,
    clean_and_save_gene_trees,
    clean_and_save_trees,
    filter_triplets_by_taxa,
    format_newick_with_precision,
    generate_triplets,
    get_taxa_from_tree,
    MetricsLogger,
    read_triplet_filter_file,
    read_tree_file,
    write_clean_trees,
    write_triplet_gene_trees_multiprocess,
)
from .triplet_processor import (
    analyze_triplet_gene_tree_file,
    write_pipeline_results,
)
from .config import ConfigError, load_orchestrator_config


DEFAULT_OUTPUT_FOLDER = "results"
DEFAULT_PROCESSES = 0
DEFAULT_MIN_SUPPORT_VALUE = 0.5
DEFAULT_DISCORDANT_TEST = "chi-square"
DEFAULT_SUMMARY_STATISTIC = "mean"
DEFAULT_ALPHA_DCT = 0.01
DEFAULT_ALPHA_KS = 0.05


def _default_output_path():
    """Default output directory path for CLI mode."""
    return str(Path.cwd() / DEFAULT_OUTPUT_FOLDER)


def _resolve_processes(processes):
    """Resolve process count where 0 means all available CPU cores."""
    if processes == 0:
        return cpu_count()
    return processes


def _resolve_parallel_mode(processes):
    """Resolve worker count and whether multiprocessing should be enabled."""
    resolved_processes = _resolve_processes(processes)
    use_multiprocessing = resolved_processes is not None and resolved_processes > 1
    return resolved_processes, use_multiprocessing


def _build_argument_parser():
    """Build the orchestrator CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="GhostParser orchestrator: run tree_parser and triplet_processor end-to-end."
    )

    parser.add_argument("-c", "--config-file", type=str, default=None, help="Path to a JSON or YAML config file")
    parser.add_argument("-st", "--species_tree", default=None, help="Path to the species tree file in Newick format")
    parser.add_argument("-gt", "--gene_trees", default=None, help="Path to the gene trees file in Newick format")
    parser.add_argument("-og", "--outgroup", default=None, help="Outgroup species identifier")
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
        default=_default_output_path(),
        help=f"Output folder (default: ./{DEFAULT_OUTPUT_FOLDER})",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=DEFAULT_PROCESSES,
        help="Number of worker processes for triplet extraction/processing (0 = all cores)",
    )
    return parser


def _cli_args_used_alongside_config(args):
    """Return names of non-config CLI args provided together with --config-file."""
    provided = []
    if args.species_tree is not None:
        provided.append("--species_tree")
    if args.gene_trees is not None:
        provided.append("--gene_trees")
    if args.outgroup is not None:
        provided.append("--outgroup")
    if args.triplet_filter is not None:
        provided.append("--triplet-filter")
    if args.output != _default_output_path():
        provided.append("--output")
    if args.processes != DEFAULT_PROCESSES:
        provided.append("--processes")
    return provided


def _resolve_runtime_args(args):
    """Resolve runtime arguments from either config-file mode or pure CLI mode."""
    if args.config_file:
        ignored_cli_args = _cli_args_used_alongside_config(args)
        if ignored_cli_args:
            print(
                "Warning: --config-file provided; CLI arguments not in config will be ignored: "
                + ", ".join(ignored_cli_args)
            )

        config = load_orchestrator_config(args.config_file)
        return argparse.Namespace(
            species_tree=config["species_tree"],
            gene_trees=config["gene_trees"],
            outgroup=config["outgroup"],
            triplet_filter=config.get("triplet_filter"),
            output=config.get("output") or _default_output_path(),
            processes=config.get("processes", DEFAULT_PROCESSES),
            min_support_value=config.get("min_support_value"),
            discordant_test=config.get("discordant_test") or DEFAULT_DISCORDANT_TEST,
            summary_statistic=config.get("summary_statistic") or DEFAULT_SUMMARY_STATISTIC,
            alpha_dct=(
                config.get("alpha_dct")
                if config.get("alpha_dct") is not None
                else DEFAULT_ALPHA_DCT
            ),
            alpha_ks=(
                config.get("alpha_ks")
                if config.get("alpha_ks") is not None
                else DEFAULT_ALPHA_KS
            ),
        )

    if not args.species_tree or not args.gene_trees or not args.outgroup:
        raise ValueError(
            "Missing required CLI arguments. Provide --species_tree, --gene_trees, and --outgroup, "
            "or use --config-file."
        )

    return argparse.Namespace(
        species_tree=args.species_tree,
        gene_trees=args.gene_trees,
        outgroup=args.outgroup,
        triplet_filter=args.triplet_filter,
        output=args.output,
        processes=args.processes,
        min_support_value=None,
        discordant_test=DEFAULT_DISCORDANT_TEST,
        summary_statistic=DEFAULT_SUMMARY_STATISTIC,
        alpha_dct=DEFAULT_ALPHA_DCT,
        alpha_ks=DEFAULT_ALPHA_KS,
    )


def main():
    """CLI entry point for orchestrated end-to-end run."""
    parser = _build_argument_parser()
    parsed_args = parser.parse_args()

    try:
        args = _resolve_runtime_args(parsed_args)
    except (ValueError, ConfigError) as exc:
        print(f"Error: {exc}")
        return

    species_tree_path = Path(args.species_tree)
    gene_trees_path = Path(args.gene_trees)

    if not species_tree_path.exists():
        print(f"Error: Species tree file not found: {args.species_tree}")
        return

    if not gene_trees_path.exists():
        print(f"Error: Gene trees file not found: {args.gene_trees}")
        return

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    species_tree_clean = str(output_dir / f"processed_{species_tree_path.name}")
    gene_trees_clean = str(output_dir / f"processed_{gene_trees_path.name}")
    metrics_filepath = str(output_dir / "metrics.txt")

    with MetricsLogger(metrics_filepath) as metrics:
        metrics.log(f"Processing species tree: {args.species_tree}")
        metrics.log(f"Processing gene trees: {args.gene_trees}")
        outgroup_taxa = _parse_outgroup_arg(args.outgroup)
        metrics.log(f"Outgroup: {', '.join(outgroup_taxa)}")
        metrics.log(f"Discordant count test: {args.discordant_test}")
        metrics.log(f"Summary statistic after KS: {args.summary_statistic}")
        metrics.log(f"DCT alpha: {args.alpha_dct}")
        metrics.log(f"KS alpha: {args.alpha_ks}")
        support_threshold = (
            args.min_support_value
            if args.min_support_value is not None
            else DEFAULT_MIN_SUPPORT_VALUE
        )
        metrics.log(f"Support threshold: {support_threshold}")
        metrics.log("")

        try:
            species_start = time.time()
            species_trees, dropped_species = clean_and_save_trees(
                str(species_tree_path),
                species_tree_clean,
                min_avg_support=support_threshold,
            )
            metrics.log(f"✓ Species tree cleaned and saved to: {species_tree_clean}")
            metrics.log(f"  Processed {len(species_trees)} tree(s)")
            if dropped_species:
                metrics.log(
                    f"  ⚠ Dropped {len(dropped_species)} tree(s) with avg support < {support_threshold}"
                )

            species_trees = read_tree_file(species_tree_clean)
            metrics.log(f"  Time taken: {time.time() - species_start:.2f}s")
        except Exception as exc:
            metrics.log(f"✗ Error processing species tree: {exc}")
            return

        try:
            triplets = []
            species_triplet_trees = {}
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
                min_avg_support=support_threshold,
            )
            metrics.log(f"\n✓ Gene trees cleaned and saved to: {gene_trees_clean}")
            metrics.log(f"  Processed {len(gene_trees)} tree(s)")
            metrics.log(f"  Rooted {rooted_count} tree(s) on outgroup taxa")
            if missing_root_indices:
                metrics.log(f"  Discarded {len(missing_root_indices)} gene tree(s) without outgroup taxa")
            if dropped_genes:
                metrics.log(
                    f"  ⚠ Dropped {len(dropped_genes)} tree(s) with avg support < {support_threshold}"
                )
            metrics.log(f"  Time taken: {time.time() - genes_start:.2f}s")
        except Exception as exc:
            metrics.log(f"✗ Error processing gene trees: {exc}")
            return

        processes, use_multiprocessing = _resolve_parallel_mode(args.processes)

        try:
            metrics.log("\n✓ Running staged triplet extraction + introgression inference...")
            infer_start = time.time()

            triplet_output_path = str(output_dir / "unique_triplets_gene_trees.txt")
            total_subtrees, triplets_with_trees, extraction_workers = write_triplet_gene_trees_multiprocess(
                triplets,
                gene_trees_clean,
                triplet_output_path,
                species_triplet_trees=species_triplet_trees,
                use_multiprocessing=use_multiprocessing,
                processes=processes,
            )

            metrics.log(f"✓ Triplet gene trees saved to: {triplet_output_path}")

            results = analyze_triplet_gene_tree_file(
                triplet_output_path,
                alpha_dct=args.alpha_dct,
                alpha_ks=args.alpha_ks,
                discordant_test=args.discordant_test,
                summary_statistic=args.summary_statistic,
                use_multiprocessing=use_multiprocessing,
                processes=processes,
            )
            final_tsv = str(output_dir / "orchestrator_triplet_results.tsv")
            write_pipeline_results(results, final_tsv)

            metrics.log(f"✓ Final orchestrated TSV saved to: {final_tsv}")
            metrics.log(f"  Workers used for extraction stage: {extraction_workers}")
            if use_multiprocessing:
                metrics.log("  Inference stage: multiprocessing enabled")
            else:
                metrics.log("  Inference stage: single-worker mode (no multiprocessing)")
            metrics.log(f"  Triplets analyzed: {len(results)}")
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
