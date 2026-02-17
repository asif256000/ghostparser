"""Triplet processing and GhostParser decision pipeline.

This module implements the sequential hypothesis-testing logic described in
GhostParser Figure 6 for rooted species triplets:

1. Classify each triplet gene tree as concordant/dis1/dis2.
2. Compute tree height statistic ``H(T)`` as average root-to-tip distance.
3. Run discordant count test (two-sided Z test, alpha=0.01).
4. If significant, run tree height test (two-sample KS, alpha=0.05).
5. If significant, compare medians to classify inflow vs ghost introgression.

Concordant topology is always the topology matching the rooted species triplet.
Frequency ranking is reported separately (including ``[highest freq]`` tags)
and does not redefine concordant/discordant roles in inference.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import random

import dendropy

from .triplet_utils import (
    ALL_TOPOLOGIES,
    TOPOLOGY_AB,
    TOPOLOGY_AC,
    TOPOLOGY_BC,
    classify_triplet_topology_string,
    rank_topologies_by_frequency,
)


Classification = str


@dataclass(frozen=True)
class TripletPipelineResult:
    """Result of running the GhostParser pipeline for one rooted species triplet."""

    triplet: tuple[str, str, str]
    species_topology: str
    con_topology: str
    dis1_topology: str
    dis2_topology: str
    top1_topology: str
    top2_topology: str
    top3_topology: str
    highest_freq_topologies: tuple[str, ...]
    concordant_diff: bool
    n_topology_ab: int
    n_topology_bc: int
    n_topology_ac: int
    n_con: int
    n_dis1: int
    n_dis2: int
    dct_p_value: float
    dct_significant: bool
    ks_p_value: float | None
    ks_statistic: float | None
    ks_significant: bool | None
    median_con: float | None
    median_dis: float | None
    classification: Classification
    skipped_trees: int = 0

    def to_dict(self):
        """Serialize all relevant triplet statistics to a dictionary."""
        a_taxon, b_taxon, c_taxon = self.triplet
        abc_mapping = f"A={a_taxon},B={b_taxon},C={c_taxon}"
        con_topology_display = f"{self.con_topology} [diff]" if self.concordant_diff else self.con_topology

        def _display(topology, base):
            tagged = base
            if topology in self.highest_freq_topologies:
                tagged = f"{tagged} [highest freq]"
            return tagged

        return {
            "triplet": self.triplet,
            "abc_mapping": abc_mapping,
            "species_topology": self.species_topology,
            "con_topology": self.con_topology,
            "con_topology_display": _display(self.con_topology, con_topology_display),
            "dis1_topology": self.dis1_topology,
            "dis2_topology": self.dis2_topology,
            "dis1_topology_display": _display(self.dis1_topology, self.dis1_topology),
            "dis2_topology_display": _display(self.dis2_topology, self.dis2_topology),
            "topology_frequency_ranking": [self.top1_topology, self.top2_topology, self.top3_topology],
            "highest_freq_topologies": list(self.highest_freq_topologies),
            "concordant_diff": self.concordant_diff,
            "topology_counts": {
                TOPOLOGY_AB: self.n_topology_ab,
                TOPOLOGY_BC: self.n_topology_bc,
                TOPOLOGY_AC: self.n_topology_ac,
            },
            "n_con": self.n_con,
            "n_dis1": self.n_dis1,
            "n_dis2": self.n_dis2,
            "dct_p_value": self.dct_p_value,
            "dct_significant": self.dct_significant,
            "ks_statistic": self.ks_statistic,
            "ks_p_value": self.ks_p_value,
            "ks_significant": self.ks_significant,
            "median_con": self.median_con,
            "median_dis": self.median_dis,
            "classification": self.classification,
            "skipped_trees": self.skipped_trees,
        }


def _distance_to_root(node):
    """Compute root-to-node distance using edge lengths (missing lengths treated as 0)."""
    distance = 0.0
    current = node
    while current is not None and current.parent_node is not None:
        edge_length = current.edge_length
        if edge_length is not None:
            distance += float(edge_length)
        current = current.parent_node
    return distance


def compute_tree_height_statistic(tree):
    """Compute H(T) as mean root-to-tip distance over the 3 triplet leaves.

    For a rooted triplet ``((X:b2,Y:b3):b4,Z:b1)``, this equals:

    ``H(T) = (b1 + b2 + b3 + 2*b4) / 3``.
    """
    leaves = [leaf for leaf in tree.leaf_node_iter() if leaf.taxon and leaf.taxon.label]
    if len(leaves) != 3:
        raise ValueError("Triplet tree must contain exactly 3 terminal taxa")

    distances = [_distance_to_root(leaf) for leaf in leaves]
    return sum(distances) / 3.0


def classify_triplet_topology(tree, species_triplet):
    """Classify a rooted triplet tree relative to ABC as ``con``, ``dis1``, ``dis2``.

    This helper keeps backward compatibility for callers/tests. It maps from
    canonical topology strings to short labels where:

    - ``con`` -> ``((A,B),C)``
    - ``dis1`` -> ``((B,C),A)``
    - ``dis2`` -> ``((A,C),B)``
    """
    topology = classify_triplet_topology_string(tree, species_triplet)
    if topology == TOPOLOGY_AB:
        return "con"
    if topology == TOPOLOGY_BC:
        return "dis1"
    if topology == TOPOLOGY_AC:
        return "dis2"
    raise ValueError("Unknown triplet topology")


def two_sided_discordant_z_test(n_dis1, n_dis2):
    """Two-sided Z test for discordant topology count imbalance."""
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    z_score = (n_dis1 - n_dis2) / math.sqrt(total)
    p_value = math.erfc(abs(z_score) / math.sqrt(2.0))
    return z_score, min(max(p_value, 0.0), 1.0)


def _ks_statistic(sample_a, sample_b):
    """Compute two-sample KS D statistic."""
    if not sample_a or not sample_b:
        return 0.0

    x = sorted(float(v) for v in sample_a)
    y = sorted(float(v) for v in sample_b)
    n_x = len(x)
    n_y = len(y)

    i = 0
    j = 0
    cdf_x = 0.0
    cdf_y = 0.0
    d_stat = 0.0

    while i < n_x and j < n_y:
        if x[i] < y[j]:
            value = x[i]
            while i < n_x and x[i] == value:
                i += 1
            cdf_x = i / n_x
        elif y[j] < x[i]:
            value = y[j]
            while j < n_y and y[j] == value:
                j += 1
            cdf_y = j / n_y
        else:
            value = x[i]
            while i < n_x and x[i] == value:
                i += 1
            while j < n_y and y[j] == value:
                j += 1
            cdf_x = i / n_x
            cdf_y = j / n_y

        d_stat = max(d_stat, abs(cdf_x - cdf_y))

    while i < n_x:
        value = x[i]
        while i < n_x and x[i] == value:
            i += 1
        cdf_x = i / n_x
        d_stat = max(d_stat, abs(cdf_x - cdf_y))

    while j < n_y:
        value = y[j]
        while j < n_y and y[j] == value:
            j += 1
        cdf_y = j / n_y
        d_stat = max(d_stat, abs(cdf_x - cdf_y))

    return d_stat


def _ks_asymptotic_p_value(d_stat, n_x, n_y, terms=200):
    """Asymptotic p-value approximation for two-sample KS test."""
    if n_x == 0 or n_y == 0:
        return 1.0
    if d_stat <= 0.0:
        return 1.0

    en = math.sqrt((n_x * n_y) / (n_x + n_y))
    lam = (en + 0.12 + 0.11 / en) * d_stat

    summation = 0.0
    for idx in range(1, terms + 1):
        term = math.exp(-2.0 * (idx**2) * (lam**2))
        summation += term if idx % 2 == 1 else -term
        if term < 1e-12:
            break

    p_value = 2.0 * summation
    return min(max(p_value, 0.0), 1.0)


def two_sample_ks_test(sample_a, sample_b):
    """Two-sample KS test returning (D, p-value)."""
    d_stat = _ks_statistic(sample_a, sample_b)
    p_value = _ks_asymptotic_p_value(d_stat, len(sample_a), len(sample_b))
    return d_stat, p_value


def _median(values):
    """Compute median of numeric iterable."""
    if not values:
        return None
    sorted_vals = sorted(float(v) for v in values)
    n = len(sorted_vals)
    mid = n // 2
    if n % 2 == 1:
        return sorted_vals[mid]
    return (sorted_vals[mid - 1] + sorted_vals[mid]) / 2.0


def _species_topology_from_newick(species_tree_newick, abc_triplet):
    """Determine species topology string for an ABC triplet."""
    if not species_tree_newick:
        return TOPOLOGY_AB

    tree = dendropy.Tree.get(data=species_tree_newick, schema="newick", preserve_underscores=True)
    return classify_triplet_topology_string(tree, abc_triplet)


def run_triplet_pipeline(
    species_triplet,
    triplet_gene_trees,
    alpha_dct=0.01,
    alpha_ks=0.05,
    species_topology=TOPOLOGY_AB,
    rng=None,
):
    """Run GhostParser Figure 6 pipeline for one rooted species triplet.

    Args:
        species_triplet: Tuple ``(A, B, C)`` where ``A`` and ``B`` are sisters
            in species tree convention.
        triplet_gene_trees: Iterable of rooted triplet gene tree Newick strings.
        alpha_dct: Significance threshold for DCT-like Z test.
        alpha_ks: Significance threshold for THT (KS test).
        species_topology: Species-tree topology string for this triplet.
        rng: Optional randomizer for tie-breaking topology ranks.

    Returns:
        ``TripletPipelineResult``.
    """
    randomizer = rng or random
    heights = {topology: [] for topology in ALL_TOPOLOGIES}
    skipped_trees = 0

    species_set = set(species_triplet)

    for newick_str in triplet_gene_trees:
        if not str(newick_str).strip():
            continue

        tree = dendropy.Tree.get(data=str(newick_str).strip(), schema="newick", preserve_underscores=True)
        labels = {leaf.taxon.label for leaf in tree.leaf_node_iter() if leaf.taxon and leaf.taxon.label}
        if labels != species_set:
            skipped_trees += 1
            continue

        try:
            topology = classify_triplet_topology_string(tree, species_triplet)
            tree_height = compute_tree_height_statistic(tree)
        except ValueError:
            skipped_trees += 1
            continue

        heights[topology].append(tree_height)

    topology_counts = {topology: len(heights[topology]) for topology in ALL_TOPOLOGIES}
    top1_topology, top2_topology, top3_topology = rank_topologies_by_frequency(topology_counts, rng=rng)

    con_topology = species_topology
    discordant_topologies = [topology for topology in ALL_TOPOLOGIES if topology != species_topology]
    count_a = topology_counts[discordant_topologies[0]]
    count_b = topology_counts[discordant_topologies[1]]
    if count_a > count_b:
        dis1_topology, dis2_topology = discordant_topologies[0], discordant_topologies[1]
    elif count_b > count_a:
        dis1_topology, dis2_topology = discordant_topologies[1], discordant_topologies[0]
    else:
        shuffled = list(discordant_topologies)
        randomizer.shuffle(shuffled)
        dis1_topology, dis2_topology = shuffled[0], shuffled[1]

    max_count = max(topology_counts.values()) if topology_counts else 0
    highest_freq_topologies = tuple(
        topology for topology in ALL_TOPOLOGIES if topology_counts[topology] == max_count
    )
    concordant_diff = con_topology not in highest_freq_topologies

    n_con = topology_counts[con_topology]
    n_dis1 = topology_counts[dis1_topology]
    n_dis2 = topology_counts[dis2_topology]

    _, dct_p_value = two_sided_discordant_z_test(n_dis1, n_dis2)
    dct_significant = dct_p_value <= alpha_dct

    if not dct_significant:
        return TripletPipelineResult(
            triplet=species_triplet,
            species_topology=species_topology,
            con_topology=con_topology,
            dis1_topology=dis1_topology,
            dis2_topology=dis2_topology,
            top1_topology=top1_topology,
            top2_topology=top2_topology,
            top3_topology=top3_topology,
            highest_freq_topologies=highest_freq_topologies,
            concordant_diff=concordant_diff,
            n_topology_ab=topology_counts[TOPOLOGY_AB],
            n_topology_bc=topology_counts[TOPOLOGY_BC],
            n_topology_ac=topology_counts[TOPOLOGY_AC],
            n_con=n_con,
            n_dis1=n_dis1,
            n_dis2=n_dis2,
            dct_p_value=dct_p_value,
            dct_significant=False,
            ks_p_value=None,
            ks_statistic=None,
            ks_significant=None,
            median_con=None,
            median_dis=None,
            classification="no_introgression",
            skipped_trees=skipped_trees,
        )

    ks_statistic, ks_p_value = two_sample_ks_test(heights[dis1_topology], heights[con_topology])
    ks_significant = ks_p_value <= alpha_ks

    if not ks_significant:
        return TripletPipelineResult(
            triplet=species_triplet,
            species_topology=species_topology,
            con_topology=con_topology,
            dis1_topology=dis1_topology,
            dis2_topology=dis2_topology,
            top1_topology=top1_topology,
            top2_topology=top2_topology,
            top3_topology=top3_topology,
            highest_freq_topologies=highest_freq_topologies,
            concordant_diff=concordant_diff,
            n_topology_ab=topology_counts[TOPOLOGY_AB],
            n_topology_bc=topology_counts[TOPOLOGY_BC],
            n_topology_ac=topology_counts[TOPOLOGY_AC],
            n_con=n_con,
            n_dis1=n_dis1,
            n_dis2=n_dis2,
            dct_p_value=dct_p_value,
            dct_significant=True,
            ks_p_value=ks_p_value,
            ks_statistic=ks_statistic,
            ks_significant=False,
            median_con=None,
            median_dis=None,
            classification="outflow_introgression",
            skipped_trees=skipped_trees,
        )

    median_con = _median(heights[con_topology])
    median_dis = _median(heights[dis1_topology])

    if median_con is None or median_dis is None:
        classification = "unresolved"
    elif median_con > median_dis:
        classification = "inflow_introgression"
    elif median_con < median_dis:
        classification = "ghost_introgression"
    else:
        classification = "unresolved"

    return TripletPipelineResult(
        triplet=species_triplet,
        species_topology=species_topology,
        con_topology=con_topology,
        dis1_topology=dis1_topology,
        dis2_topology=dis2_topology,
        top1_topology=top1_topology,
        top2_topology=top2_topology,
        top3_topology=top3_topology,
        highest_freq_topologies=highest_freq_topologies,
        concordant_diff=concordant_diff,
        n_topology_ab=topology_counts[TOPOLOGY_AB],
        n_topology_bc=topology_counts[TOPOLOGY_BC],
        n_topology_ac=topology_counts[TOPOLOGY_AC],
        n_con=n_con,
        n_dis1=n_dis1,
        n_dis2=n_dis2,
        dct_p_value=dct_p_value,
        dct_significant=True,
        ks_p_value=ks_p_value,
        ks_statistic=ks_statistic,
        ks_significant=True,
        median_con=median_con,
        median_dis=median_dis,
        classification=classification,
        skipped_trees=skipped_trees,
    )


def parse_triplet_gene_trees_file(filepath):
    """Parse triplet gene tree file produced by ``tree_parser``.

    Supports both header styles:
    - ``A,B,C<TAB>count``
    - ``A,B,C<TAB>count<TAB>species_triplet_newick``

    Returns:
        dict mapping triplet tuple ->
        ``{"count": int | None, "species_tree": str | None, "gene_trees": list[str]}``
    """
    triplet_map = {}
    current_triplet = None

    with open(filepath, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if set(line) == {"="}:
                current_triplet = None
                continue

            if "\t" in line and "," in line.split("\t", 1)[0]:
                parts = line.split("\t")
                triplet_text = parts[0].strip()
                taxa = tuple(part.strip() for part in triplet_text.split(",") if part.strip())
                if len(taxa) != 3:
                    raise ValueError(f"Invalid triplet header: {line}")

                count = None
                if len(parts) >= 2 and parts[1].strip():
                    try:
                        count = int(parts[1].strip())
                    except ValueError:
                        count = None

                species_tree = parts[2].strip() if len(parts) >= 3 and parts[2].strip() else None
                current_triplet = taxa
                triplet_map[current_triplet] = {
                    "count": count,
                    "species_tree": species_tree,
                    "gene_trees": [],
                }
                continue

            if line.endswith(";"):
                if current_triplet is None:
                    raise ValueError("Encountered tree line before any triplet header")
                triplet_map[current_triplet]["gene_trees"].append(line)

    return triplet_map


def analyze_triplet_gene_tree_file(filepath, alpha_dct=0.01, alpha_ks=0.05, rng=None):
    """Analyze all triplets from a triplet-gene-trees file."""
    triplet_map = parse_triplet_gene_trees_file(filepath)
    results = []
    for triplet, entry in triplet_map.items():
        species_topology = _species_topology_from_newick(entry.get("species_tree"), triplet)
        results.append(
            run_triplet_pipeline(
                triplet,
                entry["gene_trees"],
                alpha_dct=alpha_dct,
                alpha_ks=alpha_ks,
                species_topology=species_topology,
                rng=rng,
            )
        )
    return results


def collect_triplet_statistics(results):
    """Return triplet results as list of dictionaries for downstream use."""
    return [result.to_dict() for result in results]


def write_pipeline_statistics_json(results, output_filepath):
    """Write all triplet statistics as JSON list for downstream analysis."""
    data = collect_triplet_statistics(results)
    with open(output_filepath, "w") as out_f:
        json.dump(data, out_f, indent=2)


def write_pipeline_results(results, output_filepath):
    """Write pipeline results to TSV file."""
    header = [
        "triplet",
        "abc_mapping",
        "species_topology",
        "con_topology",
        "dis1_topology",
        "dis2_topology",
        "top1_topology",
        "top2_topology",
        "top3_topology",
        "highest_freq_topologies",
        "concordant_diff",
        "n_topology_ab",
        "n_topology_bc",
        "n_topology_ac",
        "classification",
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
        "skipped_trees",
    ]

    with open(output_filepath, "w") as out_f:
        out_f.write("\t".join(header) + "\n")
        for result in results:
            con_topology_display = (
                f"{result.con_topology} [diff]" if result.concordant_diff else result.con_topology
            )
            if result.con_topology in result.highest_freq_topologies:
                con_topology_display += " [highest freq]"

            dis1_topology_display = result.dis1_topology
            if result.dis1_topology in result.highest_freq_topologies:
                dis1_topology_display += " [highest freq]"

            dis2_topology_display = result.dis2_topology
            if result.dis2_topology in result.highest_freq_topologies:
                dis2_topology_display += " [highest freq]"

            row = [
                ",".join(result.triplet),
                f"A={result.triplet[0]},B={result.triplet[1]},C={result.triplet[2]}",
                result.species_topology,
                con_topology_display,
                dis1_topology_display,
                dis2_topology_display,
                result.top1_topology,
                result.top2_topology,
                result.top3_topology,
                ",".join(result.highest_freq_topologies),
                str(result.concordant_diff),
                str(result.n_topology_ab),
                str(result.n_topology_bc),
                str(result.n_topology_ac),
                result.classification,
                str(result.n_con),
                str(result.n_dis1),
                str(result.n_dis2),
                f"{result.dct_p_value:.12g}",
                str(result.dct_significant),
                "" if result.ks_statistic is None else f"{result.ks_statistic:.12g}",
                "" if result.ks_p_value is None else f"{result.ks_p_value:.12g}",
                "" if result.ks_significant is None else str(result.ks_significant),
                "" if result.median_con is None else f"{result.median_con:.12g}",
                "" if result.median_dis is None else f"{result.median_dis:.12g}",
                str(result.skipped_trees),
            ]
            out_f.write("\t".join(row) + "\n")


def main():
    """CLI entry point for triplet processing pipeline."""
    parser = argparse.ArgumentParser(description="Run GhostParser triplet processing pipeline (Fig. 6).")
    parser.add_argument("-i", "--input", required=True, help="Path to unique_triplets_gene_trees.txt")
    parser.add_argument("-o", "--output", default=None, help="Output TSV path (default: alongside input)")
    parser.add_argument(
        "--stats-output",
        default=None,
        help="Optional JSON output path for full per-triplet statistics",
    )
    parser.add_argument("--alpha-dct", type=float, default=0.01, help="DCT significance threshold (default: 0.01)")
    parser.add_argument("--alpha-ks", type=float, default=0.05, help="KS significance threshold (default: 0.05)")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else input_path.parent / "triplet_introgression_results.tsv"

    results = analyze_triplet_gene_tree_file(str(input_path), alpha_dct=args.alpha_dct, alpha_ks=args.alpha_ks)
    write_pipeline_results(results, str(output_path))

    stats_output = Path(args.stats_output) if args.stats_output else output_path.with_suffix(".json")
    write_pipeline_statistics_json(results, str(stats_output))

    print(f"Processed {len(results)} triplets")
    print(f"Results written to: {output_path}")
    print(f"Statistics JSON written to: {stats_output}")


if __name__ == "__main__":
    main()
