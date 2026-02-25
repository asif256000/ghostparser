"""Triplet processing and GhostParser decision pipeline.

This module implements the sequential hypothesis-testing logic described in
GhostParser Figure 6 for rooted species triplets:

1. Classify each triplet gene tree as concordant/dis1/dis2.
2. Compute tree height statistic ``H(T)`` as average root-to-tip distance.
3. Run discordant count test (Pearson chi-square test, alpha=0.01).
4. If significant, run tree height test (two-sample KS, alpha=0.05).
5. If significant, compare selected summary values to classify inflow vs ghost introgression.

Concordant topology is defined as the species-tree topology for the ABC triplet.
The two discordant topologies are frequency-ranked as ``dis1`` (more frequent)
and ``dis2`` (less frequent), with deterministic first-discordant tie-breaking.

The discordant count test is configurable: Pearson chi-square (default) or
two-proportion z-test.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
import multiprocessing as mp
from multiprocessing import cpu_count
from pathlib import Path

import dendropy
from scipy import stats

from .config import ConfigError, load_triplet_processor_config
from .triplet_utils import (
    ALL_TOPOLOGIES,
    TOPOLOGY_AB,
    TOPOLOGY_AC,
    TOPOLOGY_BC,
    classify_triplet_topology_string,
    rank_topologies_by_frequency,
)


Classification = str
SerializedTripletObservation = tuple[str, float]
DISCORDANT_TEST_CHOICES = ("chi-square", "z-test")
SUMMARY_STATISTIC_CHOICES = ("mean", "median")


@dataclass(frozen=True)
class TripletPipelineResult:
    """Result of running the GhostParser pipeline for one rooted species triplet."""

    triplet: tuple[str, str, str]
    species_tree: str | None
    species_topology: str
    con_topology: str
    dis1_topology: str
    dis2_topology: str
    top1_topology: str
    top2_topology: str
    top3_topology: str
    highest_freq_topologies: tuple[str, ...]
    most_frequent_matches_concordant: bool
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
    summary_statistic: str
    summary_con: float | None
    summary_dis: float | None
    classification: Classification
    analyzed_trees: int = 0

    def to_dict(self):
        """Serialize all relevant triplet statistics to a dictionary."""
        a_taxon, b_taxon, c_taxon = self.triplet
        abc_mapping = f"A={a_taxon},B={b_taxon},C={c_taxon}"

        return {
            "triplet": self.triplet,
            "abc_mapping": abc_mapping,
            "species_tree": self.species_tree,
            "con_topology": self.con_topology,
            "dis1_topology": self.dis1_topology,
            "dis2_topology": self.dis2_topology,
            "topology_frequency_ranking": [self.top1_topology, self.top2_topology, self.top3_topology],
            "highest_freq_topologies": list(self.highest_freq_topologies),
            "most_frequent_matches_concordant": self.most_frequent_matches_concordant,
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
            "summary_statistic": self.summary_statistic,
            "summary_con": self.summary_con,
            "summary_dis": self.summary_dis,
            "classification": self.classification,
            "analyzed_trees": self.analyzed_trees,
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


def classify_triplet_topology(
    tree,
    species_triplet,
    topology_counts,
    species_topology=TOPOLOGY_AB,
):
    """Classify topology as concordant/discordant1/discordant2.

    Discordant1 and discordant2 are frequency-ranked among discordant
    topologies. Returns both the per-tree label and whether concordant is
    most frequent overall.
    """
    _, dis1_topology, dis2_topology, most_frequent_matches_concordant = _resolve_topology_roles(
        topology_counts,
        species_topology,
    )

    topology = classify_triplet_topology_string(tree, species_triplet)
    if topology == species_topology:
        return "concordant", most_frequent_matches_concordant
    if topology == dis1_topology:
        return "discordant1", most_frequent_matches_concordant
    if topology == dis2_topology:
        return "discordant2", most_frequent_matches_concordant
    raise ValueError("Unknown triplet topology")


def _resolve_topology_roles(topology_counts, species_topology):
    """Resolve concordant/discordant roles and concordant frequency status."""
    con_topology = species_topology
    discordant_topologies = [topology for topology in ALL_TOPOLOGIES if topology != con_topology]
    if len(discordant_topologies) != 2:
        raise ValueError(f"Invalid species topology: {species_topology}")

    first_discordant, second_discordant = discordant_topologies
    first_count = int(topology_counts.get(first_discordant, 0))
    second_count = int(topology_counts.get(second_discordant, 0))

    if first_count == second_count:
        dis1_topology, dis2_topology = first_discordant, second_discordant
    elif first_count > second_count:
        dis1_topology, dis2_topology = first_discordant, second_discordant
    else:
        dis1_topology, dis2_topology = second_discordant, first_discordant

    n_con = int(topology_counts.get(con_topology, 0))
    n_dis1 = int(topology_counts.get(dis1_topology, 0))
    n_dis2 = int(topology_counts.get(dis2_topology, 0))
    most_frequent_matches_concordant = n_con >= n_dis1 and n_con >= n_dis2

    return con_topology, dis1_topology, dis2_topology, most_frequent_matches_concordant


def pearson_discordant_chi_square_test(n_dis1, n_dis2):
    """Custom Pearson chi-square test for discordant topology count imbalance.

    This matches the SciPy default setup for two categories with equal expected
    frequencies and computes the p-value analytically for df=1.
    """
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    expected = total / 2.0
    chi2_stat = ((n_dis1 - expected) ** 2) / expected + ((n_dis2 - expected) ** 2) / expected

    p_value = math.erfc(math.sqrt(chi2_stat / 2.0))
    return float(chi2_stat), float(p_value)


def _pearson_discordant_chi_square_test_scipy(n_dis1, n_dis2):
    """Backup SciPy Pearson chi-square implementation."""
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    result = stats.chisquare([n_dis1, n_dis2])
    return float(result.statistic), float(result.pvalue)


def two_proportion_discordant_z_test(n_dis1, n_dis2):
    """Two-proportion z-test for discordant count imbalance.

    Uses pooled standard error with two proportions defined over the same
    discordant-total denominator for ``dis1`` and ``dis2``.
    """
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    p_dis1 = n_dis1 / total
    p_dis2 = n_dis2 / total
    pooled = (n_dis1 + n_dis2) / (2 * total)
    standard_error = (pooled * (1.0 - pooled) * ((1.0 / total) + (1.0 / total))) ** 0.5
    if standard_error == 0.0:
        return 0.0, 1.0

    z_score = (p_dis1 - p_dis2) / standard_error
    p_value = 2.0 * stats.norm.sf(abs(z_score))
    return float(z_score), float(p_value)


def run_discordant_count_test(n_dis1, n_dis2, method="chi-square"):
    """Run selected discordant count test and return (statistic, p-value)."""
    if method == "chi-square":
        return pearson_discordant_chi_square_test(n_dis1, n_dis2)
    if method == "z-test":
        return two_proportion_discordant_z_test(n_dis1, n_dis2)
    raise ValueError(f"Unsupported discordant test method: {method}")


def two_sample_ks_test(sample_a, sample_b):
    """Custom two-sample KS test returning (D, p-value).

    Uses the standard two-sided asymptotic Kolmogorov approximation for p-value:
    Q_K(lambda) = 2 * sum_{j>=1} (-1)^(j-1) exp(-2*j^2*lambda^2),
    with finite-sample lambda correction.
    """
    if not sample_a or not sample_b:
        return 0.0, 1.0

    data_a = sorted(float(value) for value in sample_a)
    data_b = sorted(float(value) for value in sample_b)
    n1 = len(data_a)
    n2 = len(data_b)

    i = 0
    j = 0
    cdf_a = 0.0
    cdf_b = 0.0
    d_stat = 0.0

    while i < n1 and j < n2:
        a_val = data_a[i]
        b_val = data_b[j]
        if a_val <= b_val:
            while i < n1 and data_a[i] == a_val:
                i += 1
            cdf_a = i / n1
        if b_val <= a_val:
            while j < n2 and data_b[j] == b_val:
                j += 1
            cdf_b = j / n2
        d_stat = max(d_stat, abs(cdf_a - cdf_b))

    while i < n1:
        i += 1
        cdf_a = i / n1
        d_stat = max(d_stat, abs(cdf_a - cdf_b))

    while j < n2:
        j += 1
        cdf_b = j / n2
        d_stat = max(d_stat, abs(cdf_a - cdf_b))

    en = (n1 * n2) / (n1 + n2)
    if en <= 0:
        return float(d_stat), 1.0

    sqrt_en = math.sqrt(en)
    lam = (sqrt_en + 0.12 + 0.11 / sqrt_en) * d_stat

    if lam <= 0:
        p_value = 1.0
    else:
        series_sum = 0.0
        for k in range(1, 101):
            term = math.exp(-2.0 * (k**2) * (lam**2))
            if k % 2 == 1:
                series_sum += term
            else:
                series_sum -= term
            if term < 1e-12:
                break
        p_value = max(0.0, min(1.0, 2.0 * series_sum))

    return float(d_stat), float(p_value)


def _two_sample_ks_test_scipy(sample_a, sample_b):
    """Backup SciPy two-sample KS implementation."""
    if not sample_a or not sample_b:
        return 0.0, 1.0

    result = stats.ks_2samp(sample_a, sample_b, alternative="two-sided", method="auto")
    return float(result.statistic), float(result.pvalue)


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


def _mean(values):
    """Compute mean of numeric iterable."""
    if not values:
        return None
    values_float = [float(value) for value in values]
    return sum(values_float) / len(values_float)


def _species_topology_from_newick(species_tree_newick, abc_triplet):
    """Determine species topology string for an ABC triplet."""
    if not species_tree_newick:
        return TOPOLOGY_AB

    tree = dendropy.Tree.get(data=species_tree_newick, schema="newick", preserve_underscores=True)
    return classify_triplet_topology_string(tree, abc_triplet)


def _serialize_triplet_gene_trees(species_triplet, triplet_gene_trees):
    """Parse rooted triplet trees into lightweight (topology, height) observations."""
    observations: list[SerializedTripletObservation] = []
    species_set = set(species_triplet)

    for newick_str in triplet_gene_trees:
        if not str(newick_str).strip():
            continue

        tree = dendropy.Tree.get(data=str(newick_str).strip(), schema="newick", preserve_underscores=True)
        labels = {leaf.taxon.label for leaf in tree.leaf_node_iter() if leaf.taxon and leaf.taxon.label}
        if labels != species_set:
            continue

        try:
            topology = classify_triplet_topology_string(tree, species_triplet)
            tree_height = compute_tree_height_statistic(tree)
        except ValueError:
            continue

        observations.append((topology, tree_height))

    return observations


def _run_triplet_pipeline_from_observations(
    species_triplet,
    observations,
    alpha_dct=0.01,
    alpha_ks=0.05,
    discordant_test="chi-square",
    summary_statistic="mean",
    species_topology=TOPOLOGY_AB,
    species_tree_newick=None,
    rng=None,
):
    """Run GhostParser Figure 6 pipeline from lightweight serialized observations."""
    heights = {topology: [] for topology in ALL_TOPOLOGIES}
    for topology, tree_height in observations:
        heights[topology].append(tree_height)

    topology_counts = {topology: len(heights[topology]) for topology in ALL_TOPOLOGIES}
    analyzed_trees = sum(topology_counts.values())
    top1_topology, top2_topology, top3_topology = rank_topologies_by_frequency(topology_counts, rng=rng)

    con_topology, dis1_topology, dis2_topology, most_frequent_matches_concordant = _resolve_topology_roles(
        topology_counts,
        species_topology,
    )

    max_count = max(topology_counts.values())
    highest_freq_topologies = tuple(
        topology for topology in ALL_TOPOLOGIES if topology_counts[topology] == max_count
    )
    n_con = topology_counts[con_topology]
    n_dis1 = topology_counts[dis1_topology]
    n_dis2 = topology_counts[dis2_topology]

    if summary_statistic not in SUMMARY_STATISTIC_CHOICES:
        raise ValueError(
            f"Unsupported summary statistic: {summary_statistic}. "
            f"Choose one of: {', '.join(SUMMARY_STATISTIC_CHOICES)}"
        )

    _, dct_p_value = run_discordant_count_test(n_dis1, n_dis2, method=discordant_test)
    dct_significant = dct_p_value <= alpha_dct

    if not dct_significant:
        return TripletPipelineResult(
            triplet=species_triplet,
            species_tree=species_tree_newick,
            species_topology=species_topology,
            con_topology=con_topology,
            dis1_topology=dis1_topology,
            dis2_topology=dis2_topology,
            top1_topology=top1_topology,
            top2_topology=top2_topology,
            top3_topology=top3_topology,
            highest_freq_topologies=highest_freq_topologies,
            most_frequent_matches_concordant=most_frequent_matches_concordant,
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
            summary_statistic=summary_statistic,
            summary_con=None,
            summary_dis=None,
            classification="no_introgression",
            analyzed_trees=analyzed_trees,
        )

    ks_statistic, ks_p_value = two_sample_ks_test(heights[dis1_topology], heights[con_topology])
    ks_significant = ks_p_value <= alpha_ks

    if not ks_significant:
        return TripletPipelineResult(
            triplet=species_triplet,
            species_tree=species_tree_newick,
            species_topology=species_topology,
            con_topology=con_topology,
            dis1_topology=dis1_topology,
            dis2_topology=dis2_topology,
            top1_topology=top1_topology,
            top2_topology=top2_topology,
            top3_topology=top3_topology,
            highest_freq_topologies=highest_freq_topologies,
            most_frequent_matches_concordant=most_frequent_matches_concordant,
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
            summary_statistic=summary_statistic,
            summary_con=None,
            summary_dis=None,
            classification="inflow_introgression",
            analyzed_trees=analyzed_trees,
        )

    if summary_statistic == "mean":
        summary_con = _mean(heights[con_topology])
        summary_dis = _mean(heights[dis1_topology])
    else:
        summary_con = _median(heights[con_topology])
        summary_dis = _median(heights[dis1_topology])

    if summary_con is None or summary_dis is None:
        classification = "unresolved"
    elif summary_con > summary_dis:
        classification = "outflow_introgression"
    elif summary_con < summary_dis:
        classification = "ghost_introgression"
    else:
        classification = "unresolved"

    return TripletPipelineResult(
        triplet=species_triplet,
        species_tree=species_tree_newick,
        species_topology=species_topology,
        con_topology=con_topology,
        dis1_topology=dis1_topology,
        dis2_topology=dis2_topology,
        top1_topology=top1_topology,
        top2_topology=top2_topology,
        top3_topology=top3_topology,
        highest_freq_topologies=highest_freq_topologies,
        most_frequent_matches_concordant=most_frequent_matches_concordant,
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
        summary_statistic=summary_statistic,
        summary_con=summary_con,
        summary_dis=summary_dis,
        classification=classification,
        analyzed_trees=analyzed_trees,
    )


def run_triplet_pipeline(
    species_triplet,
    triplet_gene_trees,
    alpha_dct=0.01,
    alpha_ks=0.05,
    discordant_test="chi-square",
    summary_statistic="mean",
    species_topology=TOPOLOGY_AB,
    species_tree_newick=None,
    rng=None,
):
    """Run GhostParser Figure 6 pipeline for one rooted species triplet.

    Args:
        species_triplet: Tuple ``(A, B, C)`` where ``A`` and ``B`` are sisters
            in species tree convention.
        triplet_gene_trees: Iterable of rooted triplet gene tree Newick strings.
        alpha_dct: Significance threshold for DCT-like Z test.
        alpha_ks: Significance threshold for THT (KS test).
        discordant_test: Discordant count test method (`chi-square` or `z-test`).
        summary_statistic: Statistic used after KS for con/dis1 distributions (`mean` or `median`).
        species_topology: Species-tree topology string for this triplet.
        species_tree_newick: Species-tree triplet Newick string for this triplet.
        rng: Optional randomizer for tie-breaking topology ranks.

    Returns:
        ``TripletPipelineResult``.
    """
    observations = _serialize_triplet_gene_trees(species_triplet, triplet_gene_trees)
    return _run_triplet_pipeline_from_observations(
        species_triplet,
        observations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        species_topology=species_topology,
        species_tree_newick=species_tree_newick,
        rng=rng,
    )


def parse_triplet_gene_trees_file(filepath):
    """Parse triplet gene tree file produced by ``tree_parser``.

    Required header style:
    - ``A,B,C<TAB>count<TAB>species_triplet_newick``

    Returns:
        dict mapping triplet tuple ->
        ``{"count": int | None, "species_tree": str | None, "gene_trees": list[str]}``
    """
    triplet_map = {}
    current_triplet = None

    with open(filepath, "r") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n").rstrip("\r")

            if not line.strip():
                continue

            if set(line) == {"="}:
                current_triplet = None
                continue

            if "\t" in line and "," in line.split("\t", 1)[0]:
                parts = line.split("\t")
                if len(parts) != 3:
                    raise ValueError(f"Invalid triplet header format (expected 3 tab-separated fields): {line}")

                triplet_text = parts[0].strip()
                taxa = tuple(part.strip() for part in triplet_text.split(",") if part.strip())
                if len(taxa) != 3:
                    raise ValueError(f"Invalid triplet header: {line}")

                if not parts[1].strip():
                    raise ValueError(f"Invalid triplet count in header: {line}")
                try:
                    count = int(parts[1].strip())
                except ValueError as exc:
                    raise ValueError(f"Invalid triplet count in header: {line}") from exc

                if not parts[2].strip():
                    raise ValueError(f"Invalid species tree in header: {line}")
                species_tree = parts[2].strip()
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


def _get_mp_context():
    """Get a multiprocessing context that avoids fork in multi-threaded processes."""
    if hasattr(mp, "get_context"):
        methods = mp.get_all_start_methods()
        if "forkserver" in methods:
            return mp.get_context("forkserver")
        if "spawn" in methods:
            return mp.get_context("spawn")
    return mp


def _resolve_processes(processes):
    """Resolve process count where 0 means all available CPU cores."""
    if processes == 0:
        return cpu_count()
    return processes


def _analyze_triplet_entry(args):
    """Analyze one triplet map entry in a worker process."""
    triplet, entry, alpha_dct, alpha_ks, discordant_test, summary_statistic = args
    species_topology = _species_topology_from_newick(entry.get("species_tree"), triplet)
    observations = _serialize_triplet_gene_trees(triplet, entry["gene_trees"])
    return _run_triplet_pipeline_from_observations(
        triplet,
        observations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        species_topology=species_topology,
        species_tree_newick=entry.get("species_tree"),
    )


def analyze_triplet_gene_tree_file(
    filepath,
    alpha_dct=0.01,
    alpha_ks=0.05,
    discordant_test="chi-square",
    summary_statistic="mean",
    rng=None,
    use_multiprocessing=True,
    processes=None,
):
    """Analyze all triplets from a triplet-gene-trees file."""
    if discordant_test not in DISCORDANT_TEST_CHOICES:
        raise ValueError(
            f"Unsupported discordant test method: {discordant_test}. "
            f"Choose one of: {', '.join(DISCORDANT_TEST_CHOICES)}"
        )

    if summary_statistic not in SUMMARY_STATISTIC_CHOICES:
        raise ValueError(
            f"Unsupported summary statistic: {summary_statistic}. "
            f"Choose one of: {', '.join(SUMMARY_STATISTIC_CHOICES)}"
        )

    triplet_map = parse_triplet_gene_trees_file(filepath)
    items = list(triplet_map.items())
    if not items:
        return []

    resolved_processes = _resolve_processes(processes)
    worker_count = resolved_processes or cpu_count()
    worker_count = max(1, min(worker_count, len(items)))

    if use_multiprocessing and worker_count > 1 and rng is None:
        args = [
            (triplet, entry, alpha_dct, alpha_ks, discordant_test, summary_statistic)
            for triplet, entry in items
        ]
        chunksize = max(1, len(args) // (worker_count * 4))
        ctx = _get_mp_context()
        with ctx.Pool(processes=worker_count) as pool:
            return list(pool.imap(_analyze_triplet_entry, args, chunksize=chunksize))

    results = []
    for triplet, entry in items:
        species_topology = _species_topology_from_newick(entry.get("species_tree"), triplet)
        results.append(
            run_triplet_pipeline(
                triplet,
                entry["gene_trees"],
                alpha_dct=alpha_dct,
                alpha_ks=alpha_ks,
                discordant_test=discordant_test,
                summary_statistic=summary_statistic,
                species_topology=species_topology,
                species_tree_newick=entry.get("species_tree"),
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
        "species_tree",
        "con_topology",
        "dis1_topology",
        "dis2_topology",
        "highest_freq_topologies",
        "n_topology_ab",
        "n_topology_bc",
        "n_topology_ac",
        "n_con",
        "n_dis1",
        "n_dis2",
        "most_frequent_matches_concordant",
        "dct_p_value",
        "dct_significant",
        "ks_statistic",
        "ks_p_value",
        "ks_significant",
        "summary_statistic",
        "summary_con",
        "summary_dis",
        "classification",
        "analyzed_trees",
    ]

    with open(output_filepath, "w") as out_f:
        out_f.write("\t".join(header) + "\n")
        for result in results:
            row = [
                ",".join(result.triplet),
                f"A={result.triplet[0]},B={result.triplet[1]},C={result.triplet[2]}",
                "" if result.species_tree is None else result.species_tree,
                result.con_topology,
                result.dis1_topology,
                result.dis2_topology,
                ",".join(result.highest_freq_topologies),
                str(result.n_topology_ab),
                str(result.n_topology_bc),
                str(result.n_topology_ac),
                str(result.n_con),
                str(result.n_dis1),
                str(result.n_dis2),
                str(result.most_frequent_matches_concordant),
                f"{result.dct_p_value:.12g}",
                str(result.dct_significant),
                "" if result.ks_statistic is None else f"{result.ks_statistic:.12g}",
                "" if result.ks_p_value is None else f"{result.ks_p_value:.12g}",
                "" if result.ks_significant is None else str(result.ks_significant),
                result.summary_statistic,
                "" if result.summary_con is None else f"{result.summary_con:.12g}",
                "" if result.summary_dis is None else f"{result.summary_dis:.12g}",
                result.classification,
                str(result.analyzed_trees),
            ]
            out_f.write("\t".join(row) + "\n")


def main():
    """CLI entry point for triplet processing pipeline."""
    parser = _build_argument_parser()
    parsed_args = parser.parse_args()

    try:
        args = _resolve_runtime_args(parsed_args)
    except (ValueError, ConfigError) as exc:
        print(f"Error: {exc}")
        return

    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else input_path.parent / "triplet_introgression_results.tsv"

    results = analyze_triplet_gene_tree_file(
        str(input_path),
        alpha_dct=args.alpha_dct,
        alpha_ks=args.alpha_ks,
        discordant_test=args.discordant_test,
        summary_statistic=args.summary_statistic,
        use_multiprocessing=not args.no_multiprocessing,
        processes=args.processes,
    )
    write_pipeline_results(results, str(output_path))

    stats_output = Path(args.stats_output) if args.stats_output else output_path.with_suffix(".json")
    write_pipeline_statistics_json(results, str(stats_output))

    print(f"Processed {len(results)} triplets")
    print(f"Results written to: {output_path}")
    print(f"Statistics JSON written to: {stats_output}")


def _build_argument_parser():
    """Build triplet_processor CLI argument parser."""
    parser = argparse.ArgumentParser(description="Run GhostParser triplet processing pipeline (Fig. 6).")
    parser.add_argument("-c", "--config-file", default=None, help="Path to a JSON or YAML config file")
    parser.add_argument("-i", "--input", default=None, help="Path to unique_triplets_gene_trees.txt")
    parser.add_argument("-o", "--output", default=None, help="Output TSV path (default: alongside input)")
    parser.add_argument(
        "--stats-output",
        default=None,
        help="Optional JSON output path for full per-triplet statistics",
    )
    parser.add_argument("--alpha-dct", type=float, default=0.01, help="DCT significance threshold (default: 0.01)")
    parser.add_argument("--alpha-ks", type=float, default=0.05, help="KS significance threshold (default: 0.05)")
    parser.add_argument(
        "--discordant-test",
        choices=DISCORDANT_TEST_CHOICES,
        default="chi-square",
        help="Discordant count test to use (default: chi-square)",
    )
    parser.add_argument(
        "--summary-statistic",
        choices=SUMMARY_STATISTIC_CHOICES,
        default="mean",
        help="Statistic used for con/dis1 distributions after KS test (default: mean)",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=0,
        help="Number of worker processes for triplet analysis (0 = all cores)",
    )
    parser.add_argument(
        "--no-multiprocessing",
        action="store_true",
        help="Disable multiprocessing for triplet analysis",
    )
    return parser


def _cli_args_used_alongside_config(args):
    """Return names of non-config CLI args provided together with --config-file."""
    provided = []
    if args.input is not None:
        provided.append("--input")
    if args.output is not None:
        provided.append("--output")
    if args.stats_output is not None:
        provided.append("--stats-output")
    if args.alpha_dct != 0.01:
        provided.append("--alpha-dct")
    if args.alpha_ks != 0.05:
        provided.append("--alpha-ks")
    if args.discordant_test != "chi-square":
        provided.append("--discordant-test")
    if args.summary_statistic != "mean":
        provided.append("--summary-statistic")
    if args.processes != 0:
        provided.append("--processes")
    if args.no_multiprocessing:
        provided.append("--no-multiprocessing")
    return provided


def _resolve_runtime_args(args):
    """Resolve runtime arguments from config-file mode or plain CLI mode."""
    if args.config_file:
        ignored_cli_args = _cli_args_used_alongside_config(args)
        if ignored_cli_args:
            print(
                "Warning: --config-file provided; CLI arguments not in config will be ignored: "
                + ", ".join(ignored_cli_args)
            )

        config = load_triplet_processor_config(args.config_file)
        return argparse.Namespace(
            input=config["input"],
            output=config.get("output"),
            stats_output=config.get("stats_output"),
            alpha_dct=0.01 if config.get("alpha_dct") is None else config.get("alpha_dct"),
            alpha_ks=0.05 if config.get("alpha_ks") is None else config.get("alpha_ks"),
            discordant_test=config.get("discordant_test") or "chi-square",
            summary_statistic=config.get("summary_statistic") or "mean",
            processes=config.get("processes", 0),
            no_multiprocessing=bool(config.get("no_multiprocessing", False)),
        )

    if not args.input:
        raise ValueError("Missing required CLI argument --input, or use --config-file.")

    return argparse.Namespace(
        input=args.input,
        output=args.output,
        stats_output=args.stats_output,
        alpha_dct=args.alpha_dct,
        alpha_ks=args.alpha_ks,
        discordant_test=args.discordant_test,
        summary_statistic=args.summary_statistic,
        processes=args.processes,
        no_multiprocessing=args.no_multiprocessing,
    )


if __name__ == "__main__":
    main()
