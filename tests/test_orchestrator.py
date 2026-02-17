"""Tests for orchestrator module."""

from ghostparser.orchestrator import (
    _resolve_processes,
    analyze_triplets_parallel,
    write_orchestrated_results_tsv,
)


def test_resolve_processes_zero_uses_all_cores(monkeypatch):
    monkeypatch.setattr("ghostparser.orchestrator.cpu_count", lambda: 16)
    assert _resolve_processes(0) == 16
    assert _resolve_processes(None) is None
    assert _resolve_processes(4) == 4


def test_analyze_triplets_parallel_returns_dict_rows_single_worker():
    triplets = [("A", "B", "C")]
    species_triplet_trees = {("A", "B", "C"): "((A:1,B:1):1,C:1);"}
    gene_tree_newicks = [
        "((A:1,B:1):1,C:1);",
        "((B:1,C:1):1,A:1);",
        "((A:1,C:1):1,B:1);",
    ]

    rows, workers, total_subtrees, triplets_with_trees = analyze_triplets_parallel(
        triplets,
        species_triplet_trees,
        gene_tree_newicks,
        use_multiprocessing=False,
        processes=1,
    )

    assert workers == 1
    assert len(rows) == 1
    assert total_subtrees == 3
    assert triplets_with_trees == 1
    row = rows[0]
    assert row["triplet"] == ("A", "B", "C")
    assert "classification" in row
    assert "topology_counts" in row


def test_analyze_triplets_parallel_logs_triplet_gene_trees(tmp_path):
    triplets = [("A", "B", "C")]
    species_triplet_trees = {("A", "B", "C"): "((A:1,B:1):1,C:1);"}
    gene_tree_newicks = [
        "((A:1,B:1):1,C:1);",
        "((B:1,C:1):1,A:1);",
    ]
    log_file = tmp_path / "unique_triplets_gene_trees.txt"

    rows, workers, total_subtrees, triplets_with_trees = analyze_triplets_parallel(
        triplets,
        species_triplet_trees,
        gene_tree_newicks,
        use_multiprocessing=False,
        processes=1,
        log_triplet_gene_trees=True,
        triplet_log_filepath=str(log_file),
    )

    assert workers == 1
    assert len(rows) == 1
    assert total_subtrees == 2
    assert triplets_with_trees == 1
    content = log_file.read_text()
    assert "A,B,C\t2\t((A:1,B:1):1,C:1);" in content


def test_write_orchestrated_results_tsv(tmp_path):
    rows = [
        {
            "triplet": ("A", "B", "C"),
            "abc_mapping": "A=A,B=B,C=C",
            "species_topology": "((A,B),C)",
            "con_topology": "((A,B),C)",
            "con_topology_display": "((A,B),C)",
            "dis1_topology": "((B,C),A)",
            "dis1_topology_display": "((B,C),A) [highest freq]",
            "dis2_topology": "((A,C),B)",
            "dis2_topology_display": "((A,C),B)",
            "topology_frequency_ranking": ["((B,C),A)", "((A,B),C)", "((A,C),B)"],
            "highest_freq_topologies": ["((B,C),A)"],
            "concordant_diff": True,
            "topology_counts": {"((A,B),C)": 8, "((B,C),A)": 12, "((A,C),B)": 4},
            "n_con": 8,
            "n_dis1": 12,
            "n_dis2": 4,
            "dct_p_value": 0.02,
            "dct_significant": False,
            "ks_statistic": None,
            "ks_p_value": None,
            "ks_significant": None,
            "median_con": None,
            "median_dis": None,
            "classification": "no_introgression",
            "skipped_trees": 0,
        }
    ]

    out_file = tmp_path / "orchestrated.tsv"
    write_orchestrated_results_tsv(rows, str(out_file))

    content = out_file.read_text()
    assert "triplet\tabc_mapping\tspecies_topology" in content
    assert "A,B,C" in content
    assert "A=A,B=B,C=C" in content
    assert "((B,C),A)" in content
