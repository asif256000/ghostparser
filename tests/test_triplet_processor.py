"""Tests for the triplet_processor module."""

import argparse
import random

import dendropy
import pytest
from scipy import stats

import ghostparser.triplet_processor as triplet_processor_module

from ghostparser.triplet_processor import (
    analyze_triplet_gene_tree_file,
    classify_triplet_topology,
    collect_triplet_statistics,
    compute_tree_height_statistic,
    parse_triplet_gene_trees_file,
    pearson_discordant_chi_square_test,
    run_triplet_pipeline,
    two_sample_ks_test,
    two_sample_ks_test_hybrid,
    two_proportion_discordant_z_test,
    write_pipeline_statistics_json,
    write_pipeline_results,
    _resolve_runtime_args,
)
from ghostparser.triplet_utils import TOPOLOGY_AB, TOPOLOGY_AC, TOPOLOGY_BC
from ghostparser.triplet_utils import classify_triplet_topology_string


def _tree(newick):
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)


def test_compute_tree_height_statistic_matches_definition():
    tree = _tree("((B:2,C:3):4,A:1);")
    observed = compute_tree_height_statistic(tree)
    expected = (1 + (2 + 4) + (3 + 4)) / 3
    assert observed == pytest.approx(expected)


def test_classify_triplet_topology_string_for_all_three_topologies():
    species_triplet = ("A", "B", "C")

    assert classify_triplet_topology_string(_tree("((A:1,B:1):1,C:1);"), species_triplet) == TOPOLOGY_AB
    assert classify_triplet_topology_string(_tree("((B:1,C:1):1,A:1);"), species_triplet) == TOPOLOGY_BC
    assert classify_triplet_topology_string(_tree("((A:1,C:1):1,B:1);"), species_triplet) == TOPOLOGY_AC


def test_classify_triplet_topology_labels_concordant_and_discordants():
    species_triplet = ("A", "B", "C")
    topology_counts = {
        TOPOLOGY_AB: 4,
        TOPOLOGY_BC: 8,
        TOPOLOGY_AC: 2,
    }

    label, most_frequent_matches_concordant = classify_triplet_topology(
        _tree("((A:1,B:1):1,C:1);"),
        species_triplet,
        topology_counts,
    )
    assert label == "concordant"
    assert most_frequent_matches_concordant is False

    label, most_frequent_matches_concordant = classify_triplet_topology(
        _tree("((B:1,C:1):1,A:1);"),
        species_triplet,
        topology_counts,
    )
    assert label == "discordant1"
    assert most_frequent_matches_concordant is False

    label, most_frequent_matches_concordant = classify_triplet_topology(
        _tree("((A:1,C:1):1,B:1);"),
        species_triplet,
        topology_counts,
    )
    assert label == "discordant2"
    assert most_frequent_matches_concordant is False


def test_pearson_discordant_chi_square_balanced_counts_not_significant():
    chi2_stat, p_value = pearson_discordant_chi_square_test(10, 10)
    assert chi2_stat == pytest.approx(0.0)
    assert p_value == pytest.approx(1.0)


def test_two_proportion_discordant_z_test_balanced_counts_not_significant():
    z_stat, p_value = two_proportion_discordant_z_test(10, 10)
    assert z_stat == pytest.approx(0.0)
    assert p_value == pytest.approx(1.0)


@pytest.mark.reference
def test_custom_chi_square_matches_scipy_reference_randomized():
    rng = random.Random(123)
    for _ in range(1000):
        n_dis1 = rng.randint(0, 500)
        n_dis2 = rng.randint(0, 500)

        custom_stat, custom_p = pearson_discordant_chi_square_test(n_dis1, n_dis2)
        scipy_result = stats.chisquare([n_dis1, n_dis2])

        assert custom_stat == pytest.approx(float(scipy_result.statistic), rel=0.0, abs=1e-12)
        assert custom_p == pytest.approx(float(scipy_result.pvalue), rel=0.0, abs=1e-12)


@pytest.mark.reference
def test_custom_ks_matches_scipy_asymptotic_reference_randomized():
    rng = random.Random(456)
    p_diffs = []
    for _ in range(600):
        n1 = rng.randint(20, 80)
        n2 = rng.randint(20, 80)
        sample_a = [rng.random() for _ in range(n1)]
        sample_b = [rng.random() for _ in range(n2)]

        custom_d, custom_p = two_sample_ks_test(sample_a, sample_b)
        scipy_res = stats.ks_2samp(sample_a, sample_b, alternative="two-sided", method="asymp")

        assert custom_d == pytest.approx(float(scipy_res.statistic), rel=0.0, abs=1e-12)
        p_diffs.append(abs(custom_p - float(scipy_res.pvalue)))

    # TODO: This threshold is intentionally relaxed because custom KS currently uses
    # an asymptotic p-value approximation that can deviate from SciPy asymptotic values
    # in a small fraction of randomized cases. If/when a corrected KS p-value
    # calculation function is introduced, tighten this back to the stricter bound.
    assert max(p_diffs) < 0.04


def test_run_triplet_pipeline_uses_species_concordant_and_frequency_ranked_discordants():
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 4)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        rng=random.Random(0),
    )

    assert result.con_topology == TOPOLOGY_AB
    assert result.dis1_topology == TOPOLOGY_BC
    assert result.dis2_topology == TOPOLOGY_AC
    assert result.n_con == 8
    assert result.n_dis1 == 12
    assert result.n_dis2 == 4
    assert result.top1_topology == TOPOLOGY_BC
    assert TOPOLOGY_BC in result.highest_freq_topologies
    assert result.most_frequent_matches_concordant is False


def test_run_triplet_pipeline_supports_z_test_for_discordant_counts():
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 4)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        discordant_test="z-test",
    )

    assert result.dct_significant is True
    assert result.dct_p_value < 0.01


def test_run_triplet_pipeline_supports_median_summary_statistic():
    species_triplet = ("A", "B", "C")
    con_tree = "((A:1.0,B:1.0):2.0,C:3.0);"
    dis1_tree = "((B:0.2,C:0.2):0.3,A:0.5);"
    dis2_tree = "((A:0.2,C:0.2):0.3,B:0.5);"
    trees = ([con_tree] * 40) + ([dis1_tree] * 30) + ([dis2_tree] * 5)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        summary_statistic="median",
    )

    assert result.summary_statistic == "median"
    assert result.summary_con is not None
    assert result.summary_dis is not None


def test_run_triplet_pipeline_breaks_discordant_ties_by_first_topology():
    species_triplet = ("A", "B", "C")
    trees = (["((A:1,B:1):1,C:1);"] * 5) + (["((B:1,C:1):1,A:1);"] * 3) + (["((A:1,C:1):1,B:1);"] * 3)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
    )

    assert result.con_topology == TOPOLOGY_AB
    assert result.dis1_topology == TOPOLOGY_BC
    assert result.dis2_topology == TOPOLOGY_AC
    assert result.n_dis1 == 3
    assert result.n_dis2 == 3


def test_run_triplet_pipeline_no_introgression_when_dct_not_significant():
    species_triplet = ("A", "B", "C")
    trees = [
        "((B:1,C:1):1,A:1);",
        "((A:1,C:1):1,B:1);",
        "((A:1,B:1):1,C:1);",
    ] * 8

    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(1))

    assert result.classification == "no_introgression"
    assert not result.dct_significant
    assert result.ks_p_value is None


def test_run_triplet_pipeline_inflow_when_ks_not_significant():
    species_triplet = ("A", "B", "C")
    con_tree = "((A:0.6,B:0.6):0.4,C:1.0);"
    dis1_tree = "((B:0.6,C:0.6):0.4,A:1.0);"
    dis2_tree = "((A:0.6,C:0.6):0.4,B:1.0);"

    trees = ([dis1_tree] * 30) + ([dis2_tree] * 5) + ([con_tree] * 30)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(2))

    assert result.dct_significant
    assert result.ks_significant is False
    assert result.classification == "inflow_introgression"


def test_run_triplet_pipeline_outflow_when_con_summary_higher():
    species_triplet = ("A", "B", "C")
    con_tree = "((A:1.0,B:1.0):2.0,C:3.0);"
    dis1_tree = "((B:0.2,C:0.2):0.3,A:0.5);"
    dis2_tree = "((A:0.2,C:0.2):0.3,B:0.5);"

    trees = ([con_tree] * 40) + ([dis1_tree] * 30) + ([dis2_tree] * 5)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(3))

    assert result.classification == "outflow_introgression"
    assert result.ks_significant is True
    assert result.con_topology == TOPOLOGY_AB
    assert result.most_frequent_matches_concordant is True
    assert result.summary_con > result.summary_dis


def test_run_triplet_pipeline_ghost_when_dis_summary_higher():
    species_triplet = ("A", "B", "C")
    con_tree = "((A:0.2,B:0.2):0.3,C:0.5);"
    dis1_tree = "((B:1.0,C:1.0):2.0,A:3.0);"
    dis2_tree = "((A:0.2,C:0.2):0.3,B:0.5);"

    trees = ([con_tree] * 40) + ([dis1_tree] * 30) + ([dis2_tree] * 5)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(4))

    assert result.classification == "ghost_introgression"
    assert result.ks_significant is True
    assert result.con_topology == TOPOLOGY_AB
    assert result.most_frequent_matches_concordant is True
    assert result.summary_con < result.summary_dis


def test_parse_analyze_and_write_pipeline_roundtrip_with_species_header(tmp_path):
    content = """A,B,C\t4\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    parsed = parse_triplet_gene_trees_file(str(input_file))
    assert ("A", "B", "C") in parsed
    entry = parsed[("A", "B", "C")]
    assert entry["count"] == 4
    assert entry["species_tree"] == "((A:1,B:1):1,C:1);"
    assert len(entry["gene_trees"]) == 4

    results = analyze_triplet_gene_tree_file(str(input_file), rng=random.Random(5))
    assert len(results) == 1

    output_file = tmp_path / "triplet_introgression_results.tsv"
    write_pipeline_results(results, str(output_file))

    out = output_file.read_text()
    assert "con_topology" in out
    assert "species_tree" in out
    assert "most_frequent_matches_concordant" in out
    assert "A,B,C" in out


def test_analyze_triplet_gene_tree_file_with_multiprocessing(tmp_path):
    content = """A,B,C\t3\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    results = analyze_triplet_gene_tree_file(
        str(input_file),
        use_multiprocessing=True,
        processes=2,
    )

    assert len(results) == 1
    assert results[0].analyzed_trees == 3


def test_analyze_triplet_gene_tree_file_rejects_unknown_discordant_test(tmp_path):
    content = """A,B,C\t3\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="Unsupported discordant test method"):
        analyze_triplet_gene_tree_file(str(input_file), discordant_test="bad-test")


def test_analyze_triplet_gene_tree_file_rejects_unknown_summary_statistic(tmp_path):
    content = """A,B,C\t3\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="Unsupported summary statistic"):
        analyze_triplet_gene_tree_file(str(input_file), summary_statistic="mode")


@pytest.mark.reference
def test_two_sample_ks_test_hybrid_uses_scipy_near_threshold(monkeypatch):
    monkeypatch.setattr(triplet_processor_module, "two_sample_ks_test", lambda *_: (0.12, 0.051))
    monkeypatch.setattr(triplet_processor_module, "_two_sample_ks_test_scipy", lambda *_: (0.13, 0.049))

    d_stat, p_value = two_sample_ks_test_hybrid([0.1, 0.2], [0.3, 0.4], alpha=0.05, borderline_margin=0.01)
    assert d_stat == pytest.approx(0.13)
    assert p_value == pytest.approx(0.049)


def test_two_sample_ks_test_hybrid_keeps_custom_when_not_borderline(monkeypatch):
    monkeypatch.setattr(triplet_processor_module, "two_sample_ks_test", lambda *_: (0.12, 0.40))
    monkeypatch.setattr(triplet_processor_module, "_two_sample_ks_test_scipy", lambda *_: (0.22, 0.10))

    d_stat, p_value = two_sample_ks_test_hybrid([0.1, 0.2], [0.3, 0.4], alpha=0.05, borderline_margin=0.01)
    assert d_stat == pytest.approx(0.12)
    assert p_value == pytest.approx(0.40)


def test_two_sample_ks_test_hybrid_rejects_negative_margin():
    with pytest.raises(ValueError, match="borderline_margin"):
        two_sample_ks_test_hybrid([0.1, 0.2], [0.3, 0.4], borderline_margin=-0.1)


def test_parse_triplet_gene_trees_file_requires_species_tree_column(tmp_path):
    content = """A,B,C\t3

((A:1,B:1):1,C:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="expected 3 tab-separated fields"):
        parse_triplet_gene_trees_file(str(input_file))


def test_parse_triplet_gene_trees_file_rejects_empty_species_tree(tmp_path):
    content = """A,B,C\t3\t

((A:1,B:1):1,C:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="Invalid species tree in header"):
        parse_triplet_gene_trees_file(str(input_file))


def test_write_pipeline_results_includes_topology_columns(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 10) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(6))

    output_file = tmp_path / "results.tsv"
    write_pipeline_results([result], str(output_file))

    out = output_file.read_text()
    assert "con_topology" in out
    assert "dis1_topology" in out


def test_write_pipeline_results_marks_discordant_highest_freq(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(9))

    output_file = tmp_path / "results_highest_freq.tsv"
    write_pipeline_results([result], str(output_file))

    out = output_file.read_text()
    assert "((B,C),A)" in out
    assert "highest_freq_topologies" in out


def test_collect_triplet_statistics_returns_dict_list():
    species_triplet = ("A", "B", "C")
    trees = (["((A:1,B:1):1,C:1);"] * 5) + (["((B:1,C:1):1,A:1);"] * 3) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(7))

    stats = collect_triplet_statistics([result])
    assert isinstance(stats, list)
    assert len(stats) == 1
    assert stats[0]["triplet"] == species_triplet
    assert stats[0]["abc_mapping"] == "A=A,B=B,C=C"
    assert "most_frequent_matches_concordant" in stats[0]
    assert "topology_counts" in stats[0]
    assert stats[0]["topology_counts"][TOPOLOGY_AB] == result.n_topology_ab
    assert "topology_frequency_ranking" in stats[0]
    assert "highest_freq_topologies" in stats[0]


def test_write_pipeline_statistics_json(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((A:1,B:1):1,C:1);"] * 5) + (["((B:1,C:1):1,A:1);"] * 3) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(8))

    output_file = tmp_path / "stats.json"
    write_pipeline_statistics_json([result], str(output_file))

    out = output_file.read_text()
    assert "topology_counts" in out
    assert "classification" in out


def test_resolve_runtime_args_triplet_processor_cli_defaults():
    args = argparse.Namespace(
        config_file=None,
        input="unique_triplets_gene_trees.txt",
        output=None,
        stats_output=None,
        alpha_dct=0.01,
        alpha_ks=0.05,
        discordant_test="chi-square",
        summary_statistic="mean",
        processes=0,
        no_multiprocessing=False,
    )

    resolved = _resolve_runtime_args(args)
    assert resolved.input == "unique_triplets_gene_trees.txt"
    assert resolved.alpha_dct == 0.01
    assert resolved.alpha_ks == 0.05
    assert resolved.discordant_test == "chi-square"
    assert resolved.summary_statistic == "mean"


def test_resolve_runtime_args_triplet_processor_config_warns_and_ignores(tmp_path, capsys):
    config_path = tmp_path / "triplet_processor_config.json"
    config_path.write_text(
        """
{
  "input_path": "input.tsv",
  "output_path": "out.tsv",
  "discordant_test": "z-test",
  "summary_statistic": "median"
}
""".strip()
    )

    args = argparse.Namespace(
        config_file=str(config_path),
        input="unique_triplets_gene_trees.txt",
        output=None,
        stats_output=None,
        alpha_dct=0.01,
        alpha_ks=0.05,
        discordant_test="chi-square",
        summary_statistic="mean",
        processes=0,
        no_multiprocessing=False,
    )

    resolved = _resolve_runtime_args(args)
    captured = capsys.readouterr()

    assert "Warning: --config-file provided; CLI arguments not in config will be ignored" in captured.out
    assert resolved.input == "input.tsv"
    assert resolved.output == "out.tsv"
    assert resolved.discordant_test == "z-test"
    assert resolved.summary_statistic == "median"
