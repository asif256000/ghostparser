"""Tests for configuration parsing."""

import json

import pytest

from ghostparser.config import (
    ConfigError,
    load_orchestrator_config,
    load_tree_parser_config,
    load_triplet_processor_config,
)


def test_load_orchestrator_config_json(tmp_path):
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroups": ["OutA", "OutB"],
                "output_folder": "out",
                "processes": 4,
                "triplet_filter": "triplets.txt",
                "min_support_value": 0.7,
                "discordant_test": "z-test",
                "summary_statistic": "median",
                "alpha_dct": 0.02,
                "alpha_ks": 0.1,
            }
        )
    )

    config = load_orchestrator_config(str(config_path))

    assert config["species_tree"] == "species.nwk"
    assert config["gene_trees"] == "genes.nwk"
    assert config["outgroup"] == ["OutA", "OutB"]
    assert config["output"] == "out"
    assert config["processes"] == 4
    assert config["triplet_filter"] == "triplets.txt"
    assert config["min_support_value"] == 0.7
    assert config["discordant_test"] == "z-test"
    assert config["summary_statistic"] == "median"
    assert config["alpha_dct"] == 0.02
    assert config["alpha_ks"] == 0.1


def test_load_orchestrator_config_yaml(tmp_path):
    yaml = pytest.importorskip("yaml")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA,OutB",
            }
        )
    )

    config = load_orchestrator_config(str(config_path))
    assert config["outgroup"] == ["OutA", "OutB"]


def test_load_orchestrator_config_missing_required(tmp_path):
    config_path = tmp_path / "bad.json"
    config_path.write_text(json.dumps({"species_tree_path": "species.nwk"}))

    with pytest.raises(ConfigError, match="Missing required config field"):
        load_orchestrator_config(str(config_path))


def test_load_orchestrator_config_invalid_discordant_test(tmp_path):
    config_path = tmp_path / "bad_test_method.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
                "discordant_test": "invalid",
            }
        )
    )

    with pytest.raises(ConfigError, match="discordant_test"):
        load_orchestrator_config(str(config_path))


def test_load_orchestrator_config_invalid_summary_statistic(tmp_path):
    config_path = tmp_path / "bad_summary.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
                "summary_statistic": "mode",
            }
        )
    )

    with pytest.raises(ConfigError, match="summary_statistic"):
        load_orchestrator_config(str(config_path))


def test_load_tree_parser_config_json(tmp_path):
    config_path = tmp_path / "tree_parser_config.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroups": ["OutA", "OutB"],
                "output_folder": "out",
                "processes": 2,
                "triplet_filter": "triplets.txt",
                "min_support_value": 0.6,
                "no_multiprocessing": True,
            }
        )
    )

    config = load_tree_parser_config(str(config_path))

    assert config["species_tree"] == "species.nwk"
    assert config["gene_trees"] == "genes.nwk"
    assert config["outgroup"] == ["OutA", "OutB"]
    assert config["output"] == "out"
    assert config["processes"] == 2
    assert config["triplet_filter"] == "triplets.txt"
    assert config["min_support_value"] == 0.6
    assert config["no_multiprocessing"] is True


def test_load_tree_parser_config_invalid_no_multiprocessing(tmp_path):
    config_path = tmp_path / "tree_parser_bad.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
                "no_multiprocessing": "yes",
            }
        )
    )

    with pytest.raises(ConfigError, match="no_multiprocessing"):
        load_tree_parser_config(str(config_path))


def test_load_triplet_processor_config_json(tmp_path):
    config_path = tmp_path / "triplet_processor_config.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "unique_triplets_gene_trees.txt",
                "output_path": "results.tsv",
                "stats_output": "stats.json",
                "alpha_dct": 0.02,
                "alpha_ks": 0.1,
                "discordant_test": "z-test",
                "summary_statistic": "median",
                "processes": 3,
                "no_multiprocessing": False,
            }
        )
    )

    config = load_triplet_processor_config(str(config_path))

    assert config["input"] == "unique_triplets_gene_trees.txt"
    assert config["output"] == "results.tsv"
    assert config["stats_output"] == "stats.json"
    assert config["alpha_dct"] == 0.02
    assert config["alpha_ks"] == 0.1
    assert config["discordant_test"] == "z-test"
    assert config["summary_statistic"] == "median"
    assert config["processes"] == 3
    assert config["no_multiprocessing"] is False


def test_load_triplet_processor_config_missing_input(tmp_path):
    config_path = tmp_path / "triplet_processor_bad.json"
    config_path.write_text(json.dumps({"alpha_dct": 0.02}))

    with pytest.raises(ConfigError, match="input_path"):
        load_triplet_processor_config(str(config_path))
