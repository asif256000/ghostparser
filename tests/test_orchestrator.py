"""Tests for orchestrator module."""

import argparse

from ghostparser.orchestrator import (
    _resolve_parallel_mode,
    _resolve_runtime_args,
    _resolve_processes,
)


def test_resolve_processes_zero_uses_all_cores(monkeypatch):
    monkeypatch.setattr("ghostparser.orchestrator.cpu_count", lambda: 16)
    assert _resolve_processes(0) == 16
    assert _resolve_processes(None) is None
    assert _resolve_processes(4) == 4


def test_resolve_parallel_mode(monkeypatch):
    monkeypatch.setattr("ghostparser.orchestrator.cpu_count", lambda: 8)

    processes, use_multiprocessing = _resolve_parallel_mode(0)
    assert processes == 8
    assert use_multiprocessing is True

    processes, use_multiprocessing = _resolve_parallel_mode(1)
    assert processes == 1
    assert use_multiprocessing is False

    processes, use_multiprocessing = _resolve_parallel_mode(4)
    assert processes == 4
    assert use_multiprocessing is True


def test_resolve_runtime_args_cli_defaults(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = argparse.Namespace(
        config_file=None,
        species_tree="species.nwk",
        gene_trees="genes.nwk",
        outgroup="Out1,Out2",
        triplet_filter=None,
        output=str(tmp_path / "results"),
        processes=0,
        min_support_value=None,
        discordant_test=None,
        summary_statistic=None,
        stats_backend=None,
        alpha_dct=None,
        alpha_ks=None,
    )

    resolved = _resolve_runtime_args(args)

    assert resolved.species_tree == "species.nwk"
    assert resolved.gene_trees == "genes.nwk"
    assert resolved.outgroup == ["Out1", "Out2"]
    assert resolved.output == str(tmp_path / "results")
    assert resolved.processes == 0
    assert resolved.min_support_value == 0.5
    assert resolved.discordant_test == "z-test"
    assert resolved.summary_statistic == "median"
    assert resolved.stats_backend == "custom"
    assert resolved.alpha_dct == 0.01
    assert resolved.alpha_ks == 0.05


def test_resolve_runtime_args_cli_custom_processes_preserved(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = argparse.Namespace(
        config_file=None,
        species_tree="species.nwk",
        gene_trees="genes.nwk",
        outgroup="Out1,Out2",
        triplet_filter=None,
        output=str(tmp_path / "results"),
        processes=5,
        min_support_value=None,
        discordant_test=None,
        summary_statistic=None,
        stats_backend=None,
        alpha_dct=None,
        alpha_ks=None,
    )

    resolved = _resolve_runtime_args(args)
    assert resolved.processes == 5


def test_resolve_runtime_args_config_with_cli_warns_and_ignores(tmp_path, capsys):
    config_path = tmp_path / "run_config.json"
    config_path.write_text(
        '{"species_tree_path":"s.nwk","gene_trees_path":"g.nwk","outgroup":"OutA"}'
    )

    args = argparse.Namespace(
        config_file=str(config_path),
        species_tree="species.nwk",
        gene_trees=None,
        outgroup=None,
        triplet_filter=None,
        output=str(tmp_path / "results"),
        processes=7,
        min_support_value=0.9,
        discordant_test="chi-square",
        summary_statistic="mean",
        stats_backend="standard",
        alpha_dct=0.2,
        alpha_ks=0.3,
    )

    resolved = _resolve_runtime_args(args)
    captured = capsys.readouterr()

    assert "Warning: --config-file provided; CLI arguments not in config will be ignored" in captured.out
    assert resolved.species_tree == "s.nwk"
    assert resolved.gene_trees == "g.nwk"
    assert resolved.outgroup == ["OutA"]
    assert resolved.processes == 0
    assert resolved.discordant_test == "z-test"
    assert resolved.summary_statistic == "median"
    assert resolved.stats_backend == "custom"
    assert resolved.alpha_dct == 0.01
    assert resolved.alpha_ks == 0.05


def test_resolve_runtime_args_config_processes_preserved_when_set(tmp_path):
    config_path = tmp_path / "run_config_with_processes.json"
    config_path.write_text(
        '{"species_tree_path":"s.nwk","gene_trees_path":"g.nwk","outgroup":"OutA","processes":4}'
    )

    args = argparse.Namespace(
        config_file=str(config_path),
        species_tree=None,
        gene_trees=None,
        outgroup=None,
        triplet_filter=None,
        output=str(tmp_path / "results"),
        processes=0,
        min_support_value=None,
        discordant_test=None,
        summary_statistic=None,
        stats_backend=None,
        alpha_dct=None,
        alpha_ks=None,
    )

    resolved = _resolve_runtime_args(args)
    assert resolved.processes == 4
