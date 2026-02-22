"""Configuration loading helpers for GhostParser CLIs."""

from __future__ import annotations

import json
from pathlib import Path


class ConfigError(ValueError):
    """Raised when a config file is invalid or missing required fields."""


def _load_raw_config(config_file: str) -> dict:
    """Load a raw config dictionary from JSON or YAML."""
    path = Path(config_file)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {config_file}")

    suffix = path.suffix.lower()
    if suffix == ".json":
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    elif suffix in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as exc:
            raise ConfigError("YAML support requires PyYAML to be installed") from exc

        with open(path, "r", encoding="utf-8") as handle:
            payload = yaml.safe_load(handle)
    else:
        raise ConfigError("Config file must be .json, .yaml, or .yml")

    if not isinstance(payload, dict):
        raise ConfigError("Config root must be a key/value object")

    return payload


def _validate_required_string(payload: dict, key: str) -> str:
    """Validate a required non-empty string field from config payload."""
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Missing required config field: {key}")
    return value.strip()


def _parse_outgroups(value) -> list[str]:
    """Parse outgroup(s) value from config.

    Accepts either:
    - a comma-separated string
    - a list/tuple/set of strings
    """
    if isinstance(value, str):
        outgroups = [part.strip() for part in value.split(",") if part.strip()]
    elif isinstance(value, (list, tuple, set)):
        outgroups = [str(part).strip() for part in value if str(part).strip()]
    else:
        outgroups = []

    if not outgroups:
        raise ConfigError("Missing required config field: outgroup(s)")
    return outgroups


def load_orchestrator_config(config_file: str) -> dict:
    """Load and normalize orchestrator config values.

    Required keys:
    - species_tree_path
    - gene_trees_path
    - outgroup or outgroups

    Supported optional keys:
    - output_folder
    - processes
    - triplet_filter
    - min_support_value
    """
    payload = _load_raw_config(config_file)

    species_tree = _validate_required_string(payload, "species_tree_path")
    gene_trees = _validate_required_string(payload, "gene_trees_path")

    outgroups_source = payload.get("outgroups")
    if outgroups_source is None:
        outgroups_source = payload.get("outgroup")
    outgroups = _parse_outgroups(outgroups_source)

    output = payload.get("output_folder")
    if output is not None:
        if not isinstance(output, str) or not output.strip():
            raise ConfigError("Config field output_folder must be a non-empty string when provided")
        output = output.strip()

    triplet_filter = payload.get("triplet_filter")
    if triplet_filter is not None:
        if not isinstance(triplet_filter, str) or not triplet_filter.strip():
            raise ConfigError("Config field triplet_filter must be a non-empty string when provided")
        triplet_filter = triplet_filter.strip()

    processes = payload.get("processes")
    if processes is not None:
        if not isinstance(processes, int) or processes < 0:
            raise ConfigError("Config field processes must be an integer >= 0")

    min_support_value = payload.get("min_support_value")
    if min_support_value is not None:
        try:
            min_support_value = float(min_support_value)
        except (TypeError, ValueError) as exc:
            raise ConfigError("Config field min_support_value must be a numeric value") from exc

    return {
        "species_tree": species_tree,
        "gene_trees": gene_trees,
        "outgroup": outgroups,
        "triplet_filter": triplet_filter,
        "output": output,
        "processes": processes,
        "min_support_value": min_support_value,
    }
