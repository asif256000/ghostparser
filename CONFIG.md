# Configuration Guide

This document describes configuration file support for `ghostparser.orchestrator`.

## CLI Usage

Run with a config file:

```bash
python -m ghostparser.orchestrator -c sample_configs/orchestrator_minimal.yaml
```

Run with plain CLI arguments:

```bash
python -m ghostparser.orchestrator \
  -st data/asuv21/species.tree \
  -gt data/asuv21/gene_trees.tree \
  -og Ephemera_danica,Isonychia_kiangsinensis
```

If both `-c/--config-file` and other CLI args are provided, config mode is used and the other CLI args are ignored with a warning.

---

## Supported Config Formats

- `.json`
- `.yaml`
- `.yml`

---

## Config Keys

### Required

- `species_tree_path` (string)
  - Path to species tree file.
- `gene_trees_path` (string)
  - Path to gene trees file.
- `outgroup` or `outgroups`
  - `outgroup`: comma-separated string (example: `"Taxon1,Taxon2"`)
  - `outgroups`: list of strings (example: `["Taxon1", "Taxon2"]`)

### Optional

- `output_folder` (string)
  - Output directory path.
  - Default: `./results` from the current working directory.
- `processes` (integer >= 0)
  - Worker count for extraction and inference.
  - `0` means all available CPU cores.
  - `1` means single-worker execution (no multiprocessing).
- `triplet_filter` (string)
  - Path to triplet filter file (comma-separated taxa per line).
- `min_support_value` (number)
  - Support filtering threshold for species and gene tree cleaning.
  - Default behavior when omitted is equivalent to `0.5`.

---

## Sample Configs

See examples in:

- `sample_configs/orchestrator_minimal.yaml`
- `sample_configs/orchestrator_full.yaml`
- `sample_configs/orchestrator_full.json`

---

## Notes

- Paths in config are interpreted relative to the current working directory unless absolute paths are used.
- In config mode, CLI options other than `-c/--config-file` are ignored.
- Use `-p 1` (or `processes: 1`) to avoid multiprocessing while still using the same pipeline.
