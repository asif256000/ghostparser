# Configuration Guide

This document describes configuration file support for:

- `ghostparser.orchestrator`
- `ghostparser.tree_parser`
- `ghostparser.triplet_processor`

## CLI Usage

Run orchestrator with a config file:

```bash
python -m ghostparser.orchestrator -c sample_configs/orchestrator_minimal.yaml
```

Run orchestrator with plain CLI arguments:

```bash
python -m ghostparser.orchestrator \
  -st data/asuv21/species.tree \
  -gt data/asuv21/gene_trees.tree \
  -og Ephemera_danica,Isonychia_kiangsinensis
```

If both `-c/--config-file` and other CLI args are provided, config mode is used and the other CLI args are ignored with a warning.

`tree_parser` config-file mode:

```bash
python -m ghostparser.tree_parser -c sample_configs/tree_parser_minimal.yaml
```

`triplet_processor` config-file mode:

```bash
python -m ghostparser.triplet_processor -c sample_configs/triplet_processor_minimal.yaml
```

---

## Supported Config Formats

- `.json`
- `.yaml`
- `.yml`

---

## Orchestrator Config Keys

### Required

- `species_tree_path` (string)
  - Path to species tree file.
- `gene_trees_path` (string)
  - Path to gene trees file.
- `outgroups` (list of strings)
  - Example:
    ```yaml
    outgroups:
      - Taxon1
      - Taxon2
    ```

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
- `discordant_test` (string)
  - Discordant count test method used by `triplet_processor` stage.
  - Allowed values: `chi-square` (default), `z-test`.
- `summary_statistic` (string)
  - Statistic used for con/dis1 distributions after KS test.
  - Allowed values: `mean` (default), `median`.
- `alpha_dct` (number)
  - P-value threshold for the discordant count test.
  - Default: `0.01`.
- `alpha_ks` (number)
  - P-value threshold for KS tree-height test.
  - Default: `0.05`.

---

## Sample Configs

See examples in:

- `sample_configs/orchestrator_minimal.yaml`
- `sample_configs/orchestrator_full.yaml`
- `sample_configs/orchestrator_full.json`

---

## Tree Parser Config Keys

### Required

- `species_tree_path` (string)
- `gene_trees_path` (string)
- `outgroups` (list of strings)
  - Example:
    ```yaml
    outgroups:
      - Taxon1
      - Taxon2
    ```

### Optional

- `output_folder` (string)
  - Output folder relative to the input species-tree folder.
- `processes` (integer >= 0)
  - Worker count for triplet extraction (`0` = all cores).
- `triplet_filter` (string)
  - Path to optional triplet filter file.
- `min_support_value` (number)
  - Support filtering threshold for species and gene tree cleaning.
  - Default: `0.5`.
- `no_multiprocessing` (boolean)
  - `true` forces single-worker extraction.

### Tree Parser Sample Configs

- `sample_configs/tree_parser_minimal.yaml`
- `sample_configs/tree_parser_full.yaml`

---

## Triplet Processor Config Keys

### Required

- `input_path` (string)
  - Path to `unique_triplets_gene_trees.txt`.

### Optional

- `output_path` (string)
  - Output TSV path.
  - Default: `<input_dir>/triplet_introgression_results.tsv`.
- `stats_output` (string)
  - Optional JSON statistics output path.
  - Default: same path as output TSV with `.json` extension.
- `alpha_dct` (number)
  - Default: `0.01`.
- `alpha_ks` (number)
  - Default: `0.05`.
- `discordant_test` (string)
  - Allowed values: `chi-square` (default), `z-test`.
- `summary_statistic` (string)
  - Allowed values: `mean` (default), `median`.
- `processes` (integer >= 0)
  - Worker count for triplet inference (`0` = all cores).
- `no_multiprocessing` (boolean)
  - `true` forces single-worker analysis.

### Triplet Processor Sample Configs

- `sample_configs/triplet_processor_minimal.yaml`
- `sample_configs/triplet_processor_full.yaml`

---

## Notes

- Paths in config are interpreted relative to the current working directory unless absolute paths are used.
- In config mode, CLI options other than `-c/--config-file` are ignored.
- Use `-p 1` (or `processes: 1`) to avoid multiprocessing while still using the same pipeline.
