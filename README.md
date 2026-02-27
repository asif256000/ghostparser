# Ghostparser

A phylogenetic introgression pipeline centered on `ghostparser.orchestrator`. The orchestrator runs two internal stages (`tree_parser` and `triplet_processor`) and produces final triplet-level inference outputs.

## Features

- Orchestrator-first end-to-end run (`tree_parser` + `triplet_processor`)
- JSON/YAML config support with shared CLI/config normalization
- Deterministic triplet orientation (`A,B,C`) and rooted species-triplet headers
- Configurable DCT/KS thresholds and statistical backend
- Optional multiprocessing with `processes=0` meaning all cores

## Quick Start

### Installation

```bash
pip install -r requirements.txt
```

### Installation (Poetry 2.x+)

If you prefer Poetry 2.x+ instead of plain pip:

```bash
poetry install
```

Run commands with Poetry:

```bash
poetry run pytest -q
poetry run python -m ghostparser.orchestrator -c run_config.yaml
```

## Orchestrator (Primary Pipeline)

Run the full pipeline (recommended):

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup
```

Run with config file:

```bash
python -m ghostparser.orchestrator -c run_config.yaml
```

CLI mode with explicit worker count:

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup --processes 0
```

### Why this pipeline works

- `tree_parser` standardizes species/gene trees, roots on outgroup(s), and writes triplet-specific gene-tree blocks.
- `triplet_processor` applies the GhostParser statistical decision pipeline to each triplet block.
- The orchestrator coordinates both steps and writes the final results table.

GhostParser is configurable (discordant test, backend, thresholds, summary statistic), so execution follows the same core pipeline stages while allowing controlled method choices.

### Orchestrator arguments

- **Required**
    - `-st, --species-tree-path`
    - `-gt, --gene-trees-path`
    - `-og, --outgroups`
- **Config mode**
    - `-c, --config-file`
- **Common optional**
    - `--output-folder`
    - `--triplet-filter`
    - `--processes`
    - `--min-support-value`
    - `--discordant-test`
    - `--summary-statistic`
    - `--stats-backend`
    - `--alpha-dct`, `--alpha-ks`

### Primary outputs

1. `unique_triplets_gene_trees.txt`
2. `orchestrator_triplet_results.tsv`
3. `metrics.txt`

## Tree Parser (Submodule)

Use this module when you only want preprocessing + triplet extraction.

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup
```

Config mode:

```bash
python -m ghostparser.tree_parser -c sample_configs/tree_parser_minimal.yaml
```

Core behavior:

- Removes support labels and preserves branch lengths
- Roots on outgroup(s), prunes outgroup clade, and logs excluded taxa
- Writes processed trees and `unique_triplets_gene_trees.txt`

## Triplet Processor (Submodule)

Use this module when you already have `unique_triplets_gene_trees.txt` and only need inference.

```bash
python -m ghostparser.triplet_processor --input-path unique_triplets_gene_trees.txt
```

Config mode:

```bash
python -m ghostparser.triplet_processor -c sample_configs/triplet_processor_minimal.yaml
```

Core behavior:

- Runs DCT (`chi-square` or `z-test`)
- Runs KS tree-height test when DCT is significant
- Applies summary-statistic comparison (`median` or `mean`) for final classification

## Example Input/Output (Tree Parser)

The arguments below are for `tree_parser`:

**Required:**
- `--species-tree-path`: Path to the species tree file in Newick format
- `--gene-trees-path`: Path to the gene trees file in Newick format
- `--outgroups`: Outgroup species identifier(s). Use comma-separated taxa for multiple outgroups.

Short aliases:
- `-st` for `--species-tree-path`
- `-gt` for `--gene-trees-path`
- `-og` for `--outgroups`
- `-c` for `--config-file`

**Optional:**
- `--output-folder`: Output folder relative to the input data folder (default: same folder as input data).
- `--triplet-filter`: Path to a triplet filter file (comma-separated taxa per line). When provided, only those
    triplets are processed. Triplets containing taxa missing from the species tree are skipped with a warning.
- `--processes`: Number of worker processes for multiprocessing.
                    Defaults to `0` (all cores). Ignored if `--no-multiprocessing` is set.
- `--no-multiprocessing`: Disable multiprocessing. Processes triplets sequentially using a single worker.
                         Useful for debugging or on systems with limited resources.

## Output Files

The tool generates these output files:

1. **`processed_<species_tree>`** - Processed species tree with support values removed and outgroup rooting applied
2. **`processed_<gene_trees>`** - Processed gene trees with support values removed and outgroup rooting applied
3. **`unique_triplets_gene_trees.txt`** - Triplets normalized to `A,B,C` (with `A,B` as species sisters), with required header format `triplet<TAB>count<TAB>species_tree` (non-empty species subtree)
4. **`metrics.txt`** - Metrics log with warnings, timings, and counts
5. **`triplet_introgression_results.tsv`** - Final triplet-level classification results (`no_introgression`, `outflow_introgression`, `inflow_introgression`, `ghost_introgression`, or `unresolved`)

## Example

**Species tree** (`species.tree`):

```
(((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4):0.5,(TaxaD:0.6,OutGroup:0.7):0.8);
```

**Gene trees** (`genes.tree`):

```
((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
```

**Run:**

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup
```

Multiple outgroups (comma-separated):

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup1,OutGroup2
```

Triplet filter file (comma-separated taxa per line):

`triplets.txt` is an input filter file; only those triplets are processed.

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup --triplet-filter triplets.txt
```

When multiple outgroups are provided, the species tree is rooted on their most recent common ancestor (MRCA) and
the outgroup clade is pruned. Any additional taxa that fall inside the outgroup clade are excluded from triplet
generation and logged as a warning (including the full list of excluded taxa) in the metrics file.

**Output** (`unique_triplets_gene_trees.txt`):

```
TaxaA,TaxaB,TaxaC	2	((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);

((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.13,TaxaB:0.24):0.35,TaxaC:0.46):0.57;


TaxaA,TaxaC,TaxaD	1	((TaxaA:0.1,TaxaC:0.2):0.3,TaxaD:0.4);

((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
```

## Documentation

- **[Orchestrator Module Guide](ghostparser/ORCHESTRATOR.md)** - End-to-end orchestrator inputs, outputs, and result columns
- **[Configuration Guide](CONFIG.md)** - Orchestrator-first config keys and module-specific config sections
- **[Test Documentation](tests/TESTS.md)** - Orchestrator-first testing map and full test index
- **[Module Documentation](ghostparser/GHOSTPARSER.md)** - Full API details for all modules

Note: module-specific guides live under the [ghostparser](ghostparser/) folder.

## Defaults at a Glance

Core defaults are centralized and applied consistently in both CLI mode and config-file mode:

- `discordant_test`: `chi-square`
- `summary_statistic`: `median`
- `stats_backend`: `standard`
- `alpha_dct`: `0.01`
- `alpha_ks`: `0.05`
- `processes`: `0` (all available CPU cores)
- `output_folder` (orchestrator/tree_parser): `./results`
- `min_support_value`: `0.5`

Backend details:

- `custom`: uses GhostParser manual implementations for chi-square, two-proportion z-test, and KS test
- `standard`: uses SciPy for chi-square/KS and statsmodels for two-proportion z-test

When `-c/--config-file` is provided, other CLI flags are ignored with a warning.

Primary orchestrator outputs:

- `unique_triplets_gene_trees.txt` (triplet-wise extracted gene trees used as input to stage 2)
- `orchestrator_triplet_results.tsv` (final per-triplet inference/statistics table)
- `metrics.txt` (run log, counts, timings, warnings)

## Testing

Run all tests:

```bash
pytest
```

See [tests/TESTS.md](tests/TESTS.md) for detailed test documentation.