# Ghostparser

A summary-statistics based phylogenetic tool for identifying introgressions in reasonably large phylogenetic dataset. This tool processes Newick format trees (species and genes), standardizes them, and extracts triplet subtrees from gene trees for downstream analysis. In the end, it runs a triplet-based pipeline to classify introgression types for each triplet of taxa.

## Features

- **Tree Standardization**: Removes support values (bootstrap/posterior probabilities) from phylogenetic trees.
- **Quality Filtering**: Filters out trees with average support less than configured value.
- **High-Precision Branch Lengths**: Preserves branch lengths with 10 decimal places, removes trailing zeros.
- **Triplet Generation**: Generates all unique triplet combinations from taxa (nC3).
- **Triplet Subtree Extraction**: Extracts subtrees from gene trees for each triplet with proper branch length calculations.
- **Low-Memory Parallel Triplets**: Processes triplet chunks in parallel while streaming gene trees from disk to avoid large in-memory triplet maps.
- **Triplet Classification Pipeline**: Classifies introgression types for each triplet based on gene tree frequencies and tree height distributions.
- **Config File Support**: Allows running the full pipeline from a JSON/YAML config file with all parameters specified.

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

### Basic Usage

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup
```

Or run `tree_parser` from a config file:

```bash
python -m ghostparser.tree_parser -c sample_configs/tree_parser_minimal.yaml
```

### Arguments

The arguments below are for `tree_parser`:

**Required:**
- `-st, --species_tree`: Path to the species tree file in Newick format
- `-gt, --gene_trees`: Path to the gene trees file in Newick format
- `-og, --outgroup`: Outgroup species identifier(s). Use comma-separated taxa for multiple outgroups.

**Optional:**
- `-o, --output`: Output folder relative to the input data folder (default: same folder as input data).
- `-tf, --triplet-filter`: Path to a triplet filter file (comma-separated taxa per line). When provided, only those
    triplets are processed. Triplets containing taxa missing from the species tree are skipped with a warning.
- `-p, --processes`: Number of worker processes for multiprocessing (only used with `-p/--processes`).
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
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup -tf triplets.txt
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

- **[Module Documentation](ghostparser/GHOSTPARSER.md)** - Detailed API documentation for all functions
- **[Orchestrator Module Guide](ghostparser/ORCHESTRATOR.md)** - End-to-end orchestrator inputs, outputs, and result columns
- **[Test Documentation](tests/TESTS.md)** - Comprehensive test suite details and categories
- **[Configuration Guide](CONFIG.md)** - Config file options and examples for orchestrator

Note: module-specific guides live under the [ghostparser](ghostparser/) folder.

## Triplet Processing (Fig. 6)

After generating `unique_triplets_gene_trees.txt` with `tree_parser`, run:

```bash
python -m ghostparser.triplet_processor -i unique_triplets_gene_trees.txt
```

Optional thresholds and output path:

```bash
python -m ghostparser.triplet_processor \
    -i unique_triplets_gene_trees.txt \
    -o triplet_introgression_results.tsv \
    --discordant-test z-test \
    --summary-statistic median \
    --stats-backend custom \
    --alpha-dct 0.01 \
    --alpha-ks 0.05
```

Or run `triplet_processor` from a config file:

```bash
python -m ghostparser.triplet_processor -c sample_configs/triplet_processor_minimal.yaml
```

Optional multiprocessing for `triplet_processor`:

```bash
python -m ghostparser.triplet_processor \
    -i unique_triplets_gene_trees.txt \
    -p 0
```

- `-p` defaults to `0`, which uses all available CPU cores
- `--no-multiprocessing` forces single-worker analysis
- `--discordant-test` supports `z-test` (default) or `chi-square`
- `--summary-statistic` supports `median` (default) or `mean`
- `--stats-backend` supports `custom` (default) or `standard`
- `--alpha-dct` and `--alpha-ks` are configurable p-value thresholds (defaults: `0.01` and `0.05`)
- when `-c/--config-file` is provided, other CLI options are ignored with a warning

CLI mode and config-file mode share the same normalization and default-validation path.

Triplet topology convention in `triplet_processor`:

- triplets are treated as `A,B,C` where `A,B` are species sisters
- `con` is the species-matching topology (`((A,B),C)`), `dis1` is the more frequent discordant topology, and `dis2` is the less frequent discordant topology (on ties, first discordant topology is used as `dis1`)
- output includes `most_frequent_matches_concordant` (`True` when concordant frequency is not lower than either discordant frequency)
- output includes `species_tree` (the extracted species-tree Newick for that triplet)

## End-to-End Orchestrator

Run the full pipeline (tree parsing + per-triplet inference) in one command:

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup
```

Or run from a JSON/YAML config file:

```bash
python -m ghostparser.orchestrator -c run_config.yaml
```

CLI mode (without config file):

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup -p 0
```

Notes:

- `-p` defaults to `0`, which uses all available CPU cores
- `-o/--output` defaults to `./results` (current working directory)
- `-c/--config-file` can be combined with other CLI options, but CLI values are ignored with a warning
- use `-p 1` for effective single-worker execution
- orchestrator runs in two stages: first generate `unique_triplets_gene_trees.txt`, then run `triplet_processor` on that file

For all config keys and examples, see [CONFIG.md](CONFIG.md).

## Defaults at a Glance

Core defaults are centralized and applied consistently in both CLI mode and config-file mode:

- `discordant_test`: `z-test`
- `summary_statistic`: `median`
- `stats_backend`: `custom`
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