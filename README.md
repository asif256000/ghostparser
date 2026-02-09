# Ghostparser

A phylogenetic tree parser for identifying ghost introgressions. This tool processes Newick format trees, standardizes them, and extracts triplet subtrees from gene trees for downstream analysis.

## Features

- **Tree Standardization**: Removes support values (bootstrap/posterior probabilities) from phylogenetic trees
- **Quality Filtering**: Filters out trees with average support < 0.5
- **High-Precision Branch Lengths**: Preserves branch lengths with 10 decimal places, removes trailing zeros
- **Triplet Generation**: Generates all unique triplet combinations from taxa (nC3)
- **Triplet Subtree Extraction**: Extracts subtrees from gene trees for each triplet with proper branch length calculations
- **Low-Memory Parallel Triplets**: Processes triplet chunks in parallel while streaming gene trees from disk to avoid large in-memory triplet maps

## Quick Start

### Installation

```bash
pip install -r requirements.txt
```

### Basic Usage

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup
```

### Arguments

**Required:**
- `-st, --species_tree`: Path to the species tree file in Newick format
- `-gt, --gene_trees`: Path to the gene trees file in Newick format
- `-og, --outgroup`: Outgroup species identifier(s). Use comma-separated taxa for multiple outgroups.

**Optional:**
- `-o, --output`: Output folder relative to the input data folder (default: same folder as input data).
- `-tf, --triplet-filter`: Path to a triplet filter file (comma-separated taxa per line). When provided, only those
    triplets are processed. Triplets containing taxa missing from the species tree are skipped with a warning.
- `-p, --processes`: Number of worker processes for multiprocessing (only used with `-p/--processes`).
                    Defaults to cpu_count(). Ignored if `--no-multiprocessing` is set.
- `--no-multiprocessing`: Disable multiprocessing. Processes triplets sequentially using a single worker.
                         Useful for debugging or on systems with limited resources.

## Output Files

The tool generates these output files:

1. **`*_clean.tree`** - Cleaned trees with support values removed and low-quality trees filtered (all downstream processing uses these cleaned files)
2. **`*_unique_triplet_gene_trees.txt`** - Triplet combinations with their corresponding gene trees
3. **`*_metrics.txt`** - Metrics log with warnings, timings, and counts

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

**Output** (`genes_unique_triplet_gene_trees.txt`):

```
TaxaA,TaxaB,TaxaC	2

((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.13,TaxaB:0.24):0.35,TaxaC:0.46):0.57;


TaxaA,TaxaC,TaxaD	1

((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
```

## Documentation

- **[Module Documentation](ghostparser/GHOSTPARSER.md)** - Detailed API documentation for all functions
- **[Test Documentation](tests/TESTS.md)** - Comprehensive test suite details and categories

## Testing

Run all tests:

```bash
pytest
```

Run with verbose output:

```bash
pytest -v
```

See [tests/TESTS.md](tests/TESTS.md) for detailed test documentation.

## Development

Set up a development environment:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```
