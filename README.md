# Ghostparser

A phylogenetic tree parser for identifying ghost introgressions. This tool processes Newick format trees, standardizes them, and extracts triplet subtrees from gene trees for downstream analysis.

## Features

- **Tree Standardization**: Removes support values (bootstrap/posterior probabilities) from phylogenetic trees
- **Quality Filtering**: Filters out trees with average support < 0.5
- **High-Precision Branch Lengths**: Preserves branch lengths with 10 decimal places, removes trailing zeros
- **Triplet Generation**: Generates all unique triplet combinations from taxa (nC3)
- **Triplet Subtree Extraction**: Extracts subtrees from gene trees for each triplet with proper branch length calculations
- **Low-Memory Parallel Triplets**: Processes triplet chunks in parallel while streaming gene trees from disk to avoid large in-memory triplet maps
- **Triplet Classification Pipeline**: Runs GhostParser Fig. 6 DCT + THT + median decision logic for each triplet

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

1. **`processed_<species_tree>`** - Processed species tree with support values removed and outgroup rooting applied
2. **`processed_<gene_trees>`** - Processed gene trees with support values removed and outgroup rooting applied
3. **`unique_triplets_gene_trees.txt`** - Triplets normalized to `A,B,C` (with `A,B` as species sisters), with gene-tree count and species-triplet subtree in header
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
- **[Test Documentation](tests/TESTS.md)** - Comprehensive test suite details and categories

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
    --alpha-dct 0.01 \
    --alpha-ks 0.05
```

Triplet topology convention in `triplet_processor`:

- triplets are treated as `A,B,C` where `A,B` are species sisters
- inference uses species-anchored roles: `con` always matches the species topology, and `dis1`/`dis2` are the two discordant topologies
- for reporting, topologies are also ranked by frequency (`top1/top2/top3`), and any highest-frequency topology is tagged as `[highest freq]`
- if species-concordant topology is not highest-frequency, `con` is additionally tagged as `[diff]`

## End-to-End Orchestrator

Run the full pipeline (tree parsing + per-triplet inference) in one command:

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup
```

Optional parallelism control:

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup -p 0
```

Notes:

- `-p 0` uses all available CPU cores
- `--no-multiprocessing` forces single-worker execution
- `--log-triplet-gene-trees` enables debug logging of generated triplets and extracted gene trees to `unique_triplets_gene_trees.txt`

Primary orchestrator outputs:

- `orchestrator_triplet_results.tsv` (final per-triplet inference/statistics table)
- `metrics.txt` (run log, counts, timings, warnings)

Optional debug output:

- `unique_triplets_gene_trees.txt` (only when `--log-triplet-gene-trees` is used)

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
