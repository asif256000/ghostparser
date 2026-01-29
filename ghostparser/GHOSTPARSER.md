# Ghostparser Module Documentation

Detailed documentation for the `ghostparser.tree_parser` module.

## Overview

The ghostparser module processes phylogenetic trees in Newick format, standardizes them by removing support values, filters low-quality trees, and extracts triplet subtrees from gene trees for downstream ghost introgression analysis.

## Module Functions

### Tree Reading and Validation

#### `read_tree_file(filepath)`

Read Newick trees from a file with validation.

**Arguments:**

- `filepath`: Path to the Newick tree file

**Returns:**

- List of Bio.Phylo tree objects

**Raises:**

- `FileNotFoundError`: If the file does not exist
- `ValueError`: If the file contains invalid Newick format

**Example:**

```python
from ghostparser.tree_parser import read_tree_file

trees = read_tree_file("species.tree")
print(f"Read {len(trees)} trees")
```

### Support Value Processing

#### `calculate_average_support(tree)`

Calculate the average support value for a tree from internal nodes.

**Arguments:**

- `tree`: A Bio.Phylo tree object

**Returns:**

- The average support value across all internal nodes, or `None` if no support values found

**Example:**

```python
avg_support = calculate_average_support(tree)
if avg_support is not None:
    print(f"Average support: {avg_support:.4f}")
```

#### `remove_support_values(tree)`

Remove support values (bootstrap/posterior probabilities) from internal nodes.

**Arguments:**

- `tree`: A Bio.Phylo tree object

**Returns:**

- The tree with support values removed from all internal nodes

**Example:**

```python
clean_tree = remove_support_values(tree)
```

### Tree Standardization

#### `standardize_tree(tree)`

Standardize a tree by removing support values while preserving branch lengths.

**Arguments:**

- `tree`: A Bio.Phylo tree object

**Returns:**

- The standardized tree with support values removed

**Example:**

```python
standardized = standardize_tree(tree)
```

### Formatting and Writing

#### `format_newick_with_precision(tree, decimal_places=10)`

Format a tree as Newick string with custom branch length precision.

**Arguments:**

- `tree`: A Bio.Phylo tree object
- `decimal_places`: Number of decimal places for branch lengths (default: 10)

**Returns:**

- Newick format string with specified precision and trailing zeros removed

**Example:**

```python
newick = format_newick_with_precision(tree, decimal_places=10)
print(newick)
```

#### `write_clean_trees(trees, output_filepath)`

Write cleaned trees to a file in Newick format.

**Arguments:**

- `trees`: List of Bio.Phylo tree objects
- `output_filepath`: Path where trees will be written

**Example:**

```python
write_clean_trees(cleaned_trees, "output_clean.tree")
```

#### `clean_and_save_trees(input_filepath, output_filepath, support_threshold=0.5)`

Clean trees and save to file, filtering by support threshold.

**Arguments:**

- `input_filepath`: Path to input Newick file
- `output_filepath`: Path for output file
- `support_threshold`: Minimum average support to keep trees (default: 0.5)

**Returns:**

- Tuple of (cleaned_trees, dropped_indices_dict)

**Example:**

```python
trees, dropped = clean_and_save_trees("input.tree", "output.tree")
print(f"Kept {len(trees)} trees, dropped {len(dropped)}")
```

### Taxa and Triplet Operations

#### `get_taxa_from_tree(tree)`

Extract all terminal taxa names from a phylogenetic tree.

**Arguments:**

- `tree`: A Bio.Phylo tree object

**Returns:**

- A sorted list of terminal taxa names

**Example:**

```python
taxa = get_taxa_from_tree(tree)
print(f"Found taxa: {', '.join(taxa)}")
```

#### `generate_triplets(taxa_list, outgroup)`

Generate all unique triplet combinations from taxa excluding the outgroup.

**Arguments:**

- `taxa_list`: List of all taxa names
- `outgroup`: The outgroup taxon to exclude

**Returns:**

- A list of tuples, where each tuple contains 3 taxa names

**Formula:** nC3 where n = len(taxa_list) - 1 (excluding outgroup)

**Example:**

```python
triplets = generate_triplets(taxa, "OutGroup")
print(f"Generated {len(triplets)} triplets")
```

#### `write_triplets_to_file(triplets, output_filepath)`

Write triplets to a file, one triplet per line.

**Arguments:**

- `triplets`: List of triplet tuples
- `output_filepath`: Path where the triplets will be written

**Example:**

```python
write_triplets_to_file(triplets, "triplets.txt")
```

### Triplet Subtree Extraction

#### `extract_triplet_subtree(tree, triplet_taxa)`

Extract a subtree containing only the specified triplet taxa.

**Arguments:**

- `tree`: A Bio.Phylo tree object
- `triplet_taxa`: List or tuple of 3 taxa names to extract

**Returns:**

- A new Bio.Phylo tree containing only the specified taxa, or `None` if not all taxa are present

**Details:**

- Creates a deep copy of the tree to avoid modifying the original
- Prunes away all taxa not in the triplet
- Recalculates branch lengths appropriately
- Branch length aggregation is handled by Bio.Phylo during pruning: when a parent ends up with a single child, Bio.Phylo collapses that edge and folds its length into the child, so root-to-leaf distances for the retained taxa match the original tree (we do not manually sum lengths)

**Example:**

```python
triplet = ("TaxaA", "TaxaB", "TaxaC")
subtree = extract_triplet_subtree(gene_tree, triplet)
if subtree:
    print(f"Extracted subtree with {len(subtree.get_terminals())} taxa")
```

#### `process_gene_trees_for_triplets(gene_trees, triplets)`

Process gene trees to extract subtrees for each triplet.

**Arguments:**

- `gene_trees`: List of Bio.Phylo tree objects
- `triplets`: List of triplet tuples

**Returns:**

- A dictionary mapping triplets to lists of subtree Newick strings

**Details:**

- For each gene tree, attempts to extract a subtree for each triplet
- Skips gene trees that don't contain all 3 taxa from a triplet
- Same gene tree can contribute to multiple triplets

**Example:**

```python
triplet_trees = process_gene_trees_for_triplets(gene_trees, triplets)
for triplet, trees in triplet_trees.items():
    print(f"{triplet}: {len(trees)} gene trees")
```

#### `write_triplet_gene_trees(triplet_gene_trees, output_filepath)`

Write triplet gene trees to a file in the specified format.

**Arguments:**

- `triplet_gene_trees`: Dictionary mapping triplets to lists of Newick strings
- `output_filepath`: Path where the results will be written

**Format:**

```
TaxonA,TaxonB,TaxonC<TAB>count

newick_tree_1
newick_tree_2
...


NextTriplet...
```

**Example:**

```python
write_triplet_gene_trees(triplet_trees, "triplet_gene_trees.txt")
```

### Utility Functions

#### `get_clean_filename(filepath)`

Generate a clean filename with \_clean prefix.

**Arguments:**

- `filepath`: The original file path

**Returns:**

- The path with \_clean inserted before the file extension

**Example:**

```python
clean_name = get_clean_filename("species.tree")
# Returns: "species_clean.tree"
```

## Command Line Interface

### Main Entry Point

#### `main()`

Main entry point for the tree parser module with argparse CLI.

**Usage:**

```bash
python -m ghostparser.tree_parser -st SPECIES_TREE -gt GENE_TREES -og OUTGROUP
```

**Arguments:**

**Required:**
- `-st, --species_tree`: Path to species tree file in Newick format
- `-gt, --gene_trees`: Path to gene trees file in Newick format
- `-og, --outgroup`: Outgroup species identifier

**Optional:**
- `-p, --processes`: Number of worker processes for multiprocessing (only used when multiprocessing is enabled).
                    Defaults to cpu_count(). Ignored if `--no-multiprocessing` is set.
- `--no-multiprocessing`: Disable multiprocessing for triplet extraction.
                         Processes triplets sequentially on a single worker.
                         Useful for debugging or systems with limited memory.

**Output Files:**

1. `{species_tree}_clean.tree` - Cleaned species tree
2. `{gene_trees}_clean.tree` - Cleaned gene trees
3. `{gene_trees}_unique_triplet_gene_trees.txt` - Triplet gene trees

**Example:**

```bash
python -m ghostparser.tree_parser \
    -st data/species.tree \
    -gt data/genes.tree \
    -og OutGroup
```

## Workflow Details

### Complete Processing Pipeline

1. **Initialization**

   - Parse command line arguments
   - Validate input files exist
   - Generate output filenames

2. **Species Tree Processing**

   - Read species tree
   - Calculate average support
   - Remove support values
   - Filter by support threshold
   - Write cleaned tree
   - Extract taxa names

3. **Triplet Generation**

   - Remove outgroup from taxa list
   - Generate all nC3 combinations
   - Display triplet count

4. **Gene Tree Processing**

   - Read gene trees
   - Calculate average support for each
   - Remove support values
   - Filter by support threshold
   - Write cleaned trees

5. **Triplet Extraction**

   - For each gene tree and triplet combination:
     - Check if gene tree contains all 3 taxa
     - Extract subtree if all present
     - Store subtree under triplet key
   - Write results to triplet gene trees file

6. **Statistics and Reporting**
   - Total triplets generated
   - Triplets with matching gene trees
   - Total subtrees extracted
   - Average trees per triplet

## Configuration

### Support Threshold

Default: `0.5`

Trees with average support below this threshold are filtered out. This ensures only high-quality trees are used for analysis.

### Branch Length Precision

Default: `10 decimal places`

Branch lengths are preserved with high precision and trailing zeros are removed for cleaner output.

## Dependencies

- **BioPython (>=1.79)**: For Phylo module (tree parsing and manipulation)
- **Python 3.x**: Uses f-strings, pathlib, type hints

## Error Handling

The module provides robust error handling:

- **Invalid Newick Format**: Raises `ValueError` with descriptive message
- **Missing Files**: Raises `FileNotFoundError` with file path
- **Empty Trees**: Validated during reading
- **Missing Taxa**: Returns `None` for invalid triplet extractions
- **Missing Outgroup**: Warns user and shows available taxa
