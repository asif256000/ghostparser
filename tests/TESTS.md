# Test Suite Documentation

The ghostparser project includes comprehensive test coverage with 42 test cases organized into multiple categories.

## Running Tests

### Run All Tests

```bash
pytest
```

### Run Tests with Verbose Output

```bash
pytest -v
```

### Run Specific Test File

```bash
pytest tests/test_tree_parser.py -v
```

### Run Specific Test

```bash
pytest tests/test_tree_parser.py::test_extract_triplet_subtree_all_taxa_present -v
```

### Run Tests with Coverage (requires `pytest-cov`)

Install the plugin if you don't have it:

```bash
pip install pytest-cov
```

Then run:

```bash
pytest --cov=ghostparser --cov-report=html
```
# Test Catalog

This file summarizes inputs, outputs, and intent for each test function. Shared fixtures live in tests/fixtures.py.

## Shared Fixtures (tests/fixtures.py)

- simple_newick_file
	- Input file name: simple_tree.nwk
	- File content (Newick):
		- (TaxaA:0.001,(TaxaB:0.098,(((TaxaC:0.001,TaxaD:0.001):0.001,TaxaE:0.001):0.086,(TaxaF:0.001,TaxaG:0.001):0.032):0.001):0.012,OutGroup:0.558);
	- Used in tests (expected outputs):
		- test_read_tree_file_single_tree: len(trees) == 1 and trees[0] is Tree
		- test_calculate_average_support_no_values: avg_support is None
		- test_standardize_tree_preserves_branch_lengths: branch lengths unchanged within 1e-10
		- test_format_newick_with_precision_trailing_zeros: no "0000000000" and endswith ";"
		- test_format_newick_with_precision_default_places: contains "(" and ")" and endswith ";"
		- test_format_newick_with_custom_precision: endswith ";" and decimal parts length <= 5
		- test_write_clean_trees: output exists and len(output_trees) == 1
		- test_clean_and_save_trees_no_filters: cleaned == 1 and dropped == 0
		- test_clean_and_save_trees_creates_output_file: output exists and size > 0
		- test_get_taxa_from_tree: len(taxa) == 8, contains TaxaA and OutGroup, sorted
		- test_get_taxa_from_tree_correct_names: taxa == ["OutGroup", "TaxaA", "TaxaB", "TaxaC", "TaxaD", "TaxaE", "TaxaF", "TaxaG"]
		- test_integration_full_workflow: len(output_trees) == 1
		- test_integration_triplets_workflow: len(taxa) == 8, len(triplets) == 35, len(lines) == 35

- newick_with_support_file
	- Input file name: tree_with_support.nwk
	- File content (Newick with support):
		- (((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup)0.85;
	- Used in tests (expected outputs):
		- test_calculate_average_support_with_values: 0.93 < avg_support < 0.95
		- test_remove_support_values: avg_support becomes None after removal
		- test_standardize_tree_removes_support: avg_support is None

- multiple_trees_file
	- Input file name: multiple_trees.nwk
	- File content (3 lines):
		- (TaxaA:0.001,(TaxaB:0.098,(TaxaC:0.001,TaxaD:0.001):0.001):0.012,OutGroup:0.558);
		- (TaxaB:0.098,(TaxaC:0.001,TaxaD:0.001):0.001,OutGroup:0.558);
		- ((TaxaC:0.001,TaxaD:0.001):0.001,(TaxaB:0.098,TaxaA:0.001):0.012,OutGroup:0.558);
	- Used in tests (expected outputs):
		- test_read_tree_file_multiple_trees: len(trees) == 3 and all are Tree
		- test_write_clean_trees_multiple: len(output_trees) == 3

- low_support_tree_file
	- Input file name: low_support_trees.nwk
	- File content (2 lines):
		- (((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup);
		- (((TaxaC,TaxaD)0.3:0.110599,(TaxaF,TaxaG)0.2:1.860334)0.4:0.500000,OutGroup);
	- Used in tests (expected outputs):
		- test_clean_and_save_trees_filters_low_support: cleaned == 1, dropped == 1, and 2 in dropped

- simple_species_tree
	- Input file name: species.tree
	- File content (Newick):
		- (((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4):0.5,(TaxaD:0.6,OutGroup:0.7):0.8);
	- Used in tests (expected outputs):
		- test_read_tree_file: len(trees) == 1 and get_taxa_from_tree(...) is non-empty
		- test_extract_triplet_subtree_missing_taxa: subtree is None
		- test_process_gene_trees_for_triplets: len(mapping) == len(triplets) and any value length > 0
		- test_clean_and_save_trees: output exists, len(cleaned) == 1, dropped == {}
		- test_format_newick_with_precision: newick endswith ";"

- simple_gene_trees
	- Input file name: genes.tree
	- File content (3 lines):
		- ((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
		- ((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
		- ((TaxaB:0.12,TaxaC:0.23):0.34,TaxaD:0.45);
	- Used in tests (expected outputs):
		- test_process_gene_trees_for_triplets: len(mapping) == len(triplets) and any value length > 0

- triplet_comparison_cases
	- Input cases (newick_str, triplet):
		- ((A:1.0,B:1.0):2.0,C:3.0,D:4.0); with triplet (A, B, C)
		- (((A:0.1,X:0.1):0.2,(B:0.1,Y:0.1):0.2):0.3,(C:0.1,Z:0.1):0.4); with triplet (A, B, C)
		- ((A:0.5,B:0.5):0.5,(C:0.2,D:0.2):0.8); with triplet (A, B, C)
	- Used in tests (expected outputs):
		- test_triplet_branch_lengths_match: for each case and pair (A,B), (A,C), (B,C), bio_dist == approx(dendro_dist)

## tests/test_tree_parser.py

- test_read_tree_file_single_tree
	- Inputs: simple_newick_file
	- Output: len(trees) == 1 and trees[0] is Tree
	- Purpose: verify single-tree parsing.

- test_read_tree_file_multiple_trees
	- Inputs: multiple_trees_file
	- Output: len(trees) == 3 and all are Tree
	- Purpose: verify multi-tree parsing.

- test_read_tree_file_not_found
	- Inputs: nonexistent path
	- Output: raises FileNotFoundError
	- Purpose: error handling for missing file.

- test_read_tree_file_invalid_newick
	- Inputs: invalid Newick text
	- Output: raises ValueError with "Invalid Newick format"
	- Purpose: error handling for invalid format.

- test_read_tree_file_random_text
	- Inputs: random invalid text
	- Output: raises ValueError with "Invalid Newick format"
	- Purpose: error handling for invalid format.

- test_read_tree_file_empty_file
	- Inputs: empty file
	- Output: raises ValueError with "Invalid Newick format"
	- Purpose: error handling for empty file.

- test_calculate_average_support_with_values
	- Inputs: newick_with_support_file
	- Output: avg_support is not None and 0.93 < avg_support < 0.95
	- Purpose: support averaging with values.

- test_calculate_average_support_no_values
	- Inputs: simple_newick_file
	- Output: avg_support is None
	- Purpose: support averaging with no values.

- test_remove_support_values
	- Inputs: newick_with_support_file
	- Output: avg_support becomes None after removal
	- Purpose: remove support values.

- test_standardize_tree_removes_support
	- Inputs: newick_with_support_file
	- Output: avg_support is None
	- Purpose: standardization removes support.

- test_standardize_tree_preserves_branch_lengths
	- Inputs: simple_newick_file
	- Output: branch lengths unchanged within 1e-10
	- Purpose: standardization preserves branch lengths.

- test_format_newick_with_precision_trailing_zeros
	- Inputs: simple_newick_file
	- Output: no "0000000000" and endswith ";"
	- Purpose: formatting trims zeros.

- test_format_newick_with_precision_default_places
	- Inputs: simple_newick_file
	- Output: contains "(" and ")" and endswith ";"
	- Purpose: default precision formatting.

- test_format_newick_with_custom_precision
	- Inputs: simple_newick_file
	- Output: endswith ";" and each decimal part length <= 5
	- Purpose: custom precision formatting.

- test_write_clean_trees
	- Inputs: simple_newick_file
	- Output: output exists and len(output_trees) == len(trees) == 1
	- Purpose: write cleaned trees.

- test_write_clean_trees_multiple
	- Inputs: multiple_trees_file
	- Output: len(output_trees) == 3
	- Purpose: write multiple trees.

- test_clean_and_save_trees_filters_low_support
	- Inputs: low_support_tree_file
	- Output: len(cleaned) == 1, len(dropped) == 1, and 2 in dropped
	- Purpose: filter low-support trees.

- test_clean_and_save_trees_no_filters
	- Inputs: simple_newick_file
	- Output: len(cleaned) == 1 and len(dropped) == 0
	- Purpose: no filtering needed.

- test_clean_and_save_trees_creates_output_file
	- Inputs: simple_newick_file
	- Output: output exists and size > 0
	- Purpose: output file creation.

- test_get_taxa_from_tree
	- Inputs: simple_newick_file
	- Output: len(taxa) == 8, includes TaxaA and OutGroup, sorted
	- Purpose: taxa extraction.

- test_get_taxa_from_tree_correct_names
	- Inputs: simple_newick_file
	- Output: taxa == ["OutGroup", "TaxaA", "TaxaB", "TaxaC", "TaxaD", "TaxaE", "TaxaF", "TaxaG"]
	- Purpose: taxa correctness.

- test_generate_triplets_count
	- Inputs: taxa list
	- Output: len(triplets) == 4
	- Purpose: triplet count formula.

- test_generate_triplets_excludes_outgroup
	- Inputs: taxa list + outgroup
	- Output: each triplet excludes OutGroup
	- Purpose: outgroup exclusion.

- test_generate_triplets_excludes_multiple_outgroups
	- Inputs: taxa list + two outgroups
	- Output: each triplet excludes OutGroup1 and OutGroup2
	- Purpose: multi-outgroup exclusion.

- test_generate_triplets_outgroup_comma_separated
	- Inputs: taxa list + comma-separated outgroups
	- Output: each triplet excludes OutGroup1 and OutGroup2
	- Purpose: comma parsing.

- test_generate_triplets_outgroup_comma_separated_with_spaces
	- Inputs: taxa list + comma-separated outgroups with whitespace
	- Output: each triplet excludes OutGroup1 and OutGroup2
	- Purpose: whitespace-tolerant parsing.

- test_root_tree_on_multiple_outgroups_prunes_outgroup_clade
	- Inputs: species tree + outgroups
	- Output: missing == set(), pruned_tree is not None, excluded has Tanypteryx and Pantala, ingroup has Anax_walsinghami and not Tanypteryx
	- Purpose: MRCA rooting and prune.

- test_generate_triplets_content
	- Inputs: taxa list
	- Output: triplets_set == {("TaxaA", "TaxaB", "TaxaC")} and len == 1
	- Purpose: triplet content correctness.

- test_generate_triplets_large_set
	- Inputs: larger taxa list
	- Output: len(triplets) == 35 and len(set(triplets)) == 35
	- Purpose: scaling check.

- test_read_triplet_filter_file_parses_valid_and_skips_invalid
	- Inputs: triplets.txt with lines "TaxaA,TaxaB,TaxaC" and "TaxaA,TaxaB"
	- Output: triplets == [("TaxaA", "TaxaB", "TaxaC")] and invalid_lines == [(2, "TaxaA,TaxaB")]
	- Purpose: validate parsing and invalid line capture.

- test_filter_triplets_by_taxa_skips_missing_taxa
	- Inputs: triplets [(TaxaA, TaxaB, TaxaC), (TaxaA, TaxaB, TaxaX)] and taxa_set {TaxaA, TaxaB, TaxaC}
	- Output: kept == [("TaxaA", "TaxaB", "TaxaC")] and skipped == [(("TaxaA", "TaxaB", "TaxaX"), ["TaxaX"]) ]
	- Purpose: skip triplets with missing taxa.

- test_write_triplets_to_file
	- Inputs: triplet list
	- Output: file exists; lines == ["TaxaA,TaxaB,TaxaC", "TaxaA,TaxaB,TaxaD", "TaxaA,TaxaC,TaxaD"]
	- Purpose: write triplets.

- test_write_triplets_to_file_format
	- Inputs: triplet list
	- Output: content == "TaxaOne,TaxaTwo,TaxaThree"
	- Purpose: output format.

- test_write_triplets_to_file_empty
	- Inputs: empty triplets
	- Output: file exists and content == ""
	- Purpose: empty handling.

- test_get_clean_filename_simple
	- Inputs: path with extension
	- Output: "/path/to/tree_clean.nwk"
	- Purpose: name formatting.

- test_get_clean_filename_different_extension
	- Inputs: path with different ext
	- Output: "/path/to/mytrees_clean.txt"
	- Purpose: name formatting.

- test_get_clean_filename_no_extension
	- Inputs: path without ext
	- Output: "/path/to/treefile_clean"
	- Purpose: name formatting.

- test_integration_full_workflow
	- Inputs: simple_newick_file
	- Output: len(trees) == 1 and len(output_trees) == 1
	- Purpose: end-to-end flow.

- test_integration_triplets_workflow
	- Inputs: simple_newick_file
	- Output: len(taxa) == 8, len(triplets) == 35, len(lines) == 35
	- Purpose: triplet workflow.

- test_extract_triplet_subtree_all_taxa_present
	- Inputs: tree and triplet
	- Output: subtree not None, taxa == {TaxaA, TaxaB, TaxaC}, len(terminals) == 3
	- Purpose: successful extraction.

- test_extract_triplet_subtree_missing_taxa
	- Inputs: tree + missing taxa
	- Output: subtree is None
	- Purpose: missing taxa handling.

- test_extract_triplet_subtree_preserves_branch_lengths
	- Inputs: tree + triplet
	- Output: newick contains TaxaA, TaxaB, TaxaC
	- Purpose: length correctness.

- test_process_gene_trees_for_triplets
	- Inputs: gene trees + triplets
	- Output: len(mapping) == 3 and each triplet has 1 tree
	- Purpose: triplet extraction on gene trees.

- test_process_gene_trees_for_triplets_empty
	- Inputs: empty gene trees
	- Output: mapping for (TaxaX, TaxaY, TaxaZ) has length 0
	- Purpose: empty handling.

- test_write_triplet_gene_trees
	- Inputs: triplet mapping
	- Output: header "TaxaA,TaxaB,TaxaC\t2", blank line, two trees, blank line, separator line, and header "TaxaD,TaxaE,TaxaF\t1"
	- Purpose: output formatting.

- test_write_triplet_gene_trees_empty_triplet
	- Inputs: triplet with no trees
	- Output: header "TaxaA,TaxaB,TaxaC\t0" and next line blank
	- Purpose: empty triplet output.

- test_integration_full_triplet_extraction_workflow
	- Inputs: species tree + gene trees
	- Output: output exists; content contains TaxaA,TaxaB,TaxaC; TaxaA,TaxaC,TaxaD; TaxaB,TaxaC,TaxaD; and a tab
	- Purpose: full extraction integration.

- test_write_triplet_gene_trees_multiprocess_single_worker
	- Inputs: triplets + gene file
	- Output: worker_count == 1, output exists, triplets_with_trees > 0, total_subtrees > 0, content includes both triplets
	- Purpose: single-worker path.

- test_write_triplet_gene_trees_multiprocess_with_workers
	- Inputs: triplets + gene file
	- Output: worker_count == 2, output exists, triplets_with_trees > 0, total_subtrees > 0, content includes both triplets
	- Purpose: multiprocessing path.

- test_metrics_logger_context_manager
	- Inputs: metrics file
	- Output: file closed after context; content includes "Test message 1" and "Test message 2"
	- Purpose: context manager behavior.

- test_metrics_logger_file_not_opened_before_enter
	- Inputs: metrics file
	- Output: before enter file does not exist and metrics_file is None; inside context metrics_file is not None
	- Purpose: lazy open behavior.

- test_triplet_gene_trees_separator_format
	- Inputs: triplet mapping
	- Output: separator_count == 2 and separators are surrounded by blank or tree/header lines
	- Purpose: output formatting.

- test_multiprocessing_triplet_writer_handles_empty_triplets
	- Inputs: triplets + gene file
	- Output: content includes "TaxaX,TaxaY,TaxaZ\t0", triplets_with_trees == 1, total_subtrees > 0
	- Purpose: empty triplet handling.

## tests/test_dendropy.py

- test_read_tree_file
	- Inputs: simple_species_tree
	- Output: len(trees) == 1 and get_taxa_from_tree(...) is non-empty
	- Purpose: read tree validation.

- test_calculate_average_support
	- Inputs: tree with support values
	- Output: avg_support == approx(0.7)
	- Purpose: support calculation.

- test_extract_triplet_subtree_missing_taxa
	- Inputs: tree + missing taxa
	- Output: subtree is None
	- Purpose: missing taxa handling.

- test_process_gene_trees_for_triplets
	- Inputs: species tree, gene trees
	- Output: len(mapping) == len(triplets) and any value length > 0
	- Purpose: extraction pipeline.

- test_generate_triplets_multiple_outgroups
	- Inputs: taxa list + outgroups
	- Output: each triplet excludes OutGroup1 and OutGroup2
	- Purpose: multi-outgroup handling.

- test_generate_triplets_outgroup_comma_separated
	- Inputs: taxa list + comma-separated
	- Output: each triplet excludes OutGroup1 and OutGroup2
	- Purpose: parsing behavior.

- test_generate_triplets_outgroup_comma_separated_with_spaces
	- Inputs: taxa list + comma-separated with whitespace
	- Output: each triplet excludes OutGroup1 and OutGroup2
	- Purpose: whitespace-tolerant parsing.

- test_read_triplet_filter_file_parses_valid_and_skips_invalid
	- Inputs: triplets.txt with lines "TaxaA,TaxaB,TaxaC" and "TaxaA,TaxaB"
	- Output: triplets == [("TaxaA", "TaxaB", "TaxaC")] and invalid_lines == [(2, "TaxaA,TaxaB")]
	- Purpose: validate parsing and invalid line capture.

- test_filter_triplets_by_taxa_skips_missing_taxa
	- Inputs: triplets [(TaxaA, TaxaB, TaxaC), (TaxaA, TaxaB, TaxaX)] and taxa_set {TaxaA, TaxaB, TaxaC}
	- Output: kept == [("TaxaA", "TaxaB", "TaxaC")] and skipped == [(("TaxaA", "TaxaB", "TaxaX"), ["TaxaX"]) ]
	- Purpose: skip triplets with missing taxa.

- test_root_tree_on_multiple_outgroups_prunes_outgroup_clade
	- Inputs: tree + outgroups
	- Output: missing == set(), pruned_tree is not None, excluded has Tanypteryx and Pantala, ingroup has Anax_walsinghami and not Tanypteryx
	- Purpose: MRCA pruning.

- test_write_triplet_gene_trees_streaming
	- Inputs: triplets + gene file
	- Output: output exists; total_subtrees > 0; triplets_with_trees > 0; content includes "TaxaA,TaxaB,TaxaC"
	- Purpose: streaming writer.

- test_write_triplet_gene_trees_multiprocess_triplets_single_worker
	- Inputs: triplets + gene file
	- Output: worker_count == 1, output exists, total_subtrees > 0, triplets_with_trees > 0
	- Purpose: triplet-parallel single worker.

- test_write_triplet_gene_trees_multiprocess_accepts_list
	- Inputs: triplets + in-memory gene list
	- Output: worker_count == 1, output exists, total_subtrees > 0, triplets_with_trees > 0, and no chunk dirs remain
	- Purpose: list input handling.

- test_write_triplet_gene_trees_separator
	- Inputs: triplet mapping
	- Output: content includes "=" * 60
	- Purpose: output formatting.

- test_clean_and_save_trees
	- Inputs: species tree
	- Output: output exists, len(cleaned) == 1, dropped == {}
	- Purpose: cleaning pipeline.

- test_format_newick_with_precision
	- Inputs: tree
	- Output: newick endswith ";"
	- Purpose: format output.

## tests/test_triplet_equivalence.py

- test_triplet_branch_lengths_match (parametrized)
	- Inputs: triplet_comparison_cases, pairs (A,B), (A,C), (B,C)
	- Output: bio_dist == approx(dendro_dist) for each case and pair
	- Purpose: verify branch length equivalence between approaches.

## tests/test_parser.py

- test_generate_response
	- Inputs: mock prompt
	- Output: "Response to: Hello, Ghostparser!"
	- Purpose: basic parser behavior.
