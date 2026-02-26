# Test Suite Documentation

There are comprehensive test coverage with test cases organized into multiple categories.

## Running Tests

### Run All Tests

```bash
pytest
```

### Run Specific Test File

```bash
pytest tests/test_tree_parser.py
```

### Run Specific Test

```bash
pytest tests/test_tree_parser.py::test_extract_triplet_subtree_all_taxa_present
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

### Run SciPy Reference Comparison Tests

These tests are marked with `pytest.mark.reference` and validate custom statistical helpers against SciPy reference behavior.

Run only reference tests:

```bash
pytest -m reference
```

Run all tests except reference tests:

```bash
pytest -m "not reference"
```

Current reference-marked tests in `tests/test_triplet_processor.py`:

- `test_custom_chi_square_matches_scipy_reference_randomized`
- `test_custom_ks_matches_scipy_asymptotic_reference_randomized`
- `test_two_sample_ks_test_hybrid_uses_scipy_near_threshold`

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

## Triplet Utilities (tests/test_tree_parser.py)

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

- test_write_triplets_to_file
	- Inputs: triplets [(TaxaA, TaxaB, TaxaC), (TaxaA, TaxaB, TaxaD)]
	- Output: content == ["TaxaA,TaxaB,TaxaC", "TaxaA,TaxaB,TaxaD"]
	- Purpose: write triplet list file.

- test_extract_triplet_subtree_missing_taxa
	- Inputs: tree with taxa (TaxaA, TaxaB, TaxaC), triplet (TaxaA, TaxaB, TaxaX)
	- Output: subtree is None
	- Purpose: missing taxa handling.

- test_process_gene_trees_for_triplets
	- Inputs: three gene trees, three triplets
	- Output: len(mapping) == 3 and each triplet has 1 tree
	- Purpose: extraction pipeline.

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

- test_format_newick_with_precision_triplet_parser
	- Inputs: tree
	- Output: newick endswith ";"
	- Purpose: format output.

## Triplet Equivalence (tests/test_tree_parser.py)

- test_triplet_branch_lengths_match (parametrized)
	- Inputs: triplet_comparison_cases, pairs (A,B), (A,C), (B,C)
	- Output: bio_dist == approx(dendro_dist) for each case and pair
	- Purpose: verify branch length equivalence between approaches.


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

- test_clean_and_save_gene_trees_discards_missing_outgroup
	- Inputs: gene trees (one missing outgroup)
	- Output: len(cleaned) == 1, missing_indices == [2], rooted_count == 1
	- Purpose: discard gene trees without outgroup.

- test_format_newick_with_precision
	- Inputs: tree
	- Output: newick endswith ";"
	- Purpose: format output.

## Triplet Processing Pipeline (tests/test_triplet_processor.py)

- test_compute_tree_height_statistic_matches_definition
	- Inputs: rooted triplet tree with known branch lengths
	- Output: `H(T)` equals average root-to-tip distance
	- Purpose: validate corrected tree-height computation.

- test_classify_triplet_topology_for_all_three_topologies
	- Inputs: rooted trees for `con`, `dis1`, and `dis2`
	- Output: correct topology labels for species triplet `(A,B,C)`
	- Purpose: verify topology classification.

- test_pearson_discordant_chi_square_balanced_counts_not_significant
	- Inputs: equal discordant counts
	- Output: `z=0`, `p=1`
	- Purpose: validate DCT baseline behavior.

- test_run_triplet_pipeline_no_introgression_when_dct_not_significant
	- Inputs: balanced discordant trees
	- Output: `classification == no_introgression`
	- Purpose: verify early stop after DCT.

- test_run_triplet_pipeline_inflow_when_ks_not_significant
	- Inputs: significant DCT with similar `H(T)` distributions
	- Output: `classification == inflow_introgression`
	- Purpose: verify THT non-significant branch.

- test_run_triplet_pipeline_outflow_when_con_median_higher
	- Inputs: significant DCT and KS, with `median(H_con) > median(H_dis)`
	- Output: `classification == outflow_introgression`
	- Purpose: verify outflow decision rule.

- test_run_triplet_pipeline_ghost_when_dis_median_higher
	- Inputs: significant DCT and KS, with `median(H_con) < median(H_dis)`
	- Output: `classification == ghost_introgression`
	- Purpose: verify ghost decision rule.

- test_parse_analyze_and_write_pipeline_roundtrip
	- Inputs: small `unique_triplets_gene_trees.txt` fixture
	- Output: parse + analyze + write TSV succeed with expected fields
	- Purpose: validate file-oriented pipeline workflow.

