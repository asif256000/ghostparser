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

## Test Categories

### File Reading Tests

Validate Newick file parsing and error handling:

- **test_read_tree_file_single_tree**: Verifies reading a single Newick tree from file
- **test_read_tree_file_multiple_trees**: Tests reading multiple trees from one file
- **test_read_tree_file_not_found**: Ensures proper FileNotFoundError for missing files
- **test_read_tree_file_invalid_newick**: Tests error handling for malformed Newick syntax (unbalanced parentheses)
- **test_read_tree_file_random_text**: Verifies rejection of non-Newick text
- **test_read_tree_file_empty_file**: Tests handling of empty input files

### Support Value Tests

Test calculation and removal of bootstrap/posterior probabilities:

- **test_calculate_average_support_with_values**: Verifies correct calculation of average support values from internal nodes
- **test_calculate_average_support_no_values**: Tests handling when no support values are present
- **test_remove_support_values**: Ensures all confidence values are removed from internal nodes

### Tree Standardization Tests

Verify support value removal while preserving branch lengths:

- **test_standardize_tree_removes_support**: Confirms support values are removed during standardization
- **test_standardize_tree_preserves_branch_lengths**: Ensures branch lengths remain intact after standardization

### Formatting Tests

Test high-precision branch length formatting:

- **test_format_newick_with_precision_trailing_zeros**: Verifies trailing zeros are removed (0.1000000000 → 0.1)
- **test_format_newick_with_precision_default_places**: Tests default 10 decimal place precision
- **test_format_newick_with_custom_precision**: Validates custom precision settings

### Tree Writing Tests

Validate writing of cleaned trees to output files:

- **test_write_clean_trees**: Tests writing a single cleaned tree
- **test_write_clean_trees_multiple**: Verifies writing multiple trees to file

### Quality Filtering Tests

Test filtering of low-quality trees:

- **test_clean_and_save_trees_filters_low_support**: Verifies trees with avg support < 0.5 are dropped
- **test_clean_and_save_trees_no_filters**: Tests that high-quality trees pass through
- **test_clean_and_save_trees_creates_output_file**: Ensures output file is created correctly

### Taxa Extraction Tests

Verify extraction of terminal taxa names:

- **test_get_taxa_from_tree**: Tests extraction of all terminal taxa from tree
- **test_get_taxa_from_tree_correct_names**: Validates taxa names are correctly identified

### Triplet Generation Tests

Test generation of all unique triplet combinations:

- **test_generate_triplets_count**: Verifies correct number of triplets (nC3 formula)
- **test_generate_triplets_excludes_outgroup**: Ensures outgroup is excluded from triplet generation
- **test_generate_triplets_content**: Validates triplet composition
- **test_generate_triplets_large_set**: Tests performance with larger taxa sets

### Triplet File Writing Tests

Test basic triplet file output:

- **test_write_triplets_to_file**: Verifies triplets are written to file correctly
- **test_write_triplets_to_file_format**: Tests comma-separated format
- **test_write_triplets_to_file_empty**: Handles empty triplet lists

### Triplet Subtree Extraction Tests

Verify extraction of triplet subtrees from gene trees:

- **test_extract_triplet_subtree_all_taxa_present**: Tests subtree extraction when all taxa exist
- **test_extract_triplet_subtree_missing_taxa**: Ensures None is returned when taxa are missing
- **test_extract_triplet_subtree_preserves_branch_lengths**: Validates branch lengths are recalculated correctly

### Gene Tree Processing Tests

Test processing of multiple gene trees:

- **test_process_gene_trees_for_triplets**: Verifies processing of gene trees for multiple triplets
- **test_process_gene_trees_for_triplets_empty**: Tests handling when no gene trees match triplets

### Triplet Gene Trees Output Tests

Test output file format:

- **test_write_triplet_gene_trees**: Validates complete output file format with tab-separated counts
- **test_write_triplet_gene_trees_empty_triplet**: Tests writing when a triplet has no matching gene trees

### File Naming Tests

Test output filename generation:

- **test_get_clean_filename_simple**: Tests adding \_clean prefix to filenames
- **test_get_clean_filename_different_extension**: Validates handling of various file extensions
- **test_get_clean_filename_no_extension**: Tests files without extensions

### Integration Tests

End-to-end workflow tests:

- **test_integration_full_workflow**: Complete pipeline test from reading to writing cleaned trees
- **test_integration_triplets_workflow**: Full triplet generation workflow
- **test_integration_full_triplet_extraction_workflow**: Complete end-to-end test including species tree, gene trees, triplet generation, and triplet gene tree extraction

## Test Fixtures

The test suite uses several fixtures for consistent test data:

- **simple_newick_file**: A basic tree with 8 taxa including OutGroup
- **newick_with_support_file**: Tree containing bootstrap/posterior probability values
- **multiple_trees_file**: File with multiple Newick trees for batch testing
- **low_support_tree_file**: Tree with average support below 0.5 threshold

## Test Data Conventions

All test data uses generic taxa names:

- Species: TaxaA, TaxaB, TaxaC, TaxaD, TaxaE, TaxaF, TaxaG
- Outgroup: OutGroup

This ensures tests are reproducible and don't depend on specific biological datasets.
