# Test Suite Documentation

This suite provides canonical (non-redundant) coverage grouped by functional area.

## Running Tests

### Run all tests

```bash
pytest
```

### Run one test file

```bash
pytest tests/test_tree_parser.py
```

### Run one test

```bash
pytest tests/test_tree_parser.py::test_extract_triplet_subtree_all_taxa_present
```

### Run reference tests (SciPy parity checks)

```bash
pytest -m reference
```

Run non-reference tests:

```bash
pytest -m "not reference"
```

Reference-marked tests in `tests/test_triplet_processor.py`:
- `test_custom_chi_square_matches_scipy_reference_randomized`
- `test_custom_ks_matches_scipy_asymptotic_reference_randomized`
- `test_two_sample_ks_test_hybrid_uses_scipy_near_threshold`

## Shared Fixtures (`tests/fixtures.py`)

- `simple_newick_file`
  - Input fixture content:
    - `(TaxaA:0.001,(TaxaB:0.098,(((TaxaC:0.001,TaxaD:0.001):0.001,TaxaE:0.001):0.086,(TaxaF:0.001,TaxaG:0.001):0.032):0.001):0.012,OutGroup:0.558);`
  - Used in:
    - `test_read_tree_file_single_tree`
    - `test_calculate_average_support_no_values`
    - `test_standardize_tree_preserves_branch_lengths`
    - `test_format_newick_with_precision_trailing_zeros`
    - `test_format_newick_with_precision_default_places`
    - `test_format_newick_with_custom_precision`
    - `test_write_clean_trees`
    - `test_clean_and_save_trees_no_filters`
    - `test_clean_and_save_trees_creates_output_file`
    - `test_get_taxa_from_tree_correct_names`
    - `test_integration_full_workflow`
    - `test_integration_triplets_workflow`
  - Expectations:
    - parsing returns one valid tree
    - average support is `None` when no support labels are present
    - branch lengths are preserved through standardization
    - formatted Newick remains valid and precision behavior is respected
    - cleaned output files are created and non-empty
    - extracted taxa match the expected sorted taxa set
    - integration workflows produce expected counts and output rows

- `newick_with_support_file`
  - Input fixture content:
    - `(((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup)0.85;`
  - Used in:
    - `test_calculate_average_support_with_values`
    - `test_remove_support_values`
    - `test_standardize_tree_removes_support`
  - Expectations:
    - average support is computed in the expected range
    - support labels are removed by cleaning/standardization
    - tree remains parseable after support removal

- `multiple_trees_file`
  - Input fixture content:
    - 3 Newick trees (one per line) with shared taxa set and varying topology
  - Used in:
    - `test_read_tree_file_multiple_trees`
    - `test_write_clean_trees_multiple`
  - Expectations:
    - parser returns all trees in the file
    - cleaned write/read roundtrip preserves tree count

- `low_support_tree_file`
  - Input fixture content:
    - one high-support tree and one low-support tree
  - Used in:
    - `test_clean_and_save_trees_filters_low_support`
  - Expectations:
    - low-support tree is dropped
    - dropped index metadata is reported correctly

- `triplet_comparison_cases`
  - Input fixture content:
    - list of `(newick_str, triplet)` pairs for cross-library triplet checks
  - Used in:
    - `test_triplet_branch_lengths_match`
    - `test_triplet_collapse_consistency_dendropy_vs_biopython`
  - Expectations:
    - pairwise patristic distances match between DendroPy and BioPython-derived triplet subtrees
    - DendroPy triplet extraction remains numerically consistent with BioPython prune/collapse behavior

## Triplet Utilities (`tests/test_tree_parser.py`)

- `test_generate_triplets_multiple_outgroups`
  - Confirms list-form outgroups are excluded.
- `test_generate_triplets_outgroup_comma_separated_with_spaces`
  - Confirms comma-separated outgroups with whitespace are excluded.
- `test_read_triplet_filter_file_parses_valid_and_skips_invalid`
  - Confirms valid triplets are parsed and malformed lines are reported.
- `test_filter_triplets_by_taxa_skips_missing_taxa`
  - Confirms triplets with taxa outside available set are skipped.
- `test_write_triplets_to_file`
  - Confirms ordered comma-separated output formatting.
- `test_write_triplet_gene_trees_streaming`
  - Confirms streaming writer emits valid triplet sections.
- `test_write_triplet_gene_trees_multiprocess_triplets_single_worker`
  - Confirms triplet-parallel path with one worker works.
- `test_write_triplet_gene_trees_multiprocess_accepts_list`
  - Confirms list input support and chunk cleanup behavior.
- `test_format_newick_with_precision_triplet_parser`
  - Confirms formatted Newick terminates with `;`.

## Triplet Equivalence (`tests/test_tree_parser.py`)

- `test_triplet_branch_lengths_match` (parametrized)
  - Confirms pairwise patristic distances match between DendroPy-extracted subtree and BioPython distance evaluation.
- `test_triplet_collapse_consistency_dendropy_vs_biopython`
  - Confirms DendroPy triplet extraction is numerically consistent with an independent BioPython prune/collapse path for all three pairwise distances in each triplet.

## Triplet Processing Pipeline (`tests/test_triplet_processor.py`)

- `test_compute_tree_height_statistic_matches_definition`
- `test_classify_triplet_topology_string_for_all_three_topologies`
- `test_pearson_discordant_chi_square_balanced_counts_not_significant`
- `test_two_proportion_discordant_z_test_balanced_counts_not_significant`
- `test_custom_chi_square_matches_scipy_reference_randomized`
- `test_custom_ks_matches_scipy_asymptotic_reference_randomized`
- `test_run_triplet_pipeline_no_introgression_when_dct_not_significant`
- `test_run_triplet_pipeline_inflow_when_ks_not_significant`
- `test_run_triplet_pipeline_outflow_when_con_summary_higher`
- `test_run_triplet_pipeline_ghost_when_dis_summary_higher`
- `test_parse_analyze_and_write_pipeline_roundtrip_with_species_header`

These tests cover topology classification, discordant-count statistics, KS behavior, and final introgression classification outputs.

## Runtime Argument Resolution and Process Defaults

### Orchestrator (`tests/test_orchestrator.py`)

- `test_resolve_runtime_args_cli_custom_processes_preserved`
- `test_resolve_runtime_args_config_processes_preserved_when_set`
- `test_resolve_runtime_args_config_with_cli_warns_and_ignores`

### Tree Parser (`tests/test_tree_parser.py`)

- `test_resolve_runtime_args_tree_parser_cli_custom_processes_preserved`
- `test_resolve_runtime_args_tree_parser_config_warns_and_ignores`
- `test_resolve_runtime_args_tree_parser_config_defaults_processes_to_zero`

### Triplet Processor (`tests/test_triplet_processor.py`)

- `test_resolve_runtime_args_triplet_processor_cli_custom_processes_preserved`
- `test_resolve_runtime_args_triplet_processor_config_processes_preserved_when_set`
- `test_resolve_runtime_args_triplet_processor_config_defaults_processes_to_zero`

### Config loaders (`tests/test_config.py`)

- `test_load_orchestrator_config_defaults_processes_to_zero`
- `test_load_tree_parser_config_defaults_processes_to_zero`
- `test_load_triplet_processor_config_defaults_processes_to_zero`

These tests confirm omitted `processes` defaults to `0`, explicit values are preserved, and config mode precedence is enforced.
