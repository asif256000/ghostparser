# Orchestrator Module Guide

This document describes the end-to-end orchestrator in [ghostparser/orchestrator.py](ghostparser/orchestrator.py):

1. tree preprocessing/triplet extraction (from `tree_parser`)
2. per-triplet introgression inference (from `triplet_processor`)
3. final TSV reporting (`orchestrator_triplet_results.tsv`)

## CLI Input Options

Run:

```bash
python -m ghostparser.orchestrator -st <species_tree> -gt <gene_trees> -og <outgroup>
```

Or with a config file:

```bash
python -m ghostparser.orchestrator -c <config.yaml>
```

Options:

- `-st`, `--species_tree` (required)
  - Path to species tree file (Newick).
- `-gt`, `--gene_trees` (required)
  - Path to gene trees file (Newick, one tree per line).
- `-og`, `--outgroup` (required)
  - Outgroup taxon name(s). Multiple outgroups are comma-separated.
- `-c`, `--config-file` (optional)
  - Path to JSON/YAML config file.
  - When provided, other CLI options are ignored and a warning is printed.
- `-tf`, `--triplet-filter` (optional)
  - Path to triplet file (`taxon1,taxon2,taxon3` per line). If omitted, triplets are generated from species-tree ingroup taxa.
- `-o`, `--output` (optional)
  - Output folder path. Default is `./results` (from the current working directory).
- `-p`, `--processes` (optional)
  - Worker processes for both triplet extraction (`tree_parser`) and per-triplet inference (`triplet_processor`).
  - Defaults to `0`, which means all available CPU cores (`cpu_count()`).

`min_support_value` is supported via config file only (optional). If omitted, default tree-parser threshold behavior is used.

## Pipeline Summary

1. Species tree cleaning/rooting
   - Uses `clean_and_save_trees`, then roots/prunes with `_root_tree_on_outgroup`.
2. Triplet generation and normalization
   - Uses filter file or combinations from ingroup taxa.
   - Converts each triplet to A/B/C orientation where A and B are species sisters (`_build_species_triplet_metadata`).
3. Gene tree cleaning/rooting
   - Uses `clean_and_save_gene_trees`.
4. Triplet extraction file generation
  - Uses `write_triplet_gene_trees_multiprocess` from `tree_parser`.
  - Writes `unique_triplets_gene_trees.txt` in the same section format as `tree_parser`.
5. Per-triplet inference from file
  - Uses `analyze_triplet_gene_tree_file` from `triplet_processor` on `unique_triplets_gene_trees.txt`.
6. Final reporting
   - Uses `write_pipeline_results` to write `orchestrator_triplet_results.tsv`.

## Final TSV Columns (How Each Is Produced)

The final file is `orchestrator_triplet_results.tsv`.

### Identity and Mapping

- `triplet`
  - Source: `row["triplet"]` from `TripletPipelineResult.to_dict()`.
  - Method: tuple `(A, B, C)` joined as `A,B,C`.

- `abc_mapping`
  - Source: `TripletPipelineResult.to_dict()`.
  - Method: formatted as `A=<taxonA>,B=<taxonB>,C=<taxonC>` from the normalized triplet tuple.

### Topology Labels

Canonical topology strings:

- `((A,B),C)`
- `((B,C),A)`
- `((A,C),B)`

- `species_tree`
  - Source: species-triplet subtree header passed into `run_triplet_pipeline`.
  - Method: exact extracted species-triplet Newick string for that triplet.

- `con_topology`
  - Source: `run_triplet_pipeline`.
  - Method: topology that matches species tree for the triplet.

- `dis1_topology`
  - Source: `run_triplet_pipeline`.
  - Method: more frequent of the two discordant topologies.

- `dis2_topology`
  - Source: `run_triplet_pipeline`.
  - Method: less frequent of the two discordant topologies.

### Frequency Metadata

- `highest_freq_topologies`
  - Source: `run_triplet_pipeline`.
  - Method: set of topologies whose count equals the global max among the three.
  - TSV representation: comma-joined list.

- `most_frequent_matches_concordant`
  - Source: `run_triplet_pipeline`.
  - Method: `True` when `n_con >= n_dis1` and `n_con >= n_dis2`; otherwise `False`.

### Raw Topology Counts

Each gene tree in the triplet section is classified with `classify_triplet_topology_string`.

- `n_topology_ab`
  - Method: count of gene trees classified as `((A,B),C)`.

- `n_topology_bc`
  - Method: count of gene trees classified as `((B,C),A)`.

- `n_topology_ac`
  - Method: count of gene trees classified as `((A,C),B)`.

### Inference Group Counts

- `n_con`
  - Method: count for species-matching `con_topology`.

- `n_dis1`
  - Method: count for more frequent discordant topology.

- `n_dis2`
  - Method: count for less frequent discordant topology.

### DCT-like Test Outputs

- `dct_p_value`
  - Method: two-sided Z-test on discordant counts (`n_dis1`, `n_dis2`):
    - $z = (n_{dis1} - n_{dis2}) / \sqrt{n_{dis1}+n_{dis2}}$
    - $p = \mathrm{erfc}(|z|/\sqrt{2})$

- `dct_significant`
  - Method: `dct_p_value <= alpha_dct` (`alpha_dct` default `0.01`).

### Tree-Height Test Outputs

Per-gene-tree height uses:

- `compute_tree_height_statistic` = mean root-to-tip distance across 3 leaves
- for `((X:b2,Y:b3):b4,Z:b1)`: $H(T) = (b1 + b2 + b3 + 2b4)/3$

- `ks_statistic`
  - Method: two-sample KS statistic between height samples of `dis1_topology` vs `con_topology`.

- `ks_p_value`
  - Method: asymptotic KS p-value approximation from `_ks_asymptotic_p_value`.

- `ks_significant`
  - Method: `ks_p_value <= alpha_ks` (`alpha_ks` default `0.05`).

### Height Medians and Final Classification

- `median_con`
  - Method: median of heights for `con_topology`; only populated if KS is significant.

- `median_dis`
  - Method: median of heights for `dis1_topology`; only populated if KS is significant.

- `classification`
  - Decision path:
    - `no_introgression` if DCT is not significant
    - `outflow_introgression` if DCT significant but KS not significant
    - if KS significant:
      - `inflow_introgression` if `median_con > median_dis`
      - `ghost_introgression` if `median_con < median_dis`
      - `unresolved` if equal/undefined.

### Data Coverage Tracking

- `analyzed_trees`
  - Method (orchestrator output):
    - number of extracted triplet gene trees that are actually used in per-triplet inference.
  - This is the direct denominator behind per-triplet topology counts and statistics.

## Example TSV Representation

Below is an example of how one row appears in `orchestrator_triplet_results.tsv`.

Header (truncated for readability):

```tsv
triplet	abc_mapping	species_tree	con_topology	dis1_topology	dis2_topology	highest_freq_topologies	n_topology_ab	n_topology_bc	n_topology_ac	n_con	n_dis1	n_dis2	most_frequent_matches_concordant	dct_p_value	dct_significant	ks_statistic	ks_p_value	ks_significant	classification	analyzed_trees
```

Example data row:

```tsv
TaxaA,TaxaB,TaxaC	A=TaxaA,B=TaxaB,C=TaxaC	((TaxaA:1,TaxaB:1):1,TaxaC:1);	((A,B),C)	((B,C),A)	((A,C),B)	((B,C),A)	8	12	4	8	12	4	False	0.001	True	0.2	0.07	False	outflow_introgression	24
```

How to read this example quickly:

- `abc_mapping` shows exactly which taxa are A/B/C in all topology labels.
- `most_frequent_matches_concordant=False` because a discordant topology frequency is higher than concordant frequency.
- `classification = outflow_introgression` because DCT is significant, but KS is not significant.

## Additional Outputs from Orchestrator Run

- `processed_<species_tree_filename>`
- `processed_<gene_trees_filename>`
- `unique_triplets_gene_trees.txt`
- `orchestrator_triplet_results.tsv`
- `metrics.txt`

## Notes on Parallelism

- Multiprocessing is applied during triplet extraction via `write_triplet_gene_trees_multiprocess`.
- Multiprocessing is also applied during inference via `triplet_processor.analyze_triplet_gene_tree_file`.
- `-p 0` explicitly resolves to all cores via `cpu_count()`.
- Use `-p 1` for effective single-worker execution.
