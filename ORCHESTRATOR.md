# Orchestrator Module Guide

This document describes the end-to-end orchestrator in [ghostparser/orchestrator.py](ghostparser/orchestrator.py):

1. tree preprocessing/triplet extraction (from `tree_parser`)
2. per-triplet introgression inference (from `triplet_processor`)
3. final TSV aggregation (`orchestrator_triplet_results.tsv`)

## CLI Input Options

Run:

```bash
python -m ghostparser.orchestrator -st <species_tree> -gt <gene_trees> -og <outgroup>
```

Options:

- `-st`, `--species_tree` (required)
  - Path to species tree file (Newick).
- `-gt`, `--gene_trees` (required)
  - Path to gene trees file (Newick, one tree per line).
- `-og`, `--outgroup` (required)
  - Outgroup taxon name(s). Multiple outgroups are comma-separated.
- `-tf`, `--triplet-filter` (optional)
  - Path to triplet file (`taxon1,taxon2,taxon3` per line). If omitted, triplets are generated from species-tree ingroup taxa.
- `-o`, `--output` (optional)
  - Output folder relative to species-tree directory. Default is species-tree directory.
- `-p`, `--processes` (optional)
  - Worker processes for extraction + inference.
  - `0` means all available CPU cores (`cpu_count()`).
- `--no-multiprocessing` (optional flag)
  - Forces single-worker execution for extraction and inference.
- `--log-triplet-gene-trees` (optional flag)
  - Debug mode: appends generated triplets and extracted gene trees to `unique_triplets_gene_trees.txt` while processing.

## Pipeline Summary

1. Species tree cleaning/rooting
   - Uses `clean_and_save_trees`, then roots/prunes with `_root_tree_on_outgroup`.
2. Triplet generation and normalization
   - Uses filter file or combinations from ingroup taxa.
   - Converts each triplet to A/B/C orientation where A and B are species sisters (`_build_species_triplet_metadata`).
3. Gene tree cleaning/rooting
   - Uses `clean_and_save_gene_trees`.
4. In-memory triplet extraction + per-triplet inference
  - For each normalized triplet, gene-tree subtrees are extracted in memory from cleaned rooted gene trees.
  - `run_triplet_pipeline` is called immediately for that triplet (parallelized across triplets).
  - Optional debug logging writes the same section format as `tree_parser` into `unique_triplets_gene_trees.txt`.
5. Final reporting
   - Per-triplet dictionaries are written to TSV with `write_orchestrated_results_tsv`.

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

- `species_topology`
  - Source: `_species_topology_from_newick` in `triplet_processor`.
  - Method: classifies the species-triplet subtree header Newick using `classify_triplet_topology_string`.

- `con_topology`
  - Source: `run_triplet_pipeline`.
  - Method: always set equal to `species_topology` (species-anchored concordant definition).

- `con_topology_display`
  - Source: `TripletPipelineResult.to_dict()`.
  - Method:
    - starts as `con_topology`
    - appends `[diff]` if concordant topology is not among highest-frequency topologies
    - appends `[highest freq]` if concordant topology is tied/selected as highest frequency.

- `dis1_topology`
  - Source: `run_triplet_pipeline`.
  - Method: the more frequent of the two discordant topologies (relative to species topology). Tie: random order.

- `dis1_topology_display`
  - Source: `TripletPipelineResult.to_dict()`.
  - Method: `dis1_topology` with `[highest freq]` tag when applicable.

- `dis2_topology`
  - Source: `run_triplet_pipeline`.
  - Method: the other discordant topology after assigning `dis1_topology`.

- `dis2_topology_display`
  - Source: `TripletPipelineResult.to_dict()`.
  - Method: `dis2_topology` with `[highest freq]` tag when applicable.

### Frequency Ranking Metadata

- `top1_topology`, `top2_topology`, `top3_topology`
  - Source: `rank_topologies_by_frequency(topology_counts, rng=...)`.
  - Method: ranks all three topologies by descending count. Ties are shuffled randomly.

- `highest_freq_topologies`
  - Source: `run_triplet_pipeline`.
  - Method: set of topologies whose count equals the global max among the three.
  - TSV representation: comma-joined list.

- `concordant_diff`
  - Source: `run_triplet_pipeline`.
  - Method: boolean `con_topology not in highest_freq_topologies`.

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
  - Method: count for `con_topology` (species-concordant group).

- `n_dis1`
  - Method: count for `dis1_topology` (major discordant group).

- `n_dis2`
  - Method: count for `dis2_topology` (minor discordant group).

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

### Data Quality Tracking

- `skipped_trees`
  - Method: number of input triplet gene trees skipped in `run_triplet_pipeline` due to taxa mismatch or parse/classification errors.

## Example TSV Representation

Below is an example of how one row appears in `orchestrator_triplet_results.tsv`.

Header (truncated for readability):

```tsv
triplet	abc_mapping	species_topology	con_topology_display	dis1_topology_display	dis2_topology_display	top1_topology	top2_topology	top3_topology	highest_freq_topologies	n_topology_ab	n_topology_bc	n_topology_ac	n_con	n_dis1	n_dis2	dct_p_value	dct_significant	ks_statistic	ks_p_value	ks_significant	classification	skipped_trees
```

Example data row:

```tsv
TaxaA,TaxaB,TaxaC	A=TaxaA,B=TaxaB,C=TaxaC	((A,B),C)	((A,B),C) [diff]	((B,C),A) [highest freq]	((A,C),B)	((B,C),A)	((A,B),C)	((A,C),B)	((B,C),A)	8	12	4	8	12	4	0.001	True	0.2	0.07	False	outflow_introgression	0
```

How to read this example quickly:

- `abc_mapping` shows exactly which taxa are A/B/C in all topology labels.
- `con_topology_display` has `[diff]` because species-concordant topology is not highest frequency here.
- `dis1_topology_display` has `[highest freq]`, showing the dominant topology is discordant.
- `classification = outflow_introgression` because DCT is significant, but KS is not significant.

## Additional Outputs from Orchestrator Run

- `processed_<species_tree_filename>`
- `processed_<gene_trees_filename>`
- `metrics.txt`
- `orchestrator_triplet_results.tsv`

Optional (debug flag):

- `unique_triplets_gene_trees.txt` (when `--log-triplet-gene-trees` is enabled)

## Notes on Parallelism

- Triplet extraction and inference are combined in one per-triplet worker path in `analyze_triplets_parallel`.
- Multiprocessing is applied across triplets via `multiprocessing.Pool`.
- `-p 0` explicitly resolves to all cores via `cpu_count()`.
