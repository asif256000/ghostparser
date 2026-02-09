"""Test cases for the dendropy module."""

import pytest

from ghostparser.dendropy import (
    calculate_average_support,
    clean_and_save_trees,
    extract_triplet_subtree,
    filter_triplets_by_taxa,
    format_newick_with_precision,
    generate_triplets,
    get_taxa_from_tree,
    _root_tree_on_outgroup,
    process_gene_trees_for_triplets,
    read_triplet_filter_file,
    read_tree_file,
    write_triplet_gene_trees_multiprocess,
    write_triplet_gene_trees_streaming,
    write_triplet_gene_trees,
)


def test_read_tree_file(simple_species_tree):
    trees = read_tree_file(str(simple_species_tree))
    assert len(trees) == 1
    assert get_taxa_from_tree(trees[0])


def test_calculate_average_support(tmp_path):
    tree_str = "((TaxaA:0.1,TaxaB:0.2)0.8:0.3,(TaxaC:0.4,TaxaD:0.5)0.6:0.7);"
    tree_file = tmp_path / "avg_support.tree"
    tree_file.write_text(tree_str)

    trees = read_tree_file(str(tree_file))
    avg_support = calculate_average_support(trees[0])
    assert avg_support == pytest.approx(0.7)


def test_extract_triplet_subtree_missing_taxa(simple_species_tree):
    tree = read_tree_file(str(simple_species_tree))[0]
    subtree = extract_triplet_subtree(tree, ("TaxaA", "TaxaB", "TaxaX"))
    assert subtree is None


def test_process_gene_trees_for_triplets(simple_species_tree, simple_gene_trees):
    species_tree = read_tree_file(str(simple_species_tree))[0]
    taxa = get_taxa_from_tree(species_tree)
    triplets = generate_triplets(taxa, "OutGroup")

    gene_trees = read_tree_file(str(simple_gene_trees))
    triplet_gene_trees = process_gene_trees_for_triplets(gene_trees, triplets)

    assert len(triplet_gene_trees) == len(triplets)
    assert any(len(trees) > 0 for trees in triplet_gene_trees.values())


def test_generate_triplets_multiple_outgroups():
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD", "OutGroup1", "OutGroup2"]
    triplets = generate_triplets(taxa, ["OutGroup1", "OutGroup2"])

    for triplet in triplets:
        assert "OutGroup1" not in triplet
        assert "OutGroup2" not in triplet


def test_generate_triplets_outgroup_comma_separated():
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD", "OutGroup1", "OutGroup2"]
    triplets = generate_triplets(taxa, "OutGroup1,OutGroup2")

    for triplet in triplets:
        assert "OutGroup1" not in triplet
        assert "OutGroup2" not in triplet


def test_generate_triplets_outgroup_comma_separated_with_spaces():
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD", "OutGroup1", "OutGroup2"]
    triplets = generate_triplets(taxa, "OutGroup1, OutGroup2")

    for triplet in triplets:
        assert "OutGroup1" not in triplet
        assert "OutGroup2" not in triplet


def test_read_triplet_filter_file_parses_valid_and_skips_invalid(tmp_path):
    filter_file = tmp_path / "triplets.txt"
    filter_file.write_text("TaxaA,TaxaB,TaxaC\nTaxaA,TaxaB\n\n")

    triplets, invalid_lines = read_triplet_filter_file(str(filter_file))

    assert triplets == [("TaxaA", "TaxaB", "TaxaC")]
    assert invalid_lines == [(2, "TaxaA,TaxaB")]


def test_filter_triplets_by_taxa_skips_missing_taxa():
    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaB", "TaxaX")]
    taxa_set = {"TaxaA", "TaxaB", "TaxaC"}

    kept, skipped = filter_triplets_by_taxa(triplets, taxa_set)

    assert kept == [("TaxaA", "TaxaB", "TaxaC")]
    assert skipped == [(("TaxaA", "TaxaB", "TaxaX"), ["TaxaX"])]


def test_root_tree_on_multiple_outgroups_prunes_outgroup_clade(tmp_path):
    newick_str = "((((((Ischnura,Platycnemis)1.0:0.527778,Hetaerina_americana)1.0:7.650645,Cloeon)1.0:1.756242,Epiophlebia)1.0:7.650645,(Tanypteryx,Pantala)1.0:7.650645),Anax_walsinghami);"
    tree_file = tmp_path / "species.tree"
    tree_file.write_text(newick_str)

    tree = read_tree_file(str(tree_file))[0]
    pruned_tree, excluded, missing, ingroup_taxa = _root_tree_on_outgroup(
        tree, ["Tanypteryx", "Pantala"]
    )

    assert missing == set()
    assert pruned_tree is not None
    assert "Tanypteryx" in excluded
    assert "Pantala" in excluded
    assert "Anax_walsinghami" in ingroup_taxa
    assert "Tanypteryx" not in ingroup_taxa


def test_write_triplet_gene_trees_streaming(tmp_path):
    gene_trees_str = """((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
((TaxaB:0.12,TaxaC:0.23):0.34,TaxaD:0.45);
"""
    gene_file = tmp_path / "genes.tree"
    gene_file.write_text(gene_trees_str)

    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaC", "TaxaD")]
    output_file = tmp_path / "triplet_gene_trees.txt"

    total_subtrees, triplets_with_trees = write_triplet_gene_trees_streaming(
        triplets, str(gene_file), str(output_file)
    )

    assert output_file.exists()
    assert total_subtrees > 0
    assert triplets_with_trees > 0
    content = output_file.read_text()
    assert "TaxaA,TaxaB,TaxaC" in content


def test_write_triplet_gene_trees_multiprocess_triplets_single_worker(tmp_path):
    gene_trees_str = """((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
((TaxaB:0.12,TaxaC:0.23):0.34,TaxaD:0.45);
"""
    gene_file = tmp_path / "genes.tree"
    gene_file.write_text(gene_trees_str)

    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaC", "TaxaD")]
    output_file = tmp_path / "triplet_gene_trees.txt"

    total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
        triplets,
        str(gene_file),
        str(output_file),
        use_multiprocessing=False,
    )

    assert worker_count == 1
    assert output_file.exists()
    assert total_subtrees > 0
    assert triplets_with_trees > 0


def test_write_triplet_gene_trees_multiprocess_accepts_list(tmp_path):
    gene_trees_newick = [
        "((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);",
        "((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);",
    ]
    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaC", "TaxaD")]
    output_file = tmp_path / "triplet_gene_trees.txt"

    total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
        triplets,
        gene_trees_newick,
        str(output_file),
        use_multiprocessing=False,
    )

    assert worker_count == 1
    assert output_file.exists()
    assert total_subtrees > 0
    assert triplets_with_trees > 0
    chunk_dirs = list(tmp_path.glob(".*_chunks"))
    assert not chunk_dirs


def test_write_triplet_gene_trees_separator(tmp_path):
    triplet_gene_trees = {
        ("TaxaA", "TaxaB", "TaxaC"): ["(TaxaA:0.1,TaxaB:0.2,TaxaC:0.3);"] ,
        ("TaxaD", "TaxaE", "TaxaF"): ["(TaxaD:0.4,TaxaE:0.5,TaxaF:0.6);"] ,
    }

    output_file = tmp_path / "triplet_gene_trees.txt"
    write_triplet_gene_trees(triplet_gene_trees, str(output_file))

    content = output_file.read_text()
    assert "=" * 60 in content


def test_clean_and_save_trees(simple_species_tree, tmp_path):
    output_file = tmp_path / "species_clean.tree"
    cleaned, dropped = clean_and_save_trees(str(simple_species_tree), str(output_file))

    assert output_file.exists()
    assert len(cleaned) == 1
    assert dropped == {}


def test_format_newick_with_precision(simple_species_tree):
    tree = read_tree_file(str(simple_species_tree))[0]
    newick = format_newick_with_precision(tree, decimal_places=5)
    assert newick.endswith(";")
