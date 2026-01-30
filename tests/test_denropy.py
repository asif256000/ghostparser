"""Test cases for the denropy module (DendroPy-backed)."""

from pathlib import Path

import pytest

from ghostparser.denropy import (
    calculate_average_support,
    clean_and_save_trees,
    extract_triplet_subtree,
    format_newick_with_precision,
    generate_triplets,
    get_taxa_from_tree,
    process_gene_trees_for_triplets,
    read_tree_file,
    write_triplet_gene_trees,
)


@pytest.fixture()
def simple_species_tree(tmp_path):
    species_tree_str = "(((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4):0.5,(TaxaD:0.6,OutGroup:0.7):0.8);"
    species_file = tmp_path / "species.tree"
    species_file.write_text(species_tree_str)
    return species_file


@pytest.fixture()
def simple_gene_trees(tmp_path):
    gene_trees_str = """((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
((TaxaB:0.12,TaxaC:0.23):0.34,TaxaD:0.45);
"""
    gene_file = tmp_path / "genes.tree"
    gene_file.write_text(gene_trees_str)
    return gene_file


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
