import pytest


@pytest.fixture
def simple_newick_file(tmp_path):
    """Create a temporary file with a simple Newick tree."""
    newick_str = "(TaxaA:0.001,(TaxaB:0.098,(((TaxaC:0.001,TaxaD:0.001):0.001,TaxaE:0.001):0.086,(TaxaF:0.001,TaxaG:0.001):0.032):0.001):0.012,OutGroup:0.558);"
    tree_file = tmp_path / "simple_tree.nwk"
    tree_file.write_text(newick_str)
    return tree_file


@pytest.fixture
def newick_with_support_file(tmp_path):
    """Create a temporary file with Newick tree containing support values."""
    newick_str = "(((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup)0.85;"
    tree_file = tmp_path / "tree_with_support.nwk"
    tree_file.write_text(newick_str)
    return tree_file


@pytest.fixture
def multiple_trees_file(tmp_path):
    """Create a temporary file with multiple Newick trees."""
    newick_lines = [
        "(TaxaA:0.001,(TaxaB:0.098,(TaxaC:0.001,TaxaD:0.001):0.001):0.012,OutGroup:0.558);",
        "(TaxaB:0.098,(TaxaC:0.001,TaxaD:0.001):0.001,OutGroup:0.558);",
        "((TaxaC:0.001,TaxaD:0.001):0.001,(TaxaB:0.098,TaxaA:0.001):0.012,OutGroup:0.558);",
    ]
    tree_file = tmp_path / "multiple_trees.nwk"
    tree_file.write_text("\n".join(newick_lines))
    return tree_file


@pytest.fixture
def low_support_tree_file(tmp_path):
    """Create a temporary file with trees having low average support."""
    newick_lines = [
        "(((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup);",
        "(((TaxaC,TaxaD)0.3:0.110599,(TaxaF,TaxaG)0.2:1.860334)0.4:0.500000,OutGroup);",
    ]
    tree_file = tmp_path / "low_support_trees.nwk"
    tree_file.write_text("\n".join(newick_lines))
    return tree_file


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


@pytest.fixture()
def triplet_comparison_cases():
    return [
        (
            "((A:1.0,B:1.0):2.0,C:3.0,D:4.0);",
            ("A", "B", "C"),
        ),
        (
            "(((A:0.1,X:0.1):0.2,(B:0.1,Y:0.1):0.2):0.3,(C:0.1,Z:0.1):0.4);",
            ("A", "B", "C"),
        ),
        (
            "((A:0.5,B:0.5):0.5,(C:0.2,D:0.2):0.8);",
            ("A", "B", "C"),
        ),
    ]
