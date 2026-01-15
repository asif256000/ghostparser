"""Test cases for the tree_parser module."""

import tempfile
from pathlib import Path

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from ghostparser.tree_parser import (
    calculate_average_support,
    clean_and_save_trees,
    extract_triplet_subtree,
    format_newick_with_precision,
    generate_triplets,
    get_clean_filename,
    get_taxa_from_tree,
    process_gene_trees_for_triplets,
    read_tree_file,
    remove_support_values,
    standardize_tree,
    write_clean_trees,
    write_triplet_gene_trees,
    write_triplets_to_file,
)

# ============================================================================
# Fixtures
# ============================================================================


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
        "(((TaxaC,TaxaD)0.3:0.110599,(TaxaF,TaxaG)0.2:1.860334)0.4:0.500000,OutGroup);",  # Low support
    ]
    tree_file = tmp_path / "low_support_trees.nwk"
    tree_file.write_text("\n".join(newick_lines))
    return tree_file


# ============================================================================
# Tests for read_tree_file
# ============================================================================


def test_read_tree_file_single_tree(simple_newick_file):
    """Test reading a single tree from a Newick file."""
    trees = read_tree_file(str(simple_newick_file))
    assert len(trees) == 1
    assert isinstance(trees[0], Tree)


def test_read_tree_file_multiple_trees(multiple_trees_file):
    """Test reading multiple trees from a Newick file."""
    trees = read_tree_file(str(multiple_trees_file))
    assert len(trees) == 3
    assert all(isinstance(tree, Tree) for tree in trees)


def test_read_tree_file_not_found():
    """Test error handling for non-existent file."""
    with pytest.raises(FileNotFoundError):
        read_tree_file("nonexistent_file.nwk")


def test_read_tree_file_invalid_newick(tmp_path):
    """Test error handling for invalid Newick format."""
    invalid_file = tmp_path / "invalid.nwk"
    # Unbalanced parentheses should cause parsing error
    invalid_file.write_text("((TaxaA,TaxaB),TaxaC")

    with pytest.raises(ValueError, match="Invalid Newick format"):
        read_tree_file(str(invalid_file))


def test_read_tree_file_random_text(tmp_path):
    """Test error handling for random text that's not Newick format."""
    invalid_file = tmp_path / "random.nwk"
    invalid_file.write_text("(((,(")  # Unbalanced parentheses

    with pytest.raises(ValueError, match="Invalid Newick format"):
        read_tree_file(str(invalid_file))


def test_read_tree_file_empty_file(tmp_path):
    """Test error handling for empty file."""
    empty_file = tmp_path / "empty.nwk"
    empty_file.write_text("")

    with pytest.raises(ValueError, match="Invalid Newick format"):
        read_tree_file(str(empty_file))


# ============================================================================
# Tests for calculate_average_support
# ============================================================================


def test_calculate_average_support_with_values(newick_with_support_file):
    """Test calculating average support from a tree with support values."""
    trees = read_tree_file(str(newick_with_support_file))
    avg_support = calculate_average_support(trees[0])

    # Average of 0.95, 0.99, 0.98, 0.85 should be ~0.9425
    assert avg_support is not None
    assert 0.93 < avg_support < 0.95


def test_calculate_average_support_no_values(simple_newick_file):
    """Test calculating average support from a tree without support values."""
    trees = read_tree_file(str(simple_newick_file))
    avg_support = calculate_average_support(trees[0])

    # Tree has no support values, should return None
    assert avg_support is None


# ============================================================================
# Tests for remove_support_values
# ============================================================================


def test_remove_support_values(newick_with_support_file):
    """Test removing support values from a tree."""
    trees = read_tree_file(str(newick_with_support_file))
    tree = trees[0]

    # Verify tree has support values before
    original_support = calculate_average_support(tree)
    assert original_support is not None

    # Remove support values
    cleaned_tree = remove_support_values(tree)

    # Verify support values are gone
    new_support = calculate_average_support(cleaned_tree)
    assert new_support is None


# ============================================================================
# Tests for standardize_tree
# ============================================================================


def test_standardize_tree_removes_support(newick_with_support_file):
    """Test that standardize_tree removes support values."""
    trees = read_tree_file(str(newick_with_support_file))
    tree = trees[0]

    standardized = standardize_tree(tree)
    avg_support = calculate_average_support(standardized)

    assert avg_support is None


def test_standardize_tree_preserves_branch_lengths(simple_newick_file):
    """Test that standardize_tree preserves branch lengths."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    # Get original branch lengths
    original_lengths = [clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None]

    standardized = standardize_tree(tree)

    # Get standardized branch lengths
    standardized_lengths = [
        clade.branch_length for clade in standardized.find_clades() if clade.branch_length is not None
    ]

    assert len(original_lengths) == len(standardized_lengths)
    for orig, std in zip(original_lengths, standardized_lengths):
        assert abs(orig - std) < 1e-10


# ============================================================================
# Tests for format_newick_with_precision
# ============================================================================


def test_format_newick_with_precision_trailing_zeros(simple_newick_file):
    """Test that trailing zeros are removed in Newick output."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    newick_str = format_newick_with_precision(tree, decimal_places=10)

    # Should not have many zeros in a row
    assert "0000000000" not in newick_str
    # Should end with semicolon
    assert newick_str.endswith(";")


def test_format_newick_with_precision_default_places(simple_newick_file):
    """Test formatting with default decimal places."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    newick_str = format_newick_with_precision(tree)

    # Should be valid Newick format
    assert newick_str.endswith(";")
    assert "(" in newick_str
    assert ")" in newick_str


def test_format_newick_with_custom_precision(simple_newick_file):
    """Test formatting with custom decimal places."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    newick_str = format_newick_with_precision(tree, decimal_places=5)

    # Should be valid Newick format with custom precision
    assert newick_str.endswith(";")
    # Check that no more than 5 decimal places are used
    import re

    decimals = re.findall(r":0\.(\d+)", newick_str)
    for decimal_part in decimals:
        assert len(decimal_part) <= 5


# ============================================================================
# Tests for write_clean_trees
# ============================================================================


def test_write_clean_trees(simple_newick_file, tmp_path):
    """Test writing cleaned trees to a file."""
    trees = read_tree_file(str(simple_newick_file))
    output_file = tmp_path / "output_trees.nwk"

    write_clean_trees(trees, str(output_file))

    assert output_file.exists()

    # Verify the output file contains valid Newick
    output_trees = read_tree_file(str(output_file))
    assert len(output_trees) == len(trees)


def test_write_clean_trees_multiple(multiple_trees_file, tmp_path):
    """Test writing multiple cleaned trees to a file."""
    trees = read_tree_file(str(multiple_trees_file))
    output_file = tmp_path / "output_trees.nwk"

    write_clean_trees(trees, str(output_file))

    output_trees = read_tree_file(str(output_file))
    assert len(output_trees) == 3


# ============================================================================
# Tests for clean_and_save_trees
# ============================================================================


def test_clean_and_save_trees_filters_low_support(low_support_tree_file, tmp_path):
    """Test that clean_and_save_trees filters trees with low average support."""
    output_file = tmp_path / "cleaned_trees.nwk"

    cleaned, dropped = clean_and_save_trees(str(low_support_tree_file), str(output_file), min_avg_support=0.5)

    # Should keep high support tree and drop low support tree
    assert len(cleaned) == 1
    assert len(dropped) == 1
    assert 2 in dropped  # Second tree should be dropped


def test_clean_and_save_trees_no_filters(simple_newick_file, tmp_path):
    """Test clean_and_save_trees with trees that pass filter."""
    output_file = tmp_path / "cleaned_trees.nwk"

    cleaned, dropped = clean_and_save_trees(str(simple_newick_file), str(output_file), min_avg_support=0.0)

    assert len(cleaned) == 1
    assert len(dropped) == 0


def test_clean_and_save_trees_creates_output_file(simple_newick_file, tmp_path):
    """Test that clean_and_save_trees creates the output file."""
    output_file = tmp_path / "cleaned_trees.nwk"

    clean_and_save_trees(str(simple_newick_file), str(output_file))

    assert output_file.exists()
    assert output_file.stat().st_size > 0


# ============================================================================
# Tests for get_taxa_from_tree
# ============================================================================


def test_get_taxa_from_tree(simple_newick_file):
    """Test extracting taxa names from a tree."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    taxa = get_taxa_from_tree(tree)

    # Should have 8 taxa: TaxaA, TaxaB, TaxaC, TaxaD, TaxaE, TaxaF, TaxaG, OutGroup
    assert len(taxa) == 8
    assert "TaxaA" in taxa
    assert "OutGroup" in taxa
    # Taxa should be sorted
    assert taxa == sorted(taxa)


def test_get_taxa_from_tree_correct_names(simple_newick_file):
    """Test that correct taxa names are extracted."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    taxa = get_taxa_from_tree(tree)
    expected_taxa = sorted(["TaxaA", "TaxaB", "TaxaC", "TaxaD", "TaxaE", "TaxaF", "TaxaG", "OutGroup"])

    assert taxa == expected_taxa


# ============================================================================
# Tests for generate_triplets
# ============================================================================


def test_generate_triplets_count():
    """Test that correct number of triplets are generated."""
    taxa = ["A", "B", "C", "D", "E"]
    outgroup = "E"

    triplets = generate_triplets(taxa, outgroup)

    # 4C3 = 4 triplets
    assert len(triplets) == 4


def test_generate_triplets_excludes_outgroup():
    """Test that outgroup is excluded from triplets."""
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD", "OutGroup"]
    outgroup = "OutGroup"

    triplets = generate_triplets(taxa, outgroup)

    # Verify no triplet contains the outgroup
    for triplet in triplets:
        assert "OutGroup" not in triplet


def test_generate_triplets_content():
    """Test the content of generated triplets."""
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD"]
    outgroup = "TaxaD"

    triplets = generate_triplets(taxa, outgroup)
    triplets_set = set(triplets)

    expected = {
        ("TaxaA", "TaxaB", "TaxaC"),
    }

    assert len(triplets_set) == 1
    assert triplets_set == expected


def test_generate_triplets_large_set():
    """Test triplet generation with larger taxa set."""
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD", "TaxaE", "TaxaF", "TaxaG", "OutGroup"]  # 8 taxa
    outgroup = "OutGroup"

    triplets = generate_triplets(taxa, outgroup)

    # 7C3 = 35
    assert len(triplets) == 35
    # All triplets should be unique
    assert len(set(triplets)) == 35


# ============================================================================
# Tests for write_triplets_to_file
# ============================================================================


def test_write_triplets_to_file(tmp_path):
    """Test writing triplets to a file."""
    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaB", "TaxaD"), ("TaxaA", "TaxaC", "TaxaD")]
    output_file = tmp_path / "triplets.txt"

    write_triplets_to_file(triplets, str(output_file))

    assert output_file.exists()

    # Verify file content
    lines = output_file.read_text().strip().split("\n")
    assert len(lines) == 3
    assert lines[0] == "TaxaA,TaxaB,TaxaC"
    assert lines[1] == "TaxaA,TaxaB,TaxaD"
    assert lines[2] == "TaxaA,TaxaC,TaxaD"


def test_write_triplets_to_file_format(tmp_path):
    """Test the format of triplets in the file."""
    triplets = [("TaxaOne", "TaxaTwo", "TaxaThree")]
    output_file = tmp_path / "triplets.txt"

    write_triplets_to_file(triplets, str(output_file))

    content = output_file.read_text().strip()
    assert content == "TaxaOne,TaxaTwo,TaxaThree"


def test_write_triplets_to_file_empty(tmp_path):
    """Test writing empty triplet list."""
    triplets = []
    output_file = tmp_path / "triplets.txt"

    write_triplets_to_file(triplets, str(output_file))

    assert output_file.exists()
    assert output_file.read_text() == ""


# ============================================================================
# Tests for get_clean_filename
# ============================================================================


def test_get_clean_filename_simple():
    """Test generating clean filename."""
    filepath = "/path/to/tree.nwk"
    clean_filepath = get_clean_filename(filepath)

    assert clean_filepath == "/path/to/tree_clean.nwk"


def test_get_clean_filename_different_extension():
    """Test clean filename with different extension."""
    filepath = "/path/to/mytrees.txt"
    clean_filepath = get_clean_filename(filepath)

    assert clean_filepath == "/path/to/mytrees_clean.txt"


def test_get_clean_filename_no_extension():
    """Test clean filename for file without extension."""
    filepath = "/path/to/treefile"
    clean_filepath = get_clean_filename(filepath)

    assert clean_filepath == "/path/to/treefile_clean"


# ============================================================================
# Integration Tests
# ============================================================================


def test_integration_full_workflow(simple_newick_file, tmp_path):
    """Test the full workflow: read, clean, and write trees."""
    # Read tree
    trees = read_tree_file(str(simple_newick_file))
    assert len(trees) == 1

    # Standardize
    standardized = standardize_tree(trees[0])

    # Write
    output_file = tmp_path / "output.nwk"
    write_clean_trees([standardized], str(output_file))

    # Verify output
    output_trees = read_tree_file(str(output_file))
    assert len(output_trees) == 1


def test_integration_triplets_workflow(simple_newick_file, tmp_path):
    """Test the workflow: read tree, extract taxa, and generate triplets."""
    # Read tree
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    # Get taxa
    taxa = get_taxa_from_tree(tree)
    assert len(taxa) == 8

    # Generate triplets
    triplets = generate_triplets(taxa, "OutGroup")
    assert len(triplets) == 35  # 7C3

    # Write triplets
    output_file = tmp_path / "triplets.txt"
    write_triplets_to_file(triplets, str(output_file))

    assert output_file.exists()
    lines = output_file.read_text().strip().split("\n")
    assert len(lines) == 35


# ============================================================================
# Triplet Extraction Tests
# ============================================================================


def test_extract_triplet_subtree_all_taxa_present():
    """Test extracting a triplet subtree when all taxa are present."""
    from io import StringIO

    tree_str = "((TaxaA:0.1,TaxaB:0.2):0.3,(TaxaC:0.15,TaxaD:0.25):0.35);"
    tree = Phylo.read(StringIO(tree_str), "newick")

    triplet = ("TaxaA", "TaxaB", "TaxaC")
    subtree = extract_triplet_subtree(tree, triplet)

    assert subtree is not None
    subtree_taxa = {t.name for t in subtree.get_terminals()}
    assert subtree_taxa == {"TaxaA", "TaxaB", "TaxaC"}
    assert len(subtree.get_terminals()) == 3


def test_extract_triplet_subtree_missing_taxa():
    """Test extracting a triplet when not all taxa are present."""
    from io import StringIO

    tree_str = "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.15);"
    tree = Phylo.read(StringIO(tree_str), "newick")

    triplet = ("TaxaA", "TaxaB", "TaxaD")  # TaxaD is not in tree
    subtree = extract_triplet_subtree(tree, triplet)

    assert subtree is None


def test_extract_triplet_subtree_preserves_branch_lengths():
    """Test that branch lengths are preserved/adjusted correctly."""
    from io import StringIO

    tree_str = "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);"
    tree = Phylo.read(StringIO(tree_str), "newick")

    triplet = ("TaxaA", "TaxaB", "TaxaC")
    subtree = extract_triplet_subtree(tree, triplet)

    assert subtree is not None

    # Convert to Newick to verify structure
    newick = format_newick_with_precision(subtree)
    # Should have the three taxa with branch lengths
    assert "TaxaA" in newick
    assert "TaxaB" in newick
    assert "TaxaC" in newick


def test_process_gene_trees_for_triplets():
    """Test processing multiple gene trees for triplet extraction."""
    from io import StringIO

    # Create gene trees
    tree1_str = "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);"
    tree2_str = "((TaxaA:0.15,TaxaC:0.25):0.35,TaxaD:0.45);"
    tree3_str = "((TaxaB:0.12,TaxaC:0.22):0.32,TaxaD:0.42);"

    tree1 = Phylo.read(StringIO(tree1_str), "newick")
    tree2 = Phylo.read(StringIO(tree2_str), "newick")
    tree3 = Phylo.read(StringIO(tree3_str), "newick")

    gene_trees = [tree1, tree2, tree3]

    # Define triplets
    triplets = [
        ("TaxaA", "TaxaB", "TaxaC"),
        ("TaxaA", "TaxaC", "TaxaD"),
        ("TaxaB", "TaxaC", "TaxaD"),
    ]

    # Process
    triplet_gene_trees = process_gene_trees_for_triplets(gene_trees, triplets)

    # Verify results
    assert len(triplet_gene_trees) == 3

    # First triplet should have 1 tree (from tree1 only)
    assert len(triplet_gene_trees[("TaxaA", "TaxaB", "TaxaC")]) == 1

    # Second triplet should have 1 tree (from tree2 only)
    assert len(triplet_gene_trees[("TaxaA", "TaxaC", "TaxaD")]) == 1

    # Third triplet should have 1 tree (from tree3 only)
    assert len(triplet_gene_trees[("TaxaB", "TaxaC", "TaxaD")]) == 1


def test_process_gene_trees_for_triplets_empty():
    """Test processing when no gene trees match triplets."""
    from io import StringIO

    tree_str = "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);"
    tree = Phylo.read(StringIO(tree_str), "newick")

    gene_trees = [tree]
    triplets = [("TaxaX", "TaxaY", "TaxaZ")]  # None of these taxa exist

    triplet_gene_trees = process_gene_trees_for_triplets(gene_trees, triplets)

    assert len(triplet_gene_trees[("TaxaX", "TaxaY", "TaxaZ")]) == 0


def test_write_triplet_gene_trees(tmp_path):
    """Test writing triplet gene trees to file in the specified format."""
    triplet_gene_trees = {
        ("TaxaA", "TaxaB", "TaxaC"): [
            "(TaxaA:0.1,TaxaB:0.2,TaxaC:0.3);",
            "(TaxaA:0.15,TaxaB:0.25,TaxaC:0.35);",
        ],
        ("TaxaD", "TaxaE", "TaxaF"): [
            "(TaxaD:0.4,TaxaE:0.5,TaxaF:0.6);",
        ],
        ("TaxaG", "TaxaH", "TaxaI"): [],  # Empty triplet
    }

    output_file = tmp_path / "triplet_gene_trees.txt"
    write_triplet_gene_trees(triplet_gene_trees, str(output_file))

    assert output_file.exists()

    content = output_file.read_text()
    lines = content.split("\n")

    # Check first triplet header
    assert "TaxaA,TaxaB,TaxaC\t2" in lines[0]

    # Check blank line after header
    assert lines[1] == ""

    # Check first tree
    assert "TaxaA:0.1,TaxaB:0.2,TaxaC:0.3" in lines[2]

    # Check second tree
    assert "TaxaA:0.15,TaxaB:0.25,TaxaC:0.35" in lines[3]

    # Check double blank lines between triplets
    assert lines[4] == ""
    assert lines[5] == ""

    # Check second triplet header
    assert "TaxaD,TaxaE,TaxaF\t1" in lines[6]


def test_write_triplet_gene_trees_empty_triplet(tmp_path):
    """Test writing when a triplet has no gene trees."""
    triplet_gene_trees = {
        ("TaxaA", "TaxaB", "TaxaC"): [],
    }

    output_file = tmp_path / "triplet_gene_trees.txt"
    write_triplet_gene_trees(triplet_gene_trees, str(output_file))

    assert output_file.exists()

    content = output_file.read_text()
    lines = content.split("\n")

    # Should have header with count 0 and blank line
    assert "TaxaA,TaxaB,TaxaC\t0" in lines[0]
    assert lines[1] == ""
    # No trees after blank line


def test_integration_full_triplet_extraction_workflow(tmp_path):
    """Test the complete workflow from species tree to triplet gene trees."""
    # Create species tree file
    species_tree_str = "(((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4):0.5,(TaxaD:0.6,OutGroup:0.7):0.8);"
    species_file = tmp_path / "species.tree"
    species_file.write_text(species_tree_str)

    # Create gene trees file
    gene_trees_str = """((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
((TaxaB:0.12,TaxaC:0.23):0.34,TaxaD:0.45);
"""
    gene_file = tmp_path / "genes.tree"
    gene_file.write_text(gene_trees_str)

    # Read species tree and generate triplets
    species_trees = read_tree_file(str(species_file))
    taxa = get_taxa_from_tree(species_trees[0])
    triplets = generate_triplets(taxa, "OutGroup")

    # Read gene trees
    gene_trees = read_tree_file(str(gene_file))

    # Process gene trees for triplets
    triplet_gene_trees = process_gene_trees_for_triplets(gene_trees, triplets)

    # Write output
    output_file = tmp_path / "triplet_gene_trees.txt"
    write_triplet_gene_trees(triplet_gene_trees, str(output_file))

    # Verify output
    assert output_file.exists()
    content = output_file.read_text()

    # Check that we have triplets in the output
    assert "TaxaA,TaxaB,TaxaC" in content
    assert "TaxaA,TaxaC,TaxaD" in content
    assert "TaxaB,TaxaC,TaxaD" in content

    # Check that counts are present
    assert "\t" in content  # Tab-separated counts
