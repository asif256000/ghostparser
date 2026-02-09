from io import StringIO

import dendropy
import pytest
from Bio import Phylo

from ghostparser.dendropy import extract_triplet_subtree as extract_triplet_subtree_dendro
from ghostparser.tree_parser import extract_triplet_subtree as extract_triplet_subtree_bio


def _dist_to_root(node):
    dist = 0.0
    while node is not None and node.edge is not None:
        if node.edge.length is not None:
            dist += node.edge.length
        node = node.parent_node
    return dist


def _dendro_distance(tree, a, b):
    tree.is_rooted = True
    node_a = tree.find_node_with_taxon_label(a)
    node_b = tree.find_node_with_taxon_label(b)
    mrca = tree.mrca(taxon_labels=[a, b])
    return _dist_to_root(node_a) + _dist_to_root(node_b) - 2.0 * _dist_to_root(mrca)


def _bio_distance(tree, a, b):
    return tree.distance(a, b)


@pytest.mark.parametrize("a,b", [("A", "B"), ("A", "C"), ("B", "C")])
def test_triplet_branch_lengths_match(triplet_comparison_cases, a, b):
    for newick_str, triplet in triplet_comparison_cases:
        bio_tree = Phylo.read(StringIO(newick_str), "newick")
        dendro_tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)

        bio_subtree = extract_triplet_subtree_bio(bio_tree, triplet)
        dendro_subtree = extract_triplet_subtree_dendro(dendro_tree, triplet)

        assert bio_subtree is not None
        assert dendro_subtree is not None

        bio_dist = _bio_distance(bio_subtree, a, b)
        dendro_dist = _dendro_distance(dendro_subtree, a, b)

        assert bio_dist == pytest.approx(dendro_dist)
