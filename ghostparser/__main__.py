"""Entry point for the ghostparser package.

This module allows the package to be executed with `python -m ghostparser`.
For tree parsing functionality, use `python -m ghostparser.tree_parser`.
"""

import sys


def main():
    """Display usage information for the ghostparser package."""
    print("GhostParser")
    print()
    print("Available modules:")
    print("  tree_parser    - Parse and standardize Newick format trees")
    print("  dendropy       - Parse and standardize Newick format trees")
    print()
    print("Usage:")
    print("  python -m ghostparser.tree_parser -st <species_tree> -gt <gene_trees> -og <outgroup>")
    print("  python -m ghostparser.dendropy -st <species_tree> -gt <gene_trees> -og <outgroup>")
    print()
    print("Example:")
    print("  python -m ghostparser.tree_parser -st species.nwk -gt genes.nwk -og Outgroup1")
    print("  python -m ghostparser.dendropy -st species.nwk -gt genes.nwk -og Outgroup1")
    print()
    print("For more information, run:")
    print("  python -m ghostparser.tree_parser --help")
    print("  python -m ghostparser.dendropy --help")


if __name__ == "__main__":
    main()
