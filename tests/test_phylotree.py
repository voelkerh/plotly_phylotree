import pytest
from plotly.graph_objs import Figure
from phylotree import create_phylogenetic_tree
from phylotree.phylotree import _PhylogeneticTree

@pytest.mark.parametrize("newick_str", ["(A,B,(C,D)E)F;", "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;"])
def test_simple_tree(newick_str):
    fig = create_phylogenetic_tree(newick_str)
    assert isinstance(fig, Figure)

def test_display_level():
    fig = create_phylogenetic_tree("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;", 2)
    assert isinstance(fig, Figure)

def test_show_labels():
    fig = create_phylogenetic_tree("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;", show_labels=False)
    assert isinstance(fig, Figure)

def test_parse_newick_valid():
    tree = _PhylogeneticTree("(A,B,(C,D)E)F;")
    parsed_tree = tree._parse_newick("(A,B,(C,D)E)F;")
    assert parsed_tree is not None
    assert parsed_tree.root.name == "F"

def test_parse_newick_invalid():
    tree = _PhylogeneticTree("(A,B,(C,D)E)F;")
    with pytest.raises(ValueError):
        tree._parse_newick("(A,B,C,D)EF));")

def test_empty_newick():
    with pytest.raises(ValueError):
        create_phylogenetic_tree("")