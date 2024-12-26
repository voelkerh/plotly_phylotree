import pytest
from plotly.graph_objs import Figure
from phylotree import create_phylogenetic_tree

@pytest.mark.parametrize("newick_str", ["(A,B,(C,D)E)F;", "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;"])
def test_simple_tree(newick_str):
    fig = create_phylogenetic_tree(newick_str)
    assert isinstance(fig, Figure)

def test_display_level():
    fig = create_phylogenetic_tree("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;", 2)
    assert isinstance(fig, Figure)
