# Plotly PhyloTree Extension

This package extends Plotly to create phylogenetic trees from Newick strings.

## Installation

```bash
pip install plotly-phylotree
```

## Usage

```python
from phylotree import create_phylogenetic_tree

newick_str = "(A,(B,C)D)E;"
fig = create_phylogenetic_tree(newick_str)
fig.show()
```
![Alt text](/examples/output_images/basic_tree_labels.png "Basic Tree")