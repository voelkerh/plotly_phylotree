# Plotly PhyloTree Extension

This package extends Plotly to create interactive plots of phylogenetic trees from Newick strings. It can be used with plotly.py and can also be integrated into Dash Apps.

See also:
- Newick format - https://en.wikipedia.org/wiki/Newick_format
- Plotly - https://github.com/plotly/plotly.py

## Installation

```bash
pip install plotly-phylotree
```

## Usage

### 1. Basic Tree
- Labels shown by default
- No distances specified

```python
from phylotree import create_phylogenetic_tree

newick_str = "(A,(B,C)D)E;"
fig = create_phylogenetic_tree(newick_str)
fig.show()
```
![Alt text](/examples/output_images/basic_tree_labels.png "Basic Tree")

### 2. Tree with specified distances
- Labels shown by default
- Distances specified

```python
from phylotree import create_phylogenetic_tree

newick_str = "(Bovine:0.69395,(Gibbon:0.36079,(Orang-Utan:0.33636,(Gorilla:0.17147,(Chimp:0.19268,Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460)"
fig = create_phylogenetic_tree(newick_str)
fig.show()
```
![Alt text](/examples/output_images/mammals_tree_labels.png "Mammals")

### 3. Tree without labels
- Labels deactivated
- Distances specified

```python
from phylotree import create_phylogenetic_tree

newick_str = "(Bovine:0.69395,(Gibbon:0.36079,(Orang-Utan:0.33636,(Gorilla:0.17147,(Chimp:0.19268,Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460)"
fig = create_phylogenetic_tree(newick_str, show_labels=False)
fig.show()
```
![Alt text](/examples/output_images/mammals_tree_no_labels.png "Mammals without labels")

## Contributing

Suggestions for improvement are welcome. As this is a small side project, please allow some time for answers and revision.

## License

MIT License