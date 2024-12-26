from src.phylotree import create_phylogenetic_tree
import numpy as np

fig = create_phylogenetic_tree("(((((ASV7:0.04339958,ASV10:0.03011500):0.02574618,ASV8:0.02691022):0.02048023,(ASV1:0.00233017,ASV4:0.00006684):0.03608202):0.03972852,ASV9:0.00796866):0.16163670,(ASV3:0.12500333,(ASV5:0.01925861,ASV6:0.00812806):0.08085625):0.06621116,ASV2:0.47242348);", show_labels=False)
fig.show()
