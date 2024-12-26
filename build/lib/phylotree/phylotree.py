from collections import OrderedDict
from plotly.graph_objs import Figure
import numpy as np
from Bio import Phylo
from io import StringIO

def create_phylogenetic_tree(
    newick_str,
    display_level=np.inf,
):
    """
    Function that returns a phylogenetic tree Plotly figure object.

    :param (str) newick_str: Newick formatted string. Polytomy is permissible.
    :param (int) display_level: The maximum level of the tree to display. The root is at level 0.

    Example 1: Simple unrooted tree with 4 leaves

    >>> from XXX import create_phylogenetic_tree

    >>> fig = create_phylogenetic_tree("(A,B,(C,D)E)F;")
    >>> fig.show()

    Example 2: Simple rooted tree with 4 leaves and distance values

    >>> from XXX import create_phylogenetic_tree

    >>> fig = create_phylogenetic_tree("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;")
    >>> fig.show()

    Example 3: Simple rooted tree with 4 leaves and distance values, but only display the first 2 levels after the root
    
    >>> from XXX import create_phylogenetic_tree

    >>> fig = create_phylogenetic_tree("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;", 2)
    >>> fig.show()
    """

    phylogenetic_tree = _Phylogenetic_Tree(
        newick_str,
        display_level,
    )

    return Figure(
        data=phylogenetic_tree.data, layout=phylogenetic_tree.layout
    )


class _Phylogenetic_Tree(object):
    """Refer to FigureFactory.create_phylogenetic_tree() for docstring."""

    def __init__(
        self,
        newick_str,
        display_level=np.inf,
        xaxis="xaxis",
        yaxis="yaxis",
    ):

        self.labels = None
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.data = []
        self.leaves = []
        self.sign = {self.xaxis: 1, self.yaxis: 1}
        self.layout = {self.xaxis: {}, self.yaxis: {}}
        self.display_level = display_level

        self.sign[self.xaxis] = -1
        self.sign[self.yaxis] = 1

        self.min_x = float("inf")
        self.max_x = float("-inf")
        self.min_y = float("inf")
        self.max_y = float("-inf")

        # Parse the Newick string
        newick_str = newick_str.replace(" ", "_")
        handle = StringIO(newick_str)
        tree = Phylo.read(handle, "newick")

        (dd_traces, ordered_labels, leaves) = self.get_phylo_tree_traces(
            tree, display_level
        )

        self.labels = ordered_labels
        self.leaves = leaves
        self.zero_vals = list(range(len(self.leaves)))

        self.layout = self.set_figure_layout()
        self.data = dd_traces

    def set_figure_layout(self):
        """
        Sets and returns default layout object for phylogenetic tree figure.

        """
        self.layout.update(
            {
                "showlegend": False,
                "autosize": True,
                "hovermode": "closest",
            }
        )

        self.set_axis_layout(self.xaxis)
        self.set_axis_layout(self.yaxis)

        total_width = self.max_x - self.min_x
        buffer = total_width * 0.5
        self.layout[self.xaxis].update(range=[self.min_x, self.max_x + buffer])

        total_height = self.max_y - self.min_y
        buffer = total_height * 0.05
        self.layout[self.yaxis].update(range=[self.max_y + buffer, self.min_y - buffer])

        return self.layout

    def set_axis_layout(self, axis_key):
        """
        Sets and returns default axis object for phylogenetic tree figure.

        :param (str) axis_key: E.g., 'xaxis', 'xaxis1', 'yaxis', yaxis1', etc.
        :rtype (dict): An axis_key dictionary with set parameters.

        """
        axis_defaults = {
            "type": "linear",
            "ticks": "",
            "mirror": "allticks",
            "rangemode": "tozero",
            "showticklabels": False,
            "zeroline": False,
            "showgrid": False,
            "showline": False,
            "automargin": True,
        }

        if len(self.labels) != 0:
            axis_key_labels = self.xaxis
            axis_key_labels = self.yaxis
            if axis_key_labels not in self.layout:
                self.layout[axis_key_labels] = {}
            self.layout[axis_key_labels]["tickvals"] = [
                zv * self.sign[axis_key] for zv in self.zero_vals
            ]
            self.layout[axis_key_labels]["ticktext"] = self.labels
            self.layout[axis_key_labels]["tickmode"] = "array"
            self.layout[axis_key_labels]["side"] = "right"

        # Set range to invert the y-axis
        if axis_key_labels == self.yaxis:
            self.layout[axis_key_labels]["range"] = [len(self.labels) + 1, -1]

        self.layout[axis_key].update(axis_defaults)

        return self.layout[axis_key]

    def get_phylo_tree_traces(self, tree, display_level):
        """
        Calculates all the elements needed for plotting a phylogenetic tree.

        :param (Bio.Phylo.BaseTree.Tree) tree: A Biopython Tree object parsed from a Newick formatted string.
        :rtype (tuple): Contains all the traces in the following order:
            (a) trace_list: List of Plotly trace objects for phylogenetic tree
            (b) xvals: All X points of the phylogenetic tree as array of arrays
                with length 4
            (c) yvals: All Y points of the phylogenetic tree as array of arrays
                with length 4
            (d) ordered_labels: leaf labels in the order they are going to
                appear on the plot
            (e) leaves: left-to-right traversal of the leaves

        """

        trace_list = []
        ordered_labels = []
        x_positions = {}
        y_positions = {}
        node_counter = 0

        # Helper function to get node name or name unnamed internal nodes
        def get_node_name(clade):
            nonlocal node_counter
            if clade.name:
                return clade.name
            else:
                node_counter += 1
                clade.name = f"internal_{node_counter}"
                return clade.name

        # Set root node position and name
        for clade in tree.find_clades(order="level"):
            if clade.name == "root":
                tree.root = clade
                break

        root_clade = tree.root
        if not root_clade.name:
            root_clade.name = "root"
        x_positions[get_node_name(root_clade)] = 0

        # Cut out unclassified clade if exists
        unclassified = None
        for clade in tree.root.clades:
            if clade.name == "unclassified":
                unclassified = clade
                tree.root.clades.remove(unclassified)
                tree.root.clades.extend(unclassified.clades)
                unclassified.clades = []
                break

        # Trim tree to display level
        def trim_tree_to_display_level(tree, display_level):
            if display_level == np.inf:
                return tree

            def trim_clade(clade, current_level):
                if current_level >= display_level:
                    clade.clades = []
                else:
                    for child in clade.clades:
                        trim_clade(child, current_level + 1)

            trim_clade(tree.root, 0)
            return tree

        tree = trim_tree_to_display_level(tree, display_level)

        # Collect all terminal nodes (leaves) and their y-positions
        terminals = tree.get_terminals()
        for idx, terminal in enumerate(terminals):
            node_name = get_node_name(terminal)
            self.leaves.append(node_name)
            ordered_labels.append(node_name)
            y_positions[node_name] = idx
            x_positions[node_name] = sum(
                clade.branch_length if clade.branch_length else 1
                for clade in tree.get_path(terminal)
            )

        # Make sure all clades have names to be referenced
        for clade in tree.find_clades(order="level"):
            get_node_name(clade)

        # Traverse tree to collect x- and y-positions for all internal nodes
        for clade in tree.find_clades(order="postorder"):
            node_name = get_node_name(clade)
            children = clade.clades
            if not clade.is_terminal():
                child_positions = [
                    y_positions[get_node_name(child)] for child in children
                ]
                y_positions[node_name] = sum(child_positions) / len(children)
            x_positions[node_name] = sum(
                clade.branch_length if clade.branch_length else 1
                for clade in tree.get_path(clade)
            )

        # Traverse tree in level order to produce traces
        for clade in tree.find_clades(order="level"):
            node_name = get_node_name(clade)
            x_node = x_positions[node_name]
            self.min_x = min(self.min_x, x_node)
            self.max_x = max(self.max_x, x_node)
            y_node = y_positions[node_name]
            self.min_y = min(self.min_y, y_node)
            self.max_y = max(self.max_y, y_node)

            children = clade.clades

            # Add a scatter trace for the node itself to display its name
            is_leaf = clade.is_terminal()
            trace_node = dict(
                type="scatter",
                x=[x_node],
                y=[y_node],
                mode="markers+text",
                text=[node_name] if is_leaf else [],
                hoverinfo="text",
                hovertext=[node_name] if not is_leaf else [],
                textposition="middle right",
                marker=dict(color="black", size=1),
            )
            trace_list.append(trace_node)

            for child in children:
                child_name = get_node_name(child)
                x0 = x_positions[node_name]
                x1 = x0 + (child.branch_length if child.branch_length else 1)
                y0 = y_positions[node_name]
                y1 = y_positions[child_name]

                trace1 = dict(
                    type="scatter",
                    x=[x0, x0],
                    y=[y0, y1],
                    mode="lines",
                    line=dict(color="black"),
                )

                trace2 = dict(
                    type="scatter",
                    x=[x0, x1],
                    y=[y1, y1],
                    mode="lines",
                    line=dict(color="black"),
                )

                trace_list.append(trace1)
                trace_list.append(trace2)

        # Insert Trace for unclassified clade if exists
        if unclassified:
            unclassified_name = get_node_name(unclassified)
            x_positions[unclassified_name] = 0
            y_positions[unclassified_name] = y_positions[root_clade.name] - 1
            trace_unclassified = dict(
                type="scatter",
                x=[x_positions[unclassified_name]],
                y=[y_positions[unclassified_name]],
                mode="markers+text",
                text=[unclassified_name],
                textposition="middle right",
                hoverinfo="text",
                hovertext=[unclassified_name],
                marker=dict(color="black", size=1),
            )
            trace_list.append(trace_unclassified)
            ordered_labels.append(unclassified_name)
            self.leaves.append(unclassified_name)

        return trace_list, ordered_labels, self.leaves