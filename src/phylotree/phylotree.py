"""
Module for creating and visualizing phylogenetic trees using Plotly.

This module provides functionality to convert Newick-formatted strings into 
interactive phylogenetic tree visualizations using Plotly.
"""

from io import StringIO
from plotly.graph_objs import Figure
import numpy as np
from Bio import Phylo

def create_phylogenetic_tree(
    newick_str,
    display_level=np.inf,
    show_labels=True
):
    """
    Function that returns a phylogenetic tree Plotly figure object.

    :param (str) newick_str: Newick formatted string. Polytomy is
        permissible.
    :param (int) display_level: The maximum level of the tree to display.
        The root is at level 0.
    :param (boolean) show_labels: Indicates, whether or not the tree labels
        will be rendered.

    Example 1: Simple unrooted tree with 4 leaves

    >>> from phylotree import create_phylogenetic_tree
    >>> fig = create_phylogenetic_tree("(A,B,(C,D)E)F;")
    >>> fig.show()

    Example 2: Simple rooted tree with 4 leaves and distance values
    >>> from phylotree import create_phylogenetic_tree
    >>> fig = create_phylogenetic_tree("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;")
    >>> fig.show()

    Example 3: Simple rooted tree with 4 leaves and distance values,
    but only display the first 2 levels after the root

    >>> fig = create_phylogenetic_tree("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;", 2)
    ... )
    >>> fig.show()
    """

    phylogenetic_tree = _PhylogeneticTree(
        newick_str,
        display_level,
        show_labels
    )

    return Figure(data=phylogenetic_tree.data, layout=phylogenetic_tree.layout)


class _PhylogeneticTree:
    """Internal class to create and manage phylogenetic tree visualization."""

    def __init__(
        self,
        newick_str,
        display_level=np.inf,
        show_labels=True,
        xaxis="xaxis",
        yaxis="yaxis",
    ):

        self.xaxis = xaxis
        self.yaxis = yaxis
        self.sign = {self.xaxis: -1, self.yaxis: 1}
        self.min_x = self.min_y = float("inf")
        self.max_x = self.max_y = float("-inf")

        tree = self._parse_newick(newick_str)
        self.display_level = display_level
        self.show_labels = show_labels
        self.layout = {self.xaxis: {}, self.yaxis: {}}

        # Generate data and layout for figure
        self.data = self.get_phylo_tree_traces(tree, display_level)
        self.layout = self.set_figure_layout()

    def _parse_newick(self, newick_str):
        """Parse Newick string into Biopython tree object.

        :param (str) newick_str: Newick formatted string
        :return: Bio.Phylo tree object
        :rtype (Bio.Phylo.BaseTree.Tree)
        """
        try:
            newick_str = newick_str.replace(" ", "_")
            handle = StringIO(newick_str)
            return Phylo.read(handle, "newick")
        except Exception as e:
            raise ValueError(f"Invalid Newick format: {e}") from e

    def set_figure_layout(self):
        """
        Sets and returns default layout object for phylogenetic tree figure.

        - Disables the legend (`showlegend=False`).
        - Enables automatic resizing (`autosize=True`).
        - Sets the hover mode to show the closest point (`hovermode='closest'`).
        - Configures the x- and y-axis with appropriate ranges
        and buffers to ensure the tree is fully visible.
        - Inverts the y-axis.
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

        # Extend x-axis and y-axis range with buffer to make space for figure
        total_width = self.max_x - self.min_x
        right_margin_buffer = total_width * 0.5  # larger to make space for leave labels
        left_margin_buffer = total_width * 0.05
        self.layout[self.xaxis].update(
            range=[self.min_x - left_margin_buffer, self.max_x + right_margin_buffer]
        )

        total_height = self.max_y - self.min_y
        upper_lower_margin_buffer = total_height * 0.05
        self.layout[self.yaxis].update(
            range=[
                self.max_y + upper_lower_margin_buffer,
                self.min_y - upper_lower_margin_buffer,
            ]
        )

        return self.layout

    def set_axis_layout(self, axis_key):
        """
        Configures the layout of the specified axis for the phylogenetic tree figure.

        :param (str) axis_key: E.g., 'xaxis', 'xaxis1', 'yaxis', yaxis1', etc.
        :rtype (dict): An axis_key dictionary with set parameters.
        """
        axis_defaults = {
            "ticks": "",
            "tickvals": [],
            "rangemode": "tozero",
            "zeroline": False,
            "showgrid": False,
            "automargin": True,
        }

        self.layout[axis_key].update(axis_defaults)

    def get_phylo_tree_traces(self, tree, display_level):
        """
        Calculates all traces needed for plotting a phylogenetic tree.

        :param (Bio.Phylo.BaseTree.Tree) tree: A Biopython Tree object parsed from a Newick string.
        :rtype (tuple): Contains all the traces in the following order:
            (a) trace_list: List of Plotly trace objects for phylogenetic tree
            (b) xvals: All X points of the phylogenetic tree as array of arrays
                with length 4
            (c) yvals: All Y points of the phylogenetic tree as array of arrays
                with length 4
        """

        x_positions = {}
        y_positions = {}
        node_counter = [0]

        root_clade = self.set_root_node(tree, x_positions, node_counter)
        unclassified = self.cut_unclassified_clade(tree)

        if not display_level == np.inf:
            self.trim_tree_to_display_level(tree.root, 0)

        # Make sure all nodes have names to be referenced
        for clade in tree.find_clades(order="level"):
            self.get_node_name(clade, node_counter)

        self.calculate_terminal_positions(tree, x_positions, y_positions, node_counter)
        self.calculate_internal_positions(tree, x_positions, y_positions, node_counter)

        trace_list = self.generate_traces(tree, x_positions, y_positions, node_counter)
        if unclassified:
            trace_unclassified = self.generate_trace_for_unclassified(
                unclassified, x_positions, y_positions, root_clade, node_counter
            )
            trace_list.append(trace_unclassified)

        return trace_list

    def set_root_node(self, tree, x_positions, node_counter):
        for clade in tree.find_clades(order="level"):
            if clade.name == "root":
                tree.root = clade
                break
        else:
            tree.root.name = tree.root.name or "root"
        x_positions[self.get_node_name(tree.root, node_counter)] = 0
        return tree.root

    # Helper function to get node name or name unnamed internal nodes
    def get_node_name(self, clade, node_counter):
        if clade.name:
            return clade.name
        else:
            node_counter[0] += 1
            clade.name = f"internal_{node_counter[0]}"
            return clade.name

    def cut_unclassified_clade(self, tree):
        for clade in tree.root.clades:
            if clade.name == "unclassified":
                unclassified = clade
                tree.root.clades.remove(unclassified)
                tree.root.clades.extend(unclassified.clades)
                unclassified.clades = []
                return unclassified
        return None

    def trim_tree_to_display_level(self, clade, current_level):
        if current_level >= self.display_level:
            clade.clades = []
        else:
            for child in clade.clades:
                self.trim_tree_to_display_level(child, current_level + 1)

    def calculate_terminal_positions(
        self, tree, x_positions, y_positions, node_counter
    ):
        terminals = tree.get_terminals()
        for idx, terminal in enumerate(terminals):
            node_name = self.get_node_name(terminal, node_counter)
            y_positions[node_name] = idx
            x_positions[node_name] = sum(
                clade.branch_length if clade.branch_length else 1
                for clade in tree.get_path(terminal)
            )

    def calculate_internal_positions(
        self, tree, x_positions, y_positions, node_counter
    ):
        for clade in tree.find_clades(order="postorder"):
            node_name = self.get_node_name(clade, node_counter)
            children = clade.clades
            if not clade.is_terminal():
                child_positions = [
                    y_positions[self.get_node_name(child, node_counter)]
                    for child in children
                ]
                y_positions[node_name] = sum(child_positions) / len(children)
            x_positions[node_name] = sum(
                clade.branch_length if clade.branch_length else 1
                for clade in tree.get_path(clade)
            )

    def generate_traces(self, tree, x_positions, y_positions, node_counter):
        trace_list = []

        for clade in tree.find_clades(order="level"):
            node_name = self.get_node_name(clade, node_counter)
            x_node = x_positions[node_name]
            y_node = y_positions[node_name]
            self.min_x = min(self.min_x, x_node)
            self.max_x = max(self.max_x, x_node)
            self.min_y = min(self.min_y, y_node)
            self.max_y = max(self.max_y, y_node)

            # Add scatter trace for node itself to display its name
            trace_node = self.generate_scatter_trace(
                x_node, y_node, node_name, (clade.is_terminal() and self.show_labels)
            )
            trace_list.append(trace_node)

            for child in clade.clades:
                child_name = self.get_node_name(child, node_counter)
                x0 = x_positions[node_name]
                x1 = x0 + (child.branch_length if child.branch_length else 1)
                y0 = y_positions[node_name]
                y1 = y_positions[child_name]

                vertical_trace = {
                    "type": "scatter",
                    "x": [x0, x0],
                    "y": [y0, y1],
                    "mode": "lines",
                    "line": {"color": "black"},
                }

                horizontal_trace = {
                    "type": "scatter",
                    "x": [x0, x1],
                    "y": [y1, y1],
                    "mode": "lines",
                    "line": {"color": "black"},
                }
                trace_list.append(vertical_trace)
                trace_list.append(horizontal_trace)

        return trace_list

    def generate_trace_for_unclassified(
        self, unclassified, x_positions, y_positions, root_clade, node_counter
    ):
        unclassified_name = self.get_node_name(unclassified, node_counter)
        x_positions[unclassified_name] = 0
        y_positions[unclassified_name] = 2 - max(y_positions.values())
        return self.generate_scatter_trace(
            x_positions[unclassified_name],
            y_positions[unclassified_name],
            unclassified_name,
            self.show_labels,
        )

    def generate_scatter_trace(self, x, y, text, is_leaf):
        return {
            "type": "scatter",
            "x": [x],
            "y": [y],
            "mode": "markers+text",
            "text": [text] if is_leaf else [],
            "hoverinfo": "text",
            "hovertext": [text] if not is_leaf else [],
            "textposition": "middle right",
            "marker": {"color": "black", "size": 1},
        }
