"""
Visualize subgraphs and quotient graphs with Graphviz layout.
"""

import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

def draw_graph(G, node_sizes=None, node_colors=None, title="Graph", figsize=(12, 8)):
    pos = graphviz_layout(G, prog="neato")
    sizes = [node_sizes.get(n, 300) for n in G.nodes] if node_sizes else 300
    colors = [node_colors.get(n, 0.5) for n in G.nodes] if node_colors else "lightblue"

    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=colors, cmap=plt.cm.Blues)
    nx.draw_networkx_edges(G, pos, arrows=True, alpha=0.6)
    nx.draw_networkx_labels(G, pos, font_size=9)
    edge_labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels={k: f"{v:.2f}" for k, v in edge_labels.items()}, font_size=8)
    plt.title(title)
    plt.axis("off")
    plt.show()
