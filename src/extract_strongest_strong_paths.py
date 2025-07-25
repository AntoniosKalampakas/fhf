"""
Extract strongest strong paths from a fuzzy directed graph.
"""

import networkx as nx

def compute_con_F(G, u, v):
    # Compute maximum path strength from u to v
    try:
        all_paths = nx.all_simple_paths(G, u, v)
        return max(min(G[u][v]['weight'] for u, v in zip(p, p[1:])) for p in all_paths)
    except:
        return 0

def is_strong_edge(G, u, v, con_uv):
    return G[u][v]['weight'] == con_uv

def extract_strongest_strong_paths(G):
    ssp = {}
    for u in G.nodes:
        for v in G.nodes:
            if u == v:
                continue
            con_uv = compute_con_F(G, u, v)
            if con_uv == 0:
                continue
            valid_paths = []
            for path in nx.all_simple_paths(G, u, v):
                edge_weights = [G[path[i]][path[i+1]]['weight'] for i in range(len(path)-1)]
                if min(edge_weights) == con_uv and all(w == compute_con_F(G, path[i], path[i+1]) for i, w in enumerate(edge_weights)):
                    valid_paths.append(path)
            if valid_paths:
                ssp[(u, v)] = valid_paths
    return ssp
