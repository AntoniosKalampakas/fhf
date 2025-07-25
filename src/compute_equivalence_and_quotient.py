"""
Compute equivalence relations R_V and R_V^σ and build quotient graph.
"""

import networkx as nx

def compute_equivalence_classes(ssp_dict, sigma=1.0):
    nodes = set()
    for (u, v) in ssp_dict:
        nodes.add(u)
        nodes.add(v)
    parent = {v: v for v in nodes}

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    def union(u, v):
        pu, pv = find(u), find(v)
        if pu != pv:
            parent[pu] = pv

    for (u, v), paths_uv in ssp_dict.items():
        if (v, u) in ssp_dict:
            strength_uv = min(min([min([w for w in path]) for path in paths_uv]), 
                              min([min([w for w in path]) for path in ssp_dict[(v, u)]]))
            if strength_uv >= sigma:
                union(u, v)

    # group nodes into classes
    classes = {}
    for v in nodes:
        leader = find(v)
        classes.setdefault(leader, []).append(v)
    return list(classes.values())

def build_quotient_graph(G, equivalence_classes, mu):
    Q = nx.DiGraph()
    for i, C in enumerate(equivalence_classes):
        Q.add_node(i, mu=max(mu[v] for v in C))

    for i, Ci in enumerate(equivalence_classes):
        for j, Cj in enumerate(equivalence_classes):
            if i != j:
                max_weight = 0
                for u in Ci:
                    for v in Cj:
                        if G.has_edge(u, v):
                            max_weight = max(max_weight, G[u][v]['weight'])
                if max_weight > 0:
                    Q.add_edge(i, j, weight=max_weight)
    return Q
