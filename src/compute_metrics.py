"""
Compute reachability impact, relay score, clustering, and centralities.
"""

import networkx as nx


def compute_reachability(hyperop_result):
    impact = {}
    for (u, v), fset in hyperop_result.items():
        for w in fset:
            impact[w] = impact.get(w, 0) + 1
    return impact


def compute_relay_score(hyperop_result):
    relay = {}
    for (u, v), fset in hyperop_result.items():
        for w in fset:
            if w != u and w != v:
                relay[w] = relay.get(w, 0) + 1
    return relay


def compute_clustering_and_transitivity(G):
    clustering = nx.average_clustering(G.to_undirected())
    transitivity = nx.transitivity(G)
    return clustering, transitivity


def compute_centralities(G):
    return {
        "in_degree": dict(G.in_degree()),
        "out_degree": dict(G.out_degree()),
        "betweenness": nx.betweenness_centrality(G),
        "closeness": nx.closeness_centrality(G),
        "pagerank": nx.pagerank(G),
    }
