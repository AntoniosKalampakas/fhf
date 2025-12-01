#!/usr/bin/env python
# coding: utf-8
"""
fhf_pipeline.py

End-to-end pipeline for the Facebook-like fuzzy trust network.

It implements in one place:
- loading and log-normalizing the Facebook WOSN message weights
- extracting the largest strongly connected component
- restricting to the top-N weighted edges (for tractable hyperoperation)
- computing max–min connectivity strengths Con_F(u, v)
- extracting strong edges and enumerating strongest strong paths
- computing the vertex fuzzy path hyperoperation ⊙_F
- defining R_V^σ via mutual connectivity and building the quotient graph
- computing reachability impact and relay score on classes
- computing class-level centralities (aggregated from node-level):
    * betweenness, PageRank, in- and out-degree
- computing Spearman correlations (RI/RS vs centralities)
- generating:
    * Facebook_Top_Weighted_Subgraph.png
    * Quotient_Graph.png
    * Scat_RI_vs_PR.png
    * Scat_RS_vs_BC.png
    * Spear.png
- saving a CSV with all class-level metrics:
    * quotient_class_metrics.csv
"""

from __future__ import annotations

import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Hashable, Iterable, List, Mapping, Optional, Set, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------

DATA_FILE = os.path.join("data", "facebook-wosn-links-message-weight.txt")

# number of top-weighted edges used to build the working graph for
# strongest-paths, quotient, and metrics
TOP_N_EDGES = 600

# threshold for equivalence relation R_V^σ
SIGMA = 0.5

# maximum simple path length for strongest-strong enumeration
MAX_PATH_LENGTH = 5

# maximum number of paths per (u, v) to store
MAX_PATHS_PER_PAIR = 20

# random seed for layouts
RANDOM_SEED = 42

Node = Hashable
Path = List[Node]


# ---------------------------------------------------------------------
# 1. LOAD AND PREPROCESS
# ---------------------------------------------------------------------

def load_facebook_fuzzy_graph(path: str) -> nx.DiGraph:
    """
    Load the Facebook-like Social Network edge list and build a
    directed fuzzy graph with log-normalized edge weights in (0,1].
    """
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["source", "target", "weight"],
    )

    # Log-normalize weights
    df["weight_log"] = np.log1p(df["weight"])
    df["weight_log_norm"] = df["weight_log"] / df["weight_log"].max()

    G = nx.DiGraph()
    for _, row in df.iterrows():
        u = row["source"]
        v = row["target"]
        w = float(row["weight_log_norm"])
        G.add_edge(u, v, weight=w)

    return G


def largest_strongly_connected_subgraph(G: nx.DiGraph) -> nx.DiGraph:
    """
    Restrict G to its largest strongly connected component.
    """
    if nx.is_strongly_connected(G):
        return G.copy()
    largest_scc_nodes = max(nx.strongly_connected_components(G), key=len)
    return G.subgraph(largest_scc_nodes).copy()


def top_weighted_subgraph(G: nx.DiGraph, top_n: int) -> nx.DiGraph:
    """
    Take the top_n highest-weight edges and keep the largest weakly connected
    component of that subgraph. This is used as the working graph for
    strongest-path hyperoperation and quotient computations.
    """
    edges_sorted = sorted(
        G.edges(data=True),
        key=lambda x: x[2]["weight"],
        reverse=True,
    )
    top_edges = edges_sorted[:top_n]

    H = nx.DiGraph()
    for u, v, d in top_edges:
        H.add_edge(u, v, weight=d["weight"])

    if not nx.is_weakly_connected(H):
        largest_wcc_nodes = max(nx.weakly_connected_components(H), key=len)
        H = H.subgraph(largest_wcc_nodes).copy()

    return H


# ---------------------------------------------------------------------
# 2. MAX–MIN CONNECTIVITY AND STRONG EDGES
# ---------------------------------------------------------------------

def max_min_dijkstra(
    G: nx.DiGraph,
    source: Node,
    weight: str = "weight",
) -> Dict[Node, float]:
    """
    Single-source max–min (bottleneck) path strengths from `source`:

        Con_F(source, v) = max_{P: source->v} min_{e in P} λ(e).

    Uses a variant of Dijkstra's algorithm that propagates max-min capacity.
    """
    import heapq

    strengths: Dict[Node, float] = {v: 0.0 for v in G.nodes}
    strengths[source] = 1.0

    # max-heap via negative values
    heap: List[Tuple[float, Node]] = [(-1.0, source)]

    while heap:
        neg_strength, u = heapq.heappop(heap)
        current = -neg_strength

        if current < strengths[u]:
            continue

        for _, v, data in G.out_edges(u, data=True):
            w = float(data.get(weight, 0.0))
            if w <= 0.0:
                continue
            candidate = min(current, w)
            if candidate > strengths[v]:
                strengths[v] = candidate
                heapq.heappush(heap, (-candidate, v))

    return strengths


def compute_connectivity_strengths(
    G: nx.DiGraph,
    weight: str = "weight",
) -> Dict[Tuple[Node, Node], float]:
    """
    Compute Con_F(u, v) for all ordered pairs (u, v) in G.
    """
    con: Dict[Tuple[Node, Node], float] = {}
    for u in G.nodes:
        strengths = max_min_dijkstra(G, u, weight=weight)
        for v, val in strengths.items():
            con[(u, v)] = val
    return con


def extract_strong_edges(
    G: nx.DiGraph,
    con: Mapping[Tuple[Node, Node], float],
    weight: str = "weight",
) -> nx.DiGraph:
    """
    Strong-edge subgraph:

        (u, v) is strong iff λ(u, v) == Con_F(u, v) > 0.
    """
    G_str = nx.DiGraph()
    G_str.add_nodes_from(G.nodes(data=True))

    for u, v, data in G.edges(data=True):
        w = float(data.get(weight, 0.0))
        if w > 0.0 and abs(w - con.get((u, v), 0.0)) < 1e-12:
            G_str.add_edge(u, v, **data)

    return G_str


# ---------------------------------------------------------------------
# 3. STRONGEST STRONG PATHS AND HYPEROPERATION
# ---------------------------------------------------------------------

def enumerate_strongest_strong_paths(
    G_str: nx.DiGraph,
    con: Mapping[Tuple[Node, Node], float],
    max_paths_per_pair: Optional[int] = None,
    cutoff: Optional[int] = None,
    weight: str = "weight",
) -> Dict[Tuple[Node, Node], List[Path]]:
    """
    Enumerate strongest strong paths in the strong-edge subgraph.

    For each (u, v) with Con_F(u, v) > 0, we enumerate simple paths in G_str
    whose min edge weight equals Con_F(u, v).
    """
    paths: Dict[Tuple[Node, Node], List[Path]] = defaultdict(list)

    for u in G_str.nodes:
        for v in G_str.nodes:
            strength = con.get((u, v), 0.0)
            if strength <= 0.0 or u == v:
                continue

            try:
                candidate_paths = nx.all_simple_paths(
                    G_str,
                    source=u,
                    target=v,
                    cutoff=cutoff,
                )
            except nx.NetworkXNoPath:
                continue

            for P in candidate_paths:
                edge_weights = [
                    float(G_str[a][b].get(weight, 0.0))
                    for a, b in zip(P[:-1], P[1:])
                ]
                if not edge_weights:
                    continue

                path_strength = min(edge_weights)
                if abs(path_strength - strength) < 1e-12:
                    paths[(u, v)].append(P)
                    if (
                        max_paths_per_pair is not None
                        and len(paths[(u, v)]) >= max_paths_per_pair
                    ):
                        break

    return paths


@dataclass
class FuzzyHyperoperation:
    """
    Vertex fuzzy path hyperoperation ⊙_F:

      (u, v) ↦ fuzzy subset of vertices lying on strongest strong paths
                from u to v.
    """
    mu: Mapping[Node, float]
    support: Dict[Tuple[Node, Node], Set[Node]]

    def __call__(self, u: Node, v: Node) -> Dict[Node, float]:
        supp = self.support.get((u, v), set())
        return {w: float(self.mu.get(w, 0.0)) for w in supp}


def build_fuzzy_hyperoperation(
    G: nx.DiGraph,
    strongest_paths: Mapping[Tuple[Node, Node], Iterable[Path]],
    mu: Optional[Mapping[Node, float]] = None,
) -> FuzzyHyperoperation:
    """
    Construct ⊙_F from strongest strong paths.
    """
    if mu is None:
        mu = {v: 1.0 for v in G.nodes}

    support: Dict[Tuple[Node, Node], Set[Node]] = {}
    for (u, v), paths_uv in strongest_paths.items():
        nodes_on_paths: Set[Node] = set()
        for P in paths_uv:
            for w in P:
                nodes_on_paths.add(w)
        if nodes_on_paths:
            support[(u, v)] = nodes_on_paths

    return FuzzyHyperoperation(mu=mu, support=support)


# ---------------------------------------------------------------------
# 4. R_V^σ EQUIVALENCE AND QUOTIENT GRAPH
# ---------------------------------------------------------------------

def compute_equivalence_classes_from_con(
    G: nx.DiGraph,
    con: Mapping[Tuple[Node, Node], float],
    sigma: float,
) -> List[Set[Node]]:
    """
    Equivalence classes under R_V^σ:

        u ~ v  iff Con_F(u, v) ≥ σ and Con_F(v, u) ≥ σ.

    We build an undirected graph H whose edges connect such pairs
    and take its connected components as classes.
    """
    H = nx.Graph()
    H.add_nodes_from(G.nodes)

    nodes = list(G.nodes)
    n = len(nodes)
    for i in range(n):
        u = nodes[i]
        for j in range(i + 1, n):
            v = nodes[j]
            if (
                con.get((u, v), 0.0) >= sigma
                and con.get((v, u), 0.0) >= sigma
            ):
                H.add_edge(u, v)

    classes = [set(c) for c in nx.connected_components(H)]
    return classes


def build_quotient_graph(
    G: nx.DiGraph,
    classes: List[Set[Node]],
) -> Tuple[nx.DiGraph, Dict[Node, int]]:
    """
    Quotient graph:

      - each class C_i becomes a node i
      - edge weight between i and j is max of weights u->v
        with u in C_i, v in C_j, i != j.
    """
    Q = nx.DiGraph()
    node_to_class: Dict[Node, int] = {}

    for idx, cls in enumerate(classes):
        Q.add_node(idx, mu=1.0, size=len(cls))
        for v in cls:
            node_to_class[v] = idx

    for u, v, d in G.edges(data=True):
        cu = node_to_class[u]
        cv = node_to_class[v]
        if cu == cv:
            continue
        w = float(d.get("weight", 0.0))
        if Q.has_edge(cu, cv):
            if w > Q[cu][cv]["weight"]:
                Q[cu][cv]["weight"] = w
        else:
            Q.add_edge(cu, cv, weight=w)

    return Q, node_to_class


# ---------------------------------------------------------------------
# 5. METRICS: REACHABILITY IMPACT, RELAY, CENTRALITIES
# ---------------------------------------------------------------------

def compute_node_reachability_and_relay(
    hyperop: FuzzyHyperoperation,
) -> Tuple[Dict[Node, int], Dict[Node, int]]:
    """
    Node-level reachability impact and relay score from ⊙_F.

    reach[w] = number of pairs (u, v) for which w lies in u ⊙_F v
    relay[w]  = number of pairs (u, v) for which w lies in u ⊙_F v
                with w != u and w != v.
    """
    reach: Dict[Node, int] = defaultdict(int)
    relay: Dict[Node, int] = defaultdict(int)

    for (u, v), supp in hyperop.support.items():
        for w in supp:
            reach[w] += 1
            if w != u and w != v:
                relay[w] += 1

    return dict(reach), dict(relay)


def aggregate_class_metrics(
    G: nx.DiGraph,
    Q: nx.DiGraph,
    node_to_class: Dict[Node, int],
    node_reach: Dict[Node, int],
    node_relay: Dict[Node, int],
) -> pd.DataFrame:
    """
    Aggregate node-level metrics to class-level:

    For each class:
      - size = |C|
      - RI   = sum node_reach over members
      - RS   = sum node_relay over members
      - bc_max, bc_mean    (betweenness)
      - pr_max, pr_mean    (PageRank)
      - outdeg_sum, indeg_sum
    """
    # Node-level centralities on the working graph G
    bc = nx.betweenness_centrality(G)
    pr = nx.pagerank(G)
    outdeg = dict(G.out_degree())
    indeg = dict(G.in_degree())

    classes = sorted(Q.nodes())
    rows = []

    for c in classes:
        members = [v for v, cc in node_to_class.items() if cc == c]
        size = len(members)

        # reach/relay
        RI = sum(node_reach.get(v, 0) for v in members)
        RS = sum(node_relay.get(v, 0) for v in members)

        # betweenness
        bc_vals = [bc[v] for v in members]
        bc_max = max(bc_vals) if bc_vals else 0.0
        bc_mean = float(np.mean(bc_vals)) if bc_vals else 0.0

        # PageRank
        pr_vals = [pr[v] for v in members]
        pr_max = max(pr_vals) if pr_vals else 0.0
        pr_mean = float(np.mean(pr_vals)) if pr_vals else 0.0

        # degrees
        outdeg_sum = sum(outdeg.get(v, 0) for v in members)
        indeg_sum = sum(indeg.get(v, 0) for v in members)

        rows.append(
            {
                "class_id": c,
                "size": size,
                "RI": RI,
                "RS": RS,
                "bc_max": bc_max,
                "bc_mean": bc_mean,
                "pr_max": pr_max,
                "pr_mean": pr_mean,
                "outdeg_sum": outdeg_sum,
                "indeg_sum": indeg_sum,
            }
        )

    df = pd.DataFrame(rows).set_index("class_id").sort_index()
    return df


# ---------------------------------------------------------------------
# 6. PLOTS: NETWORK, QUOTIENT, SCATTER, HEATMAP
# ---------------------------------------------------------------------

def plot_top_weighted_graph(G_top: nx.DiGraph, out_file: str) -> None:
    """
    Layout and plot the top-N weighted subgraph.
    """
    plt.figure(figsize=(12, 10))
    pos = nx.spring_layout(G_top, seed=RANDOM_SEED, k=0.7)

    nx.draw_networkx_nodes(G_top, pos, node_color="skyblue", node_size=80)
    nx.draw_networkx_edges(
        G_top,
        pos,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=10,
        width=0.8,
        alpha=0.7,
    )

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


def plot_quotient_graph(Q: nx.DiGraph, df_classes: pd.DataFrame, out_file: str) -> None:
    """
    Plot the quotient graph with node size scaled by RI and color by RS.
    """
    plt.figure(figsize=(12, 10))
    pos = nx.spring_layout(Q, seed=RANDOM_SEED, k=1.2)

    classes = sorted(Q.nodes())
    RI = df_classes["RI"]
    RS = df_classes["RS"]

    # node sizes and colors
    sizes = [50 + 450 * (RI[c] / (RI.max() + 1e-9)) for c in classes]
    colors = [RS[c] for c in classes]

    nodes = nx.draw_networkx_nodes(
        Q,
        pos,
        nodelist=classes,
        node_size=sizes,
        node_color=colors,
        cmap="viridis",
    )
    nx.draw_networkx_edges(
        Q,
        pos,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=10,
        width=0.8,
        alpha=0.6,
    )
    nx.draw_networkx_labels(
        Q,
        pos,
        labels={c: f"C{c}" for c in classes},
        font_size=8,
    )

    cbar = plt.colorbar(nodes)
    cbar.set_label("Relay Score (RS)")

    plt.title("Quotient Graph F/R_V^σ (size=RI, color=RS)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


def plot_scatter_RI_vs_PR(df_classes: pd.DataFrame, out_file: str) -> None:
    plt.figure(figsize=(6, 5))
    x = df_classes["pr_max"]
    y = df_classes["RI"]

    plt.scatter(x, y, alpha=0.8)
    plt.xlabel("PageRank (max over class)")
    plt.ylabel("Reachability Impact (RI)")
    plt.title("Reachability Impact vs PageRank")
    plt.tight_layout()
    plt.savefig(out_file, dpi=600)
    plt.close()


def plot_scatter_RS_vs_BC(df_classes: pd.DataFrame, out_file: str) -> None:
    plt.figure(figsize=(6, 5))
    x = df_classes["bc_max"]
    y = df_classes["RS"]

    plt.scatter(x, y, alpha=0.8)
    plt.xlabel("Betweenness (max over class)")
    plt.ylabel("Relay Score (RS)")
    plt.title("Relay Score vs Betweenness")
    plt.tight_layout()
    plt.savefig(out_file, dpi=600)
    plt.close()


def plot_spearman_heatmap(df_classes: pd.DataFrame, out_file: str) -> None:
    """
    Spearman correlations heatmap for RI, RS vs class-level baselines.
    """
    metrics_for_RI_RS = [
        "bc_max",
        "bc_mean",
        "pr_max",
        "pr_mean",
        "outdeg_sum",
        "indeg_sum",
    ]

    cols = ["RI", "RS"] + metrics_for_RI_RS
    corr = df_classes[cols].corr(method="spearman")

    rows = ["RI", "RS"]
    heat_cols = metrics_for_RI_RS
    heat_data = corr.loc[rows, heat_cols].values

    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(heat_data, aspect="auto", cmap="viridis", origin="lower")
    fig.colorbar(im, ax=ax, label="Spearman ρ")

    ax.set_xticks(np.arange(len(heat_cols)))
    ax.set_yticks(np.arange(len(rows)))
    ax.set_xticklabels(heat_cols, rotation=45, ha="right")
    ax.set_yticklabels(rows)

    for i in range(len(rows)):
        for j in range(len(heat_cols)):
            ax.text(
                j,
                i,
                f"{heat_data[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=8,
            )

    ax.set_title("Class-level Spearman correlations (RI, RS vs baselines)")
    fig.tight_layout()
    plt.savefig(out_file, dpi=600)
    plt.close()


# ---------------------------------------------------------------------
# 7. MAIN PIPELINE
# ---------------------------------------------------------------------

def main():
    # 1) Load and restrict
    G_raw = load_facebook_fuzzy_graph(DATA_FILE)
    G_scc = largest_strongly_connected_subgraph(G_raw)
    G_work = top_weighted_subgraph(G_scc, TOP_N_EDGES)

    print(f"Loaded raw graph: {G_raw.number_of_nodes()} nodes, {G_raw.number_of_edges()} edges")
    print(f"Largest SCC: {G_scc.number_of_nodes()} nodes, {G_scc.number_of_edges()} edges")
    print(f"Working graph (top {TOP_N_EDGES} edges): {G_work.number_of_nodes()} nodes, {G_work.number_of_edges()} edges")

    # 2) Connectivity and strong edges
    con = compute_connectivity_strengths(G_work, weight="weight")
    G_str = extract_strong_edges(G_work, con, weight="weight")
    print(f"Strong-edge graph: {G_str.number_of_edges()} edges")

    # 3) Strongest strong paths and hyperoperation
    ssp = enumerate_strongest_strong_paths(
        G_str,
        con,
        max_paths_per_pair=MAX_PATHS_PER_PAIR,
        cutoff=MAX_PATH_LENGTH,
        weight="weight",
    )
    hyperop = build_fuzzy_hyperoperation(G_work, ssp, mu=None)

    # 4) Equivalence classes and quotient graph
    classes = compute_equivalence_classes_from_con(G_work, con, sigma=SIGMA)
    Q, node_to_class = build_quotient_graph(G_work, classes)

    print(f"Threshold σ={SIGMA}: {len(classes)} equivalence classes in quotient graph")
    print(f"Quotient graph: {Q.number_of_nodes()} nodes, {Q.number_of_edges()} edges")

    # 5) Node-level reachability and relay, aggregate to classes
    node_reach, node_relay = compute_node_reachability_and_relay(hyperop)
    df_classes = aggregate_class_metrics(
        G_work,
        Q,
        node_to_class,
        node_reach,
        node_relay,
    )

    # 6) Save class metrics
    df_classes.to_csv("quotient_class_metrics.csv")
    print("Saved class metrics to quotient_class_metrics.csv")

    # 7) Plots
    plot_top_weighted_graph(G_work, "Facebook_Top_Weighted_Subgraph.png")
    plot_quotient_graph(Q, df_classes, "Quotient_Graph.png")
    plot_scatter_RI_vs_PR(df_classes, "Scat_RI_vs_PR.png")
    plot_scatter_RS_vs_BC(df_classes, "Scat_RS_vs_BC.png")
    plot_spearman_heatmap(df_classes, "Spear.png")

    print("Generated figures: Facebook_Top_Weighted_Subgraph.png, Quotient_Graph.png,")
    print("Scat_RI_vs_PR.png, Scat_RS_vs_BC.png, Spear.png.")


if __name__ == "__main__":
    main()
