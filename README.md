This repository contains the complete reproducible code for the following papers:

1. Quotient-Based Structural Abstraction in Fuzzy Trust Networks via Strong Path Hyperoperations
This project develops a computational framework for:
vertex-based fuzzy path hyperoperations on directed fuzzy graphs
max–min connectivity strengths and strongest strong paths
a threshold-based equivalence relation 
quotient fuzzy graphs under mutual strongest connectivity
structural abstraction and influence analysis on a Facebook-like trust network derived from private message exchanges

The code implements:
log-normalization of trust weights
extraction of the largest strongly connected component
selection of a top-weighted subgraph
computation of ConF(u,v)
construction of the strong-edge subgraph
enumeration of strongest strong paths
the vertex fuzzy path hyperoperation 
the equivalence relation 
quotient fuzzy graph construction
Reachability Impact (RI) and Relay Score (RS)
class-level centrality aggregation
Spearman correlations and scatter plots


2. Strong-Path PageRank on Directed Fuzzy Graphs
This folder contains all code supporting the second paper, where PageRank is computed on quotient fuzzy graphs using strongest-path semantics.
It includes:
construction of coarse and fine quotient graphs
computation of class-level PageRank:
classical PageRank (PR)
weighted PageRank (WPR)
strong-path PageRank (SP-PR)
vertex membership extraction
clamped edge-membership construction
numerical solution of the strong-path PageRank system
scripts for generating all figures in the paper:
Figure 1: quotient graph visualization
Figure 2: PR vs WPR vs SP-PR bar chart
Figure 3: heatmap of differences SP-PR vs PR
All numerical values, PageRank vectors, and figures are produced directly from the repository code.
