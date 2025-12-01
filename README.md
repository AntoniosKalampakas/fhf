This repository contains the complete reproducible code for the paper:

"Quotient-Based Structural Abstraction in Fuzzy Trust Networks via Strong Path Hyperoperations"  
by Antonios Kalampakas 

The project develops a computational framework for:
- vertex-based fuzzy path hyperoperations on directed fuzzy graphs
- max–min connectivity strengths and strongest strong paths
- a threshold-based equivalence relation R_V^σ
- quotient fuzzy graphs under mutual strongest connectivity
- structural abstraction and influence analysis on a Facebook-like trust network derived from private message exchanges

## Features

The code implements:
- log-normalization of trust weights
- extraction of the largest strongly connected component
- selection of a top-weighted subgraph
- computation of Con_F(u, v)
- construction of the strong-edge subgraph
- enumeration of strongest strong paths
- the vertex fuzzy path hyperoperation u ⊙_F v
- the equivalence relation R_V^σ
- quotient fuzzy graph construction
- Reachability Impact (RI) and Relay Score (RS)
- class-level centrality aggregation
- Spearman correlations and scatter plots

All results used in the paper are produced automatically.

## Repository structure

- fhf_pipeline.py – complete end-to-end pipeline
- data/ – input Facebook-like weighted edges
- src/ – optional modular implementation
- outputs/ – generated figures and metrics (created automatically)

## How to reproduce all results

From the repository root, run:

```bash
python fhf_pipeline.py
