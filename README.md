# Fuzzy Hypergraph Modeling of Directed Trust Networks

This repository contains the code used in the paper:

**"Hypercompositional Modeling of Directed Fuzzy Networks via Strong Path Structures"**  
by Antonios Kalampakas (2025)

It implements vertex-based fuzzy path hyperoperations, threshold-based equivalence relations, and quotient fuzzy graphs on a Facebook-like trust network derived from private message exchanges.

## 🔍 Purpose

The repository demonstrates:
- How to compute strongest strong paths in a fuzzy directed graph
- How to evaluate the fuzzy path hyperoperation \( \odot_F \)
- How to define and compute equivalence relations \( R_V \), \( R_V^\sigma \)
- How to generate quotient fuzzy graphs
- How to calculate structural metrics such as reachability impact, relay score, clustering coefficient, and centralities

## 📁 Dataset

The dataset is included in `data/facebook-wosn-links-message-weight.txt` and is based on:
> T. Opsahl and P. Panzarasa (2009). "Clustering in Weighted Networks."  
> Dataset retrieved from [https://toreopsahl.com/datasets/](https://toreopsahl.com/datasets/)  
> Licensed for research use with attribution.

## ▶️ How to Run

1. Install dependencies:
```bash
pip install -r requirements.txt
