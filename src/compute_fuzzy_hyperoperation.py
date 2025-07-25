"""
Compute the fuzzy path hyperoperation ⊛_F from strongest strong paths.
"""

def fuzzy_hyperoperation(ssp_dict, mu):
    hyperop_result = {}
    for (u, v), paths in ssp_dict.items():
        fuzzy_set = {}
        for path in paths:
            for node in path:
                fuzzy_set[node] = max(fuzzy_set.get(node, 0), mu.get(node, 0))
        hyperop_result[(u, v)] = fuzzy_set
    return hyperop_result
