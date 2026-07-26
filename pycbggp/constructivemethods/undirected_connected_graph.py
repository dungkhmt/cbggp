"""
Undirected Connected Graph Generator. Author: Nguyen Ngoc Tuan Anh
See ALGORITHMS.md for algorithm details.
"""

import sys
import os
import random

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph
from DSU import DSU


def check_feasibility(nb_nodes: int, nb_edges: int):
    """Returns (is_feasible: bool, error_message: str)."""
    V, E = nb_nodes, nb_edges
    if V <= 0:
        return False, f"V={V} must be >= 1"
    if E < 0:
        return False, f"E={E} must be >= 0"
    if V == 1:
        return (True, "") if E == 0 else (False, f"V=1 requires E=0, got E={E}")
    if E < V - 1:
        return False, f"E={E} < V-1={V-1}: not enough edges for connectivity"
    if E > V * (V - 1) // 2:
        return False, f"E={E} > max_edges={V*(V-1)//2}"
    return True, ""


def generate_undirected_connected_graph(nb_nodes: int, nb_edges: int):
    """
    Generate a random undirected connected graph with exactly V nodes and E edges.

    Args:
        nb_nodes: V >= 1
        nb_edges: V-1 <= E <= V*(V-1)/2

    Returns:
        Graph, or None if infeasible.
    """
    V, E = nb_nodes, nb_edges

    ok, msg = check_feasibility(V, E)
    if not ok:
        print(f"[ConnectedGraph] Infeasible: {msg}")
        return None

    if V == 1:
        return Graph(1)

    # Phase 1: Random spanning tree via DSU
    all_edges = [(u, v) for u in range(V) for v in range(u + 1, V)]
    random.shuffle(all_edges)

    dsu = DSU(V)
    tree_edges, extra_edges = [], []
    for u, v in all_edges:
        (tree_edges if dsu.union(u, v) else extra_edges).append((u, v))

    # Phase 2: Fill E-(V-1) extra edges (already shuffled)
    selected = tree_edges + extra_edges[:E - (V - 1)]

    # Phase 3: Random vertex permutation
    perm = list(range(V))
    random.shuffle(perm)

    G = Graph(V)
    for u, v in selected:
        G.AddEdge(perm[u], perm[v])
    return G
