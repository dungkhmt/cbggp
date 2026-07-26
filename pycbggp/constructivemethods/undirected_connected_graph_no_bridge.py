"""
Undirected Connected Graph — No Bridge (2-Edge-Connected). Author: Nguyen Ngoc Tuan Anh
See ALGORITHMS.md for algorithm details.
"""

import sys
import os
import random

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph


def check_feasibility(nb_nodes: int, nb_edges: int):
    """Returns (is_feasible: bool, error_message: str)."""
    V, E = nb_nodes, nb_edges
    if V < 3:
        return (True, "") if V == 1 and E == 0 else (
            False, f"V={V} < 3: 2-edge-connected requires V >= 3")
    if E < V:
        return False, f"E={E} < V={V}: need at least a Hamiltonian cycle (E >= V)"
    if E > V * (V - 1) // 2:
        return False, f"E={E} > max_edges={V*(V-1)//2}"
    return True, ""


def generate_undirected_connected_graph_no_bridge(nb_nodes: int, nb_edges: int):
    """
    Generate a random 2-edge-connected (no bridge) graph with V nodes and E edges.

    Args:
        nb_nodes: V >= 3
        nb_edges: V <= E <= V*(V-1)/2

    Returns:
        Graph, or None if infeasible.
    """
    V, E = nb_nodes, nb_edges

    ok, msg = check_feasibility(V, E)
    if not ok:
        print(f"[NoBridgeGraph] Infeasible: {msg}")
        return None

    if V == 1:
        return Graph(1)

    # Phase 1: Hamiltonian cycle on random permutation
    perm = list(range(V))
    random.shuffle(perm)

    edge_set = set()
    for i in range(V):
        u, v = perm[i], perm[(i + 1) % V]
        edge_set.add((min(u, v), max(u, v)))

    # Phase 2: Add E-V extra random edges
    if E - V > 0:
        candidates = [
            (u, v) for u in range(V) for v in range(u + 1, V)
            if (u, v) not in edge_set
        ]
        random.shuffle(candidates)
        for u, v in candidates[:E - V]:
            edge_set.add((u, v))

    G = Graph(V)
    for u, v in edge_set:
        G.AddEdge(u, v)
    return G
