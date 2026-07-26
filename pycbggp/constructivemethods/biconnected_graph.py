"""
Biconnected (2-Vertex-Connected) Undirected Graph Generator. Author: Nguyen Ngoc Tuan Anh
See ALGORITHMS.md for algorithm details.

Used for both:
  - gen_undirected_connected_graph_no_articulation_point
  - gen_undirected_connected_graph_no_bridge_no_articulation_point
(biconnected ⟺ no articulation point ⟹ no bridge)
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
        return False, f"V={V} < 3: biconnected requires V >= 3"
    if E < V:
        return False, f"E={E} < V={V}: biconnected requires E >= V"
    if E > V * (V - 1) // 2:
        return False, f"E={E} > max_edges={V*(V-1)//2}"
    return True, ""


def generate_biconnected_graph(nb_nodes: int, nb_edges: int):
    """
    Generate a random biconnected graph with V nodes and E edges.

    Args:
        nb_nodes: V >= 3
        nb_edges: V <= E <= V*(V-1)/2

    Returns:
        Graph, or None if infeasible.
    """
    V, E = nb_nodes, nb_edges

    ok, msg = check_feasibility(V, E)
    if not ok:
        print(f"[BiconnectedGraph] Infeasible: {msg}")
        return None

    edge_set = set()

    def add_edge(u, v):
        edge_set.add((min(u, v), max(u, v)))

    # Step 1: Triangle base
    add_edge(0, 1); add_edge(1, 2); add_edge(2, 0)

    in_graph = {0, 1, 2}
    free_vertices = list(range(3, V))
    random.shuffle(free_vertices)

    # Step 2: Open ear decomposition
    i = 0
    while i < len(free_vertices):
        remaining = len(free_vertices) - i
        ear_size = random.randint(1, min(remaining, 3))
        ear_verts = free_vertices[i:i + ear_size]
        i += ear_size

        u, w = random.sample(list(in_graph), 2)
        path = [u] + ear_verts + [w]
        for j in range(len(path) - 1):
            add_edge(path[j], path[j + 1])
        for v in ear_verts:
            in_graph.add(v)

    # Step 3: Fill extra edges
    need = E - len(edge_set)
    if need > 0:
        candidates = [
            (u, v) for u in range(V) for v in range(u + 1, V)
            if (u, v) not in edge_set
        ]
        random.shuffle(candidates)
        for u, v in candidates[:need]:
            edge_set.add((u, v))

    # Random vertex permutation
    perm = list(range(V))
    random.shuffle(perm)

    G = Graph(V)
    for u, v in edge_set:
        G.AddEdge(perm[u], perm[v])
    return G
