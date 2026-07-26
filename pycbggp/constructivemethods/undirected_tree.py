"""
Random Undirected Tree Generator. Author: Nguyen Ngoc Tuan Anh
See ALGORITHMS.md for algorithm details (Prüfer sequence, uniform distribution).
"""

import sys
import os
import random
import heapq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph


def generate_undirected_tree(nb_nodes: int):
    """
    Generate a uniformly random labeled tree with V nodes (Prüfer sequence decode).

    Args:
        nb_nodes: V >= 1

    Returns:
        Graph, or None if V < 1.
    """
    V = nb_nodes

    if V < 1:
        print(f"[Tree] Infeasible: V={V} must be >= 1")
        return None
    if V == 1:
        return Graph(1)
    if V == 2:
        G = Graph(2); G.AddEdge(0, 1); return G

    # Step 1: Random Prüfer sequence of length V-2
    prufer = [random.randint(0, V - 1) for _ in range(V - 2)]

    # Step 2: Decode
    degree = [1] * V
    for node in prufer:
        degree[node] += 1

    leaves = []
    for v in range(V):
        if degree[v] == 1:
            heapq.heappush(leaves, v)

    edges = []
    for p in prufer:
        leaf = heapq.heappop(leaves)
        edges.append((leaf, p))
        degree[leaf] -= 1
        degree[p] -= 1
        if degree[p] == 1:
            heapq.heappush(leaves, p)

    last_two = [v for v in range(V) if degree[v] == 1]
    assert len(last_two) == 2
    edges.append((last_two[0], last_two[1]))

    # Step 3: Random vertex permutation
    perm = list(range(V))
    random.shuffle(perm)

    G = Graph(V)
    for u, v in edges:
        G.AddEdge(perm[u], perm[v])
    return G
