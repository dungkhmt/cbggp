"""
Undirected Tree — Bounded Degree and Diameter. Author: Nguyen Ngoc Tuan Anh
See ALGORITHMS.md for algorithm details (slot-based depth-constrained attachment).
"""

import sys
import os
import random
from collections import deque

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph


def _max_subtree_size(depth: int, max_children: int) -> int:
    """Max nodes in a subtree of height <= depth with at most max_children per node."""
    if depth <= 0 or max_children <= 0:
        return 1
    if max_children == 1:
        return depth + 1
    return (max_children ** (depth + 1) - 1) // (max_children - 1)


def check_feasibility(nb_nodes: int, ub_deg: int, ub_diameter: int):
    """Returns (is_feasible: bool, error_message: str)."""
    V, D, diam = nb_nodes, ub_deg, ub_diameter

    if V < 1:
        return False, f"V={V} must be >= 1"
    if V == 1:
        return True, ""
    if V == 2:
        return (True, "") if D >= 1 and diam >= 1 else (False, "V=2 needs D>=1, diam>=1")
    if D < 1:
        return False, f"ub_deg={D} < 1"
    if diam < 1:
        return False, f"ub_diameter={diam} < 1 for V > 1"

    short = diam // 2
    long  = diam - short
    max_n = (
        1 + D * _max_subtree_size(short - 1, D - 1)
        if diam % 2 == 0 else
        1 + _max_subtree_size(long - 1, D - 1) + (D - 1) * _max_subtree_size(short - 1, D - 1)
    )

    if V > max_n:
        return False, f"V={V} > max possible nodes={max_n} for ub_deg={D}, ub_diameter={diam}"
    return True, ""


def _bfs_farthest(adj, start):
    n = len(adj)
    dist = [-1] * n
    dist[start] = 0
    q = deque([start])
    farthest = start
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                q.append(v)
                if dist[v] > dist[farthest]:
                    farthest = v
    return farthest, dist


def get_diameter(adj) -> int:
    """Compute tree diameter via 2-BFS. O(V)."""
    far1, _ = _bfs_farthest(adj, 0)
    far2, dist = _bfs_farthest(adj, far1)
    return dist[far2]


def generate_undirected_tree_bounded_diameter_degree(
        nb_nodes: int, ub_deg: int, ub_diameter: int):
    """
    Generate a random tree with V nodes, degree <= ub_deg, diameter <= ub_diameter.

    Args:
        nb_nodes:    V >= 1
        ub_deg:      max degree >= 1
        ub_diameter: max diameter >= 0

    Returns:
        Graph, or None if infeasible.
    """
    V, D, diam = nb_nodes, ub_deg, ub_diameter

    ok, msg = check_feasibility(V, D, diam)
    if not ok:
        print(f"[TreeBounded] Infeasible: {msg}")
        return None

    if V == 1:
        return Graph(1)
    if V == 2:
        G = Graph(2); G.AddEdge(0, 1); return G

    short_depth = diam // 2
    long_depth  = diam - short_depth

    edges = []

    # avail: list of [parent_id, child_depth_budget, slots_remaining]
    # Root gets 1 long-arm slot + (D-1) short-arm slots.
    avail = []
    if diam % 2 == 1:
        avail.append([0, long_depth, 1])
        if D > 1:
            avail.append([0, short_depth, D - 1])
    else:
        avail.append([0, short_depth, D])

    next_node = 1
    remaining = V - 1

    while remaining > 0:
        if not avail:
            print(f"[TreeBounded] Error: no available parents, {remaining} nodes unplaced.")
            return None

        total_cap = sum(s * _max_subtree_size(b - 1, D - 1) for (_, b, s) in avail)
        total_slots = sum(s for (_, _, s) in avail)

        # Pick a slot weighted by slot count
        r = random.randint(0, total_slots - 1)
        idx, cumsum = 0, 0
        for i, (_, _, s) in enumerate(avail):
            cumsum += s
            if r < cumsum:
                idx = i; break

        parent_id, child_depth, slots = avail[idx]
        cap_per_child = _max_subtree_size(child_depth - 1, D - 1)
        other_cap = total_cap - slots * cap_per_child

        shortfall = remaining - other_cap
        min_take = max(1, -(-shortfall // max(cap_per_child, 1))) if shortfall > 0 else 1
        max_take = min(slots, remaining)
        min_take = min(min_take, max_take)

        nb_take = min(random.randint(min_take, max_take), remaining)

        for _ in range(nb_take):
            if remaining == 0:
                break
            v = next_node; next_node += 1; remaining -= 1
            edges.append((parent_id, v))
            if child_depth - 1 > 0 and D - 1 > 0:
                avail.append([v, child_depth - 1, D - 1])

        avail[idx][2] -= nb_take
        if avail[idx][2] <= 0:
            avail.pop(idx)

    # Verify diameter (should never fail)
    adj = [[] for _ in range(V)]
    for u, v in edges:
        adj[u].append(v); adj[v].append(u)
    assert get_diameter(adj) <= diam

    # Random vertex permutation
    perm = list(range(V))
    random.shuffle(perm)

    G = Graph(V)
    for u, v in edges:
        G.AddEdge(perm[u], perm[v])
    return G
