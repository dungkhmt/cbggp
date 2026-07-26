"""
Constrained Graph Heuristic Methods (MCMC, Simulated Annealing).
See ALGORITHMS.md for constructive details and MCMC/SA details.

The constructive algorithm (generate_constructive), Tarjan analysis, and
feasibility check are in constructivemethods/constrained_graph.py.
This file adds randomized search on top for cases where diversity matters.
"""

import sys
import os
import random
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph
from constructivemethods.constrained_graph import (
    generate_constructive,
    tarjan_analysis,
    check_feasibility,
    verify_graph,
    _build_graph_from_edges,
)


# MCMC Edge Rewiring

def generate_mcmc(nb_nodes, nb_edges, nb_connected_components,
                  nb_bridges, nb_articulation_points, nb_iterations=1000):
    """
    Generate a random graph satisfying (V, E, C, B, A) via MCMC edge rewiring.

    Starts from a constructive graph, then performs random double-edge swaps
    within the same biconnected component. Each swap is accepted only if
    (C, B, A) are preserved.

    Args:
        nb_iterations: Number of swap attempts (default 1000).

    Returns:
        Graph, or None if constructive base fails.
    """
    V, E, C, B, A = (nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points)

    G = generate_constructive(V, E, C, B, A)
    if G is None:
        return None

    for _ in range(nb_iterations):
        analysis = tarjan_analysis(G)
        bccs = analysis['biconnected_components']

        non_trivial = [bcc for bcc in bccs if len(bcc) >= 4]
        if not non_trivial:
            break

        bcc = random.choice(non_trivial)
        eid1, eid2 = random.sample(list(bcc), 2)
        e1, e2 = G.edges[eid1], G.edges[eid2]
        a, b = e1.fromNode, e1.toNode
        c, d = e2.fromNode, e2.toNode

        if len({a, b, c, d}) < 4:
            continue

        cur_edges = {(min(e.fromNode, e.toNode), max(e.fromNode, e.toNode))
                     for e in G.edges if e.used == 1}
        old_k1 = (min(a, b), max(a, b))
        old_k2 = (min(c, d), max(c, d))

        swap_options = [((a, c), (b, d)), ((a, d), (b, c))]
        random.shuffle(swap_options)

        for (u1, v1), (u2, v2) in swap_options:
            if u1 == v1 or u2 == v2:
                continue
            k1 = (min(u1, v1), max(u1, v1))
            k2 = (min(u2, v2), max(u2, v2))
            if k1 in cur_edges or k2 in cur_edges:
                continue

            new_edges = (cur_edges - {old_k1, old_k2}) | {k1, k2}
            G_new = Graph(V)
            for u, v in new_edges:
                G_new.AddEdge(u, v)

            na = tarjan_analysis(G_new)
            if (na['nb_connected_components'] == C
                    and na['nb_bridges'] == B
                    and na['nb_articulation_points'] == A):
                G = G_new
                break

    return G


# Simulated Annealing

def _random_initial_graph(V, E, C):
    """Build a random starting graph with V nodes, E edges, C components."""
    edge_set = set()
    if V == 0 or E == 0:
        return edge_set

    vertices = list(range(V))
    random.shuffle(vertices)

    comp_sizes = [1] * C
    for _ in range(V - C):
        comp_sizes[random.randint(0, C - 1)] += 1

    offset = 0
    comp_ranges = []
    for size in comp_sizes:
        comp = vertices[offset:offset + size]
        comp_ranges.append(comp)
        for i in range(1, size):
            u, v = comp[i - 1], comp[i]
            edge_set.add((min(u, v), max(u, v)))
        offset += size

    attempts = 0
    while len(edge_set) < E and attempts < E * 20:
        attempts += 1
        comp = random.choice(comp_ranges)
        if len(comp) < 2:
            continue
        u, v = random.sample(comp, 2)
        edge_set.add((min(u, v), max(u, v)))

    return edge_set


def generate_simulated_annealing(nb_nodes, nb_edges, nb_connected_components,
                                 nb_bridges, nb_articulation_points,
                                 max_iterations=10000, T_init=100.0,
                                 alpha=0.995, w_C=10, w_B=3, w_A=1):
    """
    Heuristic optimization via Simulated Annealing.

    Penalty: F(G) = w_C*|C_curr-C| + w_B*|B_curr-B| + w_A*|A_curr-A|
    Goal: F = 0.  Accepts worse states with Boltzmann probability P = exp(-ΔF/T).

    Args:
        max_iterations: SA iteration budget.
        T_init, alpha: Initial temperature and cooling rate.
        w_C, w_B, w_A: Penalty weights.

    Returns:
        Graph if F=0 achieved, None otherwise.
    """
    V, E, C, B, A = (nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points)

    ok, msg = check_feasibility(V, E, C, B, A)
    if not ok:
        print(f"[SA] Infeasible: {msg}")
        return None

    # Try constructive first
    G_try = generate_constructive(V, E, C, B, A)
    if G_try is not None:
        info = tarjan_analysis(G_try)
        if (info['nb_connected_components'] == C
                and info['nb_bridges'] == B
                and info['nb_articulation_points'] == A):
            return G_try

    edge_set = _random_initial_graph(V, E, C)
    all_verts = list(range(V))
    while len(edge_set) < E:
        u, v = random.sample(all_verts, 2)
        key = (min(u, v), max(u, v))
        if key not in edge_set:
            edge_set.add(key)

    def fitness(es):
        Gt = Graph(V)
        for u, v in es:
            Gt.AddEdge(u, v)
        info = tarjan_analysis(Gt)
        return (w_C * abs(info['nb_connected_components'] - C)
                + w_B * abs(info['nb_bridges'] - B)
                + w_A * abs(info['nb_articulation_points'] - A))

    current_fit = fitness(edge_set)
    if current_fit == 0:
        return _build_graph_from_edges(V, edge_set, shuffle=False)

    best_set = set(edge_set)
    best_fit = current_fit
    T = T_init

    for _ in range(max_iterations):
        if current_fit == 0:
            break

        new_set = set(edge_set)
        old_edge = random.choice(list(new_set))
        new_set.discard(old_edge)

        added = False
        for _ in range(50):
            u = random.randint(0, V - 1)
            v = random.randint(0, V - 1)
            if u == v:
                continue
            key = (min(u, v), max(u, v))
            if key not in new_set:
                new_set.add(key)
                added = True
                break

        if not added:
            T *= alpha
            continue

        new_fit = fitness(new_set)
        delta = new_fit - current_fit

        if delta <= 0 or random.random() < math.exp(-delta / max(T, 1e-10)):
            edge_set = new_set
            current_fit = new_fit
            if current_fit < best_fit:
                best_fit = current_fit
                best_set = set(edge_set)

        T *= alpha

    if best_fit == 0:
        return _build_graph_from_edges(V, best_set, shuffle=False)

    print(f"[SA] Did not converge after {max_iterations} iterations. "
          f"Best fitness = {best_fit}")
    return None


# Unified dispatcher (used by api.py gen_undirected_graph)

def generate_graph(nb_nodes, nb_edges, nb_connected_components,
                   nb_bridges, nb_articulation_points,
                   method='constructive', **kwargs):
    """
    Generate a graph satisfying (V, E, C, B, A) using the chosen method.

    method:
        'constructive' — 4-phase constructive algorithm (default, exact)
        'mcmc'         — MCMC edge rewiring (diverse random graphs)
                         kwargs: nb_iterations=1000
        'sa'           — Simulated Annealing (heuristic fallback)
                         kwargs: max_iterations=10000, T_init=100.0, alpha=0.995

    Returns:
        Graph, or None.
    """
    if method == 'constructive':
        return generate_constructive(nb_nodes, nb_edges, nb_connected_components,
                                     nb_bridges, nb_articulation_points)
    elif method == 'mcmc':
        return generate_mcmc(nb_nodes, nb_edges, nb_connected_components,
                             nb_bridges, nb_articulation_points,
                             kwargs.get('nb_iterations', 1000))
    elif method == 'sa':
        return generate_simulated_annealing(
            nb_nodes, nb_edges, nb_connected_components,
            nb_bridges, nb_articulation_points,
            max_iterations=kwargs.get('max_iterations', 10000),
            T_init=kwargs.get('T_init', 100.0),
            alpha=kwargs.get('alpha', 0.995),
            w_C=kwargs.get('w_C', 10),
            w_B=kwargs.get('w_B', 3),
            w_A=kwargs.get('w_A', 1),
        )
    else:
        raise ValueError(f"Invalid method: {method!r}. Choose 'constructive', 'mcmc', or 'sa'.")
