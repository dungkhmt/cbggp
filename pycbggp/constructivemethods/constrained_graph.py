"""
Constrained Undirected Graph Generator (V, E, C, B, A). Author: Nguyen Ngoc Tuan Anh
See ALGORITHMS.md for algorithm details (Tarjan analysis + 4-phase constructive).
"""

import sys
import os
import random

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph


# Tarjan-Hopcroft structural analysis  O(V+E)

def tarjan_analysis(G):
    """
    Compute structural properties of graph G via iterative Tarjan DFS.

    Returns:
        dict with keys: 'nb_connected_components', 'nb_bridges',
        'nb_articulation_points', 'bridges', 'articulation_points',
        'biconnected_components'
    """
    n = G.n
    adj = [[] for _ in range(n)]
    for e in G.edges:
        if e.used == 1:
            adj[e.fromNode].append((e.toNode, e.id))
            adj[e.toNode].append((e.fromNode, e.id))

    dfn = [0] * n
    low = [0] * n
    child_count = [0] * n
    timer = [0]
    bridges = set()
    aps = set()
    bccs = []
    nb_cc = 0
    edge_stack = []

    for start in range(n):
        if dfn[start] != 0:
            continue
        nb_cc += 1
        timer[0] += 1
        dfn[start] = low[start] = timer[0]
        call_stack = [(start, -1, 0)]

        while call_stack:
            u, p_eid, pos = call_stack[-1]

            if pos < len(adj[u]):
                call_stack[-1] = (u, p_eid, pos + 1)
                v, eid = adj[u][pos]

                if eid == p_eid:
                    continue

                if dfn[v] == 0:
                    child_count[u] += 1
                    edge_stack.append(eid)
                    timer[0] += 1
                    dfn[v] = low[v] = timer[0]
                    call_stack.append((v, eid, 0))
                else:
                    if dfn[v] < dfn[u]:
                        edge_stack.append(eid)
                        low[u] = min(low[u], dfn[v])
            else:
                call_stack.pop()

                if call_stack:
                    parent_u = call_stack[-1][0]
                    low[parent_u] = min(low[parent_u], low[u])

                    if low[u] > dfn[parent_u]:
                        bridges.add(p_eid)

                    if low[u] >= dfn[parent_u]:
                        component = set()
                        while edge_stack and edge_stack[-1] != p_eid:
                            component.add(edge_stack.pop())
                        if edge_stack:
                            component.add(edge_stack.pop())
                        if component:
                            bccs.append(component)

                        if len(call_stack) >= 2:
                            aps.add(parent_u)
                        else:
                            if child_count[parent_u] >= 2:
                                aps.add(parent_u)
                else:
                    if edge_stack:
                        bccs.append(set(edge_stack))
                        edge_stack.clear()

    return {
        'nb_connected_components': nb_cc,
        'nb_bridges': len(bridges),
        'nb_articulation_points': len(aps),
        'bridges': bridges,
        'articulation_points': aps,
        'biconnected_components': bccs,
    }


# Feasibility validation

def check_feasibility(nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points):
    """
    Check whether (V, E, C, B, A) is feasible for an undirected graph.

    Returns:
        (is_feasible: bool, error_message: str)
    """
    V, E, C, B, A = (nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points)

    if V < 0 or E < 0 or C < 0 or B < 0 or A < 0:
        return False, "All parameters must be >= 0"

    if V == 0:
        if E != 0 or C != 0 or B != 0 or A != 0:
            return False, "V=0 requires E=C=B=A=0"
        return True, ""

    if C == 0:
        return False, "C must be >= 1 when V >= 1"
    if V < C:
        return False, f"V={V} < C={C}: not enough vertices for components"

    V_main = V - C + 1
    E_min = V - C
    E_max = V_main * (V_main - 1) // 2

    if E < E_min:
        return False, f"E={E} < E_min={E_min} (need at least V-C edges)"
    if E > E_max:
        return False, f"E={E} > E_max={E_max} (exceeds complete graph capacity)"

    if B > E:
        return False, f"B={B} > E={E}: bridges cannot exceed edges"
    if B > V - C:
        return False, f"B={B} > V-C={V - C}: exceeds bridge upper bound"
    if E == V - C and B != E:
        return False, f"Forest (E=V-C={E}) requires B=E, but B={B}"

    if V_main <= 1:
        if A != 0:
            return False, f"V_main={V_main} requires A=0"
    elif V_main == 2:
        if A != 0:
            return False, f"V_main=2 requires A=0"
        if E > 0 and B != 1:
            return False, f"V_main=2, E>0 requires B=1"
    else:
        if A > V_main - 2:
            return False, f"A={A} > V_main-2={V_main - 2}: exceeds AP upper bound"
        if A == 0 and B > 0:
            return False, "A=0 with B>0 is infeasible when V_main >= 3"
        if B >= 1 and A < 1:
            return False, "B>=1 requires A>=1 when V_main >= 3"
        if 0 < E - B < 3:
            return False, f"Non-bridge edges E-B={E - B} must be 0 or >= 3"

        if B == 0 and A >= 1:
            if V_main < 2 * A + 3:
                return False, f"V_main={V_main} < 2A+3={2*A+3}: insufficient for B=0, A={A}"
            if E < A + V_main:
                return False, f"E={E} < A+V_main={A+V_main}: insufficient for B=0, A={A}"
        elif B == 1 and A >= 1:
            if V_main < 2 + 2 * A:
                return False, f"V_main={V_main} < 2+2A={2+2*A}: insufficient for B=1, A={A}"
        elif B >= 2 and A > B - 1:
            A_extra = A - (B - 1)
            min_v = B + 1 + 2 * A_extra
            if V_main < min_v:
                return False, f"V_main={V_main} < {min_v}: insufficient for B={B}, A={A}"

    return True, ""


# Helpers for constructive algorithm

def _chain_triangles(start_vertex, count, vertex_idx, edge_set, blocks, ap_set):
    """
    Append 'count' chained triangles starting from start_vertex.
    Each triangle adds 1 articulation point and 0 bridges.
    Returns (new vertex_idx, last vertex of the chain).
    """
    current = start_vertex
    for _ in range(count):
        a, b = vertex_idx, vertex_idx + 1
        edge_set.add((min(current, a), max(current, a)))
        edge_set.add((min(a, b), max(a, b)))
        edge_set.add((min(b, current), max(b, current)))
        blocks.append({current, a, b})
        ap_set.add(current)
        current = b
        vertex_idx += 2
    return vertex_idx, current


def _build_graph_from_edges(V, edge_set, shuffle=True):
    """Build a Graph from an edge set, optionally with a random vertex permutation."""
    G = Graph(V)
    if V == 0:
        return G
    perm = list(range(V))
    edges = list(edge_set)
    if shuffle:
        random.shuffle(perm)
        random.shuffle(edges)
    for u, v in edges:
        G.AddEdge(perm[u], perm[v])
    return G


# 4-Phase Constructive Algorithm

def generate_constructive(nb_nodes, nb_edges, nb_connected_components,
                          nb_bridges, nb_articulation_points):
    """
    Generate an undirected graph satisfying exactly (V, E, C, B, A).

    4 phases:
      Phase 1: C-1 isolated vertices (extra components).
      Phase 2: Bridge-AP backbone assembly (caterpillar/path/triangle chains).
      Phase 3: Vertex padding via edge subdivision inside blocks.
      Phase 4: Edge densification within biconnected blocks.

    Args:
        nb_nodes, nb_edges, nb_connected_components, nb_bridges, nb_articulation_points

    Returns:
        Graph, or None if infeasible or construction fails.
    """
    V, E, C, B, A = (nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points)

    ok, msg = check_feasibility(V, E, C, B, A)
    if not ok:
        print(f"[Constructive] Infeasible: {msg}")
        return None

    if V == 0:
        return Graph(0)

    V_main = V - C + 1

    edge_set   = set()
    bridge_set = set()
    blocks     = []
    ap_set     = set()

    def add_edge(u, v, is_bridge=False):
        key = (min(u, v), max(u, v))
        if key in edge_set:
            return False
        edge_set.add(key)
        if is_bridge:
            bridge_set.add(key)
        return True

    def has_edge(u, v):
        return (min(u, v), max(u, v)) in edge_set

    def remove_edge(u, v):
        key = (min(u, v), max(u, v))
        edge_set.discard(key)
        bridge_set.discard(key)

    vertex_idx = 0

    # --- Special cases ---
    if V_main <= 1:
        return _build_graph_from_edges(V, edge_set)
    if V_main == 2:
        add_edge(0, 1, is_bridge=True)
        return _build_graph_from_edges(V, edge_set)

    # --- Phase 2a: Bridge backbone ---
    if B == 1:
        add_edge(0, 1, is_bridge=True)
        vertex_idx = 2
    elif B >= 2:
        if A <= B - 1:
            if A == 1:
                for i in range(B):
                    add_edge(0, i + 1, is_bridge=True)
                ap_set.add(0)
                vertex_idx = B + 1
            else:
                for i in range(A - 1):
                    add_edge(i, i + 1, is_bridge=True)
                vertex_idx = A
                for i in range(B - A + 1):
                    add_edge(0, vertex_idx, is_bridge=True)
                    vertex_idx += 1
                for i in range(A - 1):
                    ap_set.add(i)
        else:
            for i in range(B):
                add_edge(i, i + 1, is_bridge=True)
            vertex_idx = B + 1
            for i in range(1, B):
                ap_set.add(i)

    A_curr  = len(ap_set)
    A_extra = A - A_curr

    # --- Phase 2b: Cycle blocks for additional APs ---
    if B == 0 and A == 0:
        for i in range(V_main):
            add_edge(i, (i + 1) % V_main)
        blocks.append(set(range(V_main)))
        vertex_idx = V_main

    elif B == 0 and A >= 1:
        add_edge(0, 1); add_edge(1, 2); add_edge(2, 0)
        blocks.append({0, 1, 2})
        vertex_idx = 3
        vertex_idx, _ = _chain_triangles(2, A, vertex_idx, edge_set, blocks, ap_set)
        A_extra = 0

    elif A_extra > 0:
        if B == 1:
            chain_sources = [0, 1]
        elif B >= 2 and A_curr <= B - 1:
            chain_sources = [v for v in range(vertex_idx) if v not in ap_set]
        else:
            chain_sources = [0, B]

        remaining = A_extra
        src_idx = 0
        last_chain_end = None
        while remaining > 0:
            if src_idx < len(chain_sources):
                start = chain_sources[src_idx]; src_idx += 1
            elif last_chain_end is not None:
                start = last_chain_end
            else:
                break
            count = 1 if src_idx < len(chain_sources) and remaining > 1 else remaining
            vertex_idx, last_chain_end = _chain_triangles(
                start, count, vertex_idx, edge_set, blocks, ap_set)
            remaining -= count

    # --- Phase 2c: Ensure at least one non-bridge block exists ---
    V_rem = V_main - vertex_idx
    E_rem = E - len(edge_set)
    if (V_rem > 0 or E_rem > 0) and len(blocks) == 0:
        if not ap_set:
            return None
        ap_v = next(iter(ap_set))
        a, b = vertex_idx, vertex_idx + 1
        add_edge(ap_v, a); add_edge(a, b); add_edge(b, ap_v)
        blocks.append({ap_v, a, b})
        vertex_idx += 2

    # --- Phase 3: Vertex padding via edge subdivision ---
    V_rem = V_main - vertex_idx
    if V_rem < 0:
        return None
    for _ in range(V_rem):
        subdivided = False
        for bi in range(len(blocks)):
            if subdivided:
                break
            block = blocks[bi]
            verts = list(block)
            for i in range(len(verts)):
                if subdivided:
                    break
                for j in range(i + 1, len(verts)):
                    u, v = verts[i], verts[j]
                    key = (min(u, v), max(u, v))
                    if key in edge_set and key not in bridge_set:
                        x = vertex_idx; vertex_idx += 1
                        remove_edge(u, v)
                        add_edge(u, x); add_edge(x, v)
                        block.add(x)
                        subdivided = True
                        break
        if not subdivided:
            return None

    # --- Phase 4: Edge densification within blocks ---
    E_rem = E - len(edge_set)
    if E_rem < 0:
        return None
    for bi in range(len(blocks)):
        if E_rem <= 0:
            break
        verts = list(blocks[bi])
        for i in range(len(verts)):
            if E_rem <= 0:
                break
            for j in range(i + 1, len(verts)):
                if E_rem <= 0:
                    break
                u, v = verts[i], verts[j]
                if not has_edge(u, v):
                    add_edge(u, v)
                    E_rem -= 1

    if E_rem > 0:
        return None

    return _build_graph_from_edges(V, edge_set)


# Verification utility

def verify_graph(G, nb_nodes, nb_edges, nb_connected_components,
                 nb_bridges, nb_articulation_points):
    """
    Verify whether G satisfies (V, E, C, B, A) using Tarjan's algorithm.

    Returns:
        (is_valid: bool, details: dict)
    """
    V, E, C, B, A = (nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points)
    actual_V = G.n
    actual_E = sum(1 for e in G.edges if e.used == 1)
    info = tarjan_analysis(G)
    actual_C = info['nb_connected_components']
    actual_B = info['nb_bridges']
    actual_A = info['nb_articulation_points']

    is_valid = (actual_V == V and actual_E == E and actual_C == C
                and actual_B == B and actual_A == A)
    return is_valid, {
        'expected': {'V': V, 'E': E, 'C': C, 'B': B, 'A': A},
        'actual':   {'V': actual_V, 'E': actual_E, 'C': actual_C,
                     'B': actual_B, 'A': actual_A},
        'is_valid': is_valid,
    }
