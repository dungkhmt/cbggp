"""
Constrained Structural Graph Generator (V, E, C, B, A)
=======================================================
Generates an undirected graph G = (V, E) satisfying exactly 5 structural constraints:
  V = nb_nodes:                  Number of vertices
  E = nb_edges:                  Number of edges
  C = nb_connected_components:   Number of connected components
  B = nb_bridges:                Number of bridge edges (cut edges)
  A = nb_articulation_points:    Number of articulation points (cut vertices)

Implemented algorithms:
  1. Tarjan-Hopcroft structural verification (O(V+E))
  2. Feasibility Validation
  3. Classical Constructive Algorithm (4-phase)
  4. Randomized MCMC Edge Rewiring
  5. Simulated Annealing heuristic optimization

References:
  - Hopcroft & Tarjan (1973). Efficient algorithms for graph manipulation.
  - Block-Cut Tree Decomposition.
"""

import sys
import os
import random
import math

# Add parent directory to sys.path for CBSGG import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph


# ============================================================
# SECTION 1: TARJAN-HOPCROFT STRUCTURAL VERIFICATION
# ============================================================

def tarjan_analysis(G):
    """
    Analyze graph G using iterative Tarjan's DFS to compute:
      - Number of connected components
      - Number of bridges and the set of bridge edges
      - Number of articulation points and the set of cut vertices
      - Biconnected components

    Args:
        G: Graph object from CBSGG.py

    Returns:
        dict with 'nb_connected_components', 'nb_bridges',
        'nb_articulation_points', 'bridges', 'articulation_points',
        'biconnected_components'

    Time complexity: O(V + E)
    """
    n = G.n
    # Build adjacency list from Graph (only used edges)
    adj = [[] for _ in range(n)]
    for e in G.edges:
        if e.used == 1:
            adj[e.fromNode].append((e.toNode, e.id))
            adj[e.toNode].append((e.fromNode, e.id))

    dfn = [0] * n          # Discovery time
    low = [0] * n          # Lowest reachable ancestor time
    child_count = [0] * n  # Number of children in DFS tree
    timer = [0]
    bridges = set()        # Set of bridge edge IDs
    aps = set()            # Set of articulation point vertices
    bccs = []              # List of biconnected components (each = set of edge IDs)
    nb_cc = 0
    edge_stack = []        # Edge stack for BCC extraction

    for start in range(n):
        if dfn[start] != 0:
            continue
        nb_cc += 1
        timer[0] += 1
        dfn[start] = low[start] = timer[0]

        # Iterative DFS: each frame = (u, parent_edge_id, adj_position)
        call_stack = [(start, -1, 0)]

        while call_stack:
            u, p_eid, pos = call_stack[-1]

            if pos < len(adj[u]):
                # Visit next neighbor
                call_stack[-1] = (u, p_eid, pos + 1)
                v, eid = adj[u][pos]

                if eid == p_eid:
                    continue  # Skip parent edge

                if dfn[v] == 0:
                    # Tree edge: vertex v not yet visited
                    child_count[u] += 1
                    edge_stack.append(eid)
                    timer[0] += 1
                    dfn[v] = low[v] = timer[0]
                    call_stack.append((v, eid, 0))
                else:
                    # Back edge: v already visited before u
                    if dfn[v] < dfn[u]:
                        edge_stack.append(eid)
                        low[u] = min(low[u], dfn[v])
            else:
                # All neighbors processed -> "return" from DFS
                call_stack.pop()

                if call_stack:
                    parent_u = call_stack[-1][0]
                    low[parent_u] = min(low[parent_u], low[u])

                    # Bridge criterion: low[u] > dfn[parent_u]
                    if low[u] > dfn[parent_u]:
                        bridges.add(p_eid)

                    # Extract BCC when low[u] >= dfn[parent_u]
                    if low[u] >= dfn[parent_u]:
                        component = set()
                        while edge_stack and edge_stack[-1] != p_eid:
                            component.add(edge_stack.pop())
                        if edge_stack:
                            component.add(edge_stack.pop())
                        if component:
                            bccs.append(component)

                        # Articulation point criterion
                        if len(call_stack) >= 2:
                            # parent_u is not root -> articulation point
                            aps.add(parent_u)
                        else:
                            # parent_u is root -> AP only if >= 2 children
                            if child_count[parent_u] >= 2:
                                aps.add(parent_u)
                else:
                    # u is root, remaining edges on stack = last BCC
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


# ============================================================
# SECTION 2: FEASIBILITY VALIDATION
# ============================================================

def check_feasibility(nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points):
    """
    Check whether parameter tuple (V, E, C, B, A) is feasible.

    Returns:
        (is_feasible: bool, error_message: str)
    """
    V, E, C, B, A = (nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points)

    # --- Non-negative constraints ---
    if V < 0 or E < 0 or C < 0 or B < 0 or A < 0:
        return False, "All parameters must be >= 0"

    # --- Empty graph ---
    if V == 0:
        if E != 0 or C != 0 or B != 0 or A != 0:
            return False, "V=0 requires E=C=B=A=0"
        return True, ""

    if C == 0:
        return False, "C must be >= 1 when V >= 1"

    # --- Basic constraints ---
    if V < C:
        return False, f"V={V} < C={C}: not enough vertices for components"

    V_main = V - C + 1  # Main component size (resource consolidation)
    E_min = V - C       # Edge lower bound (forest)
    E_max = V_main * (V_main - 1) // 2  # Edge upper bound (complete graph)

    if E < E_min:
        return False, f"E={E} < E_min={E_min} (need at least V-C edges)"
    if E > E_max:
        return False, f"E={E} > E_max={E_max} (exceeds complete graph capacity)"

    # --- Bridge constraints ---
    if B > E:
        return False, f"B={B} > E={E}: bridges cannot exceed edges"
    if B > V - C:
        return False, f"B={B} > V-C={V - C}: exceeds bridge upper bound"
    if E == V - C and B != E:
        return False, f"Forest (E=V-C={E}) requires B=E, but B={B}"

    # --- Articulation point constraints ---
    if V_main <= 1:
        if A != 0:
            return False, f"V_main={V_main} requires A=0"
    elif V_main == 2:
        if A != 0:
            return False, f"V_main=2 requires A=0"
        if E > 0 and B != 1:
            return False, f"V_main=2, E>0 requires B=1"
    else:
        # V_main >= 3
        if A > V_main - 2:
            return False, f"A={A} > V_main-2={V_main - 2}: exceeds AP upper bound"
        if A == 0 and B > 0:
            return False, "A=0 with B>0 is infeasible when V_main >= 3"
        if B >= 1 and A < 1:
            return False, "B>=1 requires A>=1 when V_main >= 3"

        # --- Non-bridge edge constraints ---
        if 0 < E - B < 3:
            return False, f"Non-bridge edges E-B={E - B} must be 0 or >= 3"

        # --- Minimum vertex constraints for structure ---
        if B == 0 and A >= 1:
            # Need A+1 cycle blocks: minimum 2A+3 vertices
            if V_main < 2 * A + 3:
                return False, (f"V_main={V_main} < 2A+3={2 * A + 3}: "
                               f"not enough vertices for B=0, A={A}")
            # Edge lower bound: each block minimum = cycle
            if E < A + V_main:
                return False, (f"E={E} < A+V_main={A + V_main}: "
                               f"not enough edges for B=0, A={A}")
        elif B == 1 and A >= 1:
            if V_main < 2 + 2 * A:
                return False, (f"V_main={V_main} < 2+2A={2 + 2 * A}: "
                               f"not enough vertices for B=1, A={A}")
        elif B >= 2 and A > B - 1:
            A_extra = A - (B - 1)
            min_v = B + 1 + 2 * A_extra
            if V_main < min_v:
                return False, (f"V_main={V_main} < {min_v}: "
                               f"not enough vertices for B={B}, A={A}")

    return True, ""


# ============================================================
# SECTION 3: CLASSICAL CONSTRUCTIVE ALGORITHM (4-PHASE)
# ============================================================

def _chain_triangles(start_vertex, count, vertex_idx, edge_set, blocks, ap_set):
    """
    Chain 'count' triangles starting from start_vertex.
    Each triangle creates 1 new articulation point, 0 new bridges.

    Returns: (new vertex_idx, last vertex of the chain)
    """
    current = start_vertex
    for _ in range(count):
        a = vertex_idx
        b = vertex_idx + 1
        edge_set.add((min(current, a), max(current, a)))
        edge_set.add((min(a, b), max(a, b)))
        edge_set.add((min(b, current), max(b, current)))
        blocks.append({current, a, b})
        ap_set.add(current)
        current = b
        vertex_idx += 2
    return vertex_idx, current


def _build_graph_from_edges(V, edge_set, shuffle=True):
    """Build a Graph object from an edge set, with shuffled vertex IDs."""
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


def generate_constructive(nb_nodes, nb_edges, nb_connected_components,
                          nb_bridges, nb_articulation_points):
    """
    Classical Constructive Algorithm — generates an undirected graph satisfying
    exactly (V, E, C, B, A) in near-linear time.

    4 phases:
      Phase 1: Component initialization (C-1 isolated vertices)
      Phase 2: Bridge-Articulation Point backbone assembly
      Phase 3: Vertex padding (Edge subdivision)
      Phase 4: Edge densification within blocks

    Returns: Graph or None if infeasible
    """
    V, E, C, B, A = (nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points)

    ok, msg = check_feasibility(V, E, C, B, A)
    if not ok:
        print(f"[Constructive] Infeasible: {msg}")
        return None

    if V == 0:
        return Graph(0)

    V_main = V - C + 1  # Main component vertex count

    # Build graph using in-memory edge set
    edge_set = set()   # {(min(u,v), max(u,v))}
    bridge_set = set()  # Set of bridge edges
    blocks = []         # List of non-bridge blocks (sets of vertices)
    ap_set = set()      # Set of articulation points

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

    vertex_idx = 0  # Next available vertex index

    # === Special cases ===
    if V_main <= 1:
        return _build_graph_from_edges(V, edge_set)

    if V_main == 2:
        add_edge(0, 1, is_bridge=True)
        return _build_graph_from_edges(V, edge_set)

    # === V_main >= 3: Start 4-phase construction ===

    # ──────────────────────────────────────────────
    # PHASE 2a: Bridge Backbone
    # ──────────────────────────────────────────────
    if B == 0:
        pass  # No bridges needed, handled in Phase 2b
    elif B == 1:
        # Single bridge: v0 - v1
        add_edge(0, 1, is_bridge=True)
        vertex_idx = 2
    elif B >= 2:
        if A <= B - 1:
            # Caterpillar: spine + pendant edges for exactly A APs from bridges
            if A == 1:
                # Star: v0 center, B leaves
                for i in range(B):
                    add_edge(0, i + 1, is_bridge=True)
                ap_set.add(0)
                vertex_idx = B + 1
            else:
                # Path + pendants: A-1 path edges + (B-A+1) pendant edges from v0
                for i in range(A - 1):
                    add_edge(i, i + 1, is_bridge=True)
                vertex_idx = A
                for i in range(B - A + 1):
                    add_edge(0, vertex_idx, is_bridge=True)
                    vertex_idx += 1
                # APs: v0 (hub), v1..v_{A-2} (internal path vertices)
                for i in range(A - 1):
                    ap_set.add(i)
        else:
            # A > B-1: use pure path, supplement APs with cycle blocks
            for i in range(B):
                add_edge(i, i + 1, is_bridge=True)
            vertex_idx = B + 1
            # APs from path: v1..v_{B-1}
            for i in range(1, B):
                ap_set.add(i)

    A_curr = len(ap_set)
    A_extra = A - A_curr

    # ──────────────────────────────────────────────
    # PHASE 2b: Cycle blocks for additional APs
    # ──────────────────────────────────────────────
    if B == 0 and A == 0:
        # Biconnected graph: build Hamiltonian cycle
        for i in range(V_main):
            add_edge(i, (i + 1) % V_main)
        blocks.append(set(range(V_main)))
        vertex_idx = V_main

    elif B == 0 and A >= 1:
        # Chain of A+1 triangles
        add_edge(0, 1)
        add_edge(1, 2)
        add_edge(2, 0)
        blocks.append({0, 1, 2})
        vertex_idx = 3
        vertex_idx, _ = _chain_triangles(2, A, vertex_idx,
                                         edge_set, blocks, ap_set)
        A_extra = 0

    elif A_extra > 0:
        # Chain triangles from bridge tree leaves
        if B == 1:
            chain_sources = [0, 1]
        elif B >= 2 and A_curr <= B - 1:
            # Caterpillar leaves: find non-AP vertices
            chain_sources = [v for v in range(vertex_idx) if v not in ap_set]
        else:
            # Path leaves: v0 and vB
            chain_sources = [0, B]

        remaining = A_extra
        src_idx = 0
        last_chain_end = None

        while remaining > 0:
            if src_idx < len(chain_sources):
                start = chain_sources[src_idx]
                src_idx += 1
            elif last_chain_end is not None:
                start = last_chain_end
            else:
                break

            # Distribute: take 1 if more sources remain, otherwise take all
            count = 1 if src_idx < len(chain_sources) and remaining > 1 else remaining
            vertex_idx, last_chain_end = _chain_triangles(
                start, count, vertex_idx, edge_set, blocks, ap_set)
            remaining -= count

        A_extra = 0

    # ──────────────────────────────────────────────
    # PHASE 2c: Ensure non-bridge blocks exist
    # ──────────────────────────────────────────────
    V_rem = V_main - vertex_idx
    E_rem = E - len(edge_set)

    if (V_rem > 0 or E_rem > 0) and len(blocks) == 0:
        if ap_set:
            ap_v = next(iter(ap_set))
            a, b = vertex_idx, vertex_idx + 1
            add_edge(ap_v, a)
            add_edge(a, b)
            add_edge(b, ap_v)
            blocks.append({ap_v, a, b})
            vertex_idx += 2
        else:
            return None

    # ──────────────────────────────────────────────
    # PHASE 3: Vertex Padding (Edge Subdivision)
    # ──────────────────────────────────────────────
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
                        x = vertex_idx
                        vertex_idx += 1
                        remove_edge(u, v)
                        add_edge(u, x)
                        add_edge(x, v)
                        block.add(x)
                        subdivided = True
                        break
        if not subdivided:
            return None

    # ──────────────────────────────────────────────
    # PHASE 4: Edge Densification within blocks
    # ──────────────────────────────────────────────
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
        return None  # Not enough block capacity

    return _build_graph_from_edges(V, edge_set)


# ============================================================
# SECTION 4: RANDOMIZED MCMC EDGE REWIRING
# ============================================================

def generate_mcmc(nb_nodes, nb_edges, nb_connected_components,
                  nb_bridges, nb_articulation_points, nb_iterations=1000):
    """
    Generate diverse random graphs satisfying (V, E, C, B, A) using
    Markov Chain Monte Carlo (MCMC) with property-preserving edge swaps.

    Starts from a constructive graph, performs double edge swaps within
    the same biconnected component, verified by Tarjan after each step.

    Returns: Graph or None
    """
    V, E, C, B, A = (nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points)

    G = generate_constructive(V, E, C, B, A)
    if G is None:
        return None

    for _ in range(nb_iterations):
        analysis = tarjan_analysis(G)
        bridge_eids = analysis['bridges']
        bccs = analysis['biconnected_components']

        # Find non-trivial BCCs (>= 4 edges for swap)
        non_trivial = [bcc for bcc in bccs if len(bcc) >= 4]
        if not non_trivial:
            break

        bcc = random.choice(non_trivial)
        bcc_eids = list(bcc)
        eid1, eid2 = random.sample(bcc_eids, 2)
        e1 = G.edges[eid1]
        e2 = G.edges[eid2]
        a, b = e1.fromNode, e1.toNode
        c, d = e2.fromNode, e2.toNode

        if len({a, b, c, d}) < 4:
            continue

        # Build current edge set
        cur_edges = set()
        for e in G.edges:
            if e.used == 1:
                cur_edges.add((min(e.fromNode, e.toNode),
                               max(e.fromNode, e.toNode)))

        old_k1 = (min(a, b), max(a, b))
        old_k2 = (min(c, d), max(c, d))

        # Try two cross-swap options
        swap_options = [((a, c), (b, d)), ((a, d), (b, c))]
        random.shuffle(swap_options)

        for (ne1, ne2) in swap_options:
            u1, v1 = ne1
            u2, v2 = ne2
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
            if (na['nb_connected_components'] == C and
                    na['nb_bridges'] == B and
                    na['nb_articulation_points'] == A):
                G = G_new
                break

    return G


# ============================================================
# SECTION 5: HEURISTIC OPTIMIZATION — SIMULATED ANNEALING
# ============================================================

def _random_initial_graph(V, E, C):
    """
    Create a random initial graph with V vertices, E edges, C components.
    Used as starting point for Simulated Annealing.
    """
    edge_set = set()
    if V == 0 or E == 0:
        return edge_set

    vertices = list(range(V))
    random.shuffle(vertices)

    # Distribute vertices among C components
    comp_sizes = [1] * C
    for _ in range(V - C):
        comp_sizes[random.randint(0, C - 1)] += 1

    # Build spanning tree for each component
    offset = 0
    comp_ranges = []
    for size in comp_sizes:
        comp = vertices[offset:offset + size]
        comp_ranges.append(comp)
        for i in range(1, size):
            u, v = comp[i - 1], comp[i]
            edge_set.add((min(u, v), max(u, v)))
        offset += size

    # Add random intra-component edges
    attempts = 0
    max_attempts = E * 20
    while len(edge_set) < E and attempts < max_attempts:
        attempts += 1
        comp = random.choice(comp_ranges)
        if len(comp) < 2:
            continue
        u, v = random.sample(comp, 2)
        key = (min(u, v), max(u, v))
        if key not in edge_set:
            edge_set.add(key)

    return edge_set


def generate_simulated_annealing(nb_nodes, nb_edges, nb_connected_components,
                                 nb_bridges, nb_articulation_points,
                                 max_iterations=10000, T_init=100.0,
                                 alpha=0.995, w_C=10, w_B=3, w_A=1):
    """
    Heuristic optimization using Simulated Annealing (SA).

    Penalty Fitness Function:
        F(G) = w_C * |C_curr - C| + w_B * |B_curr - B| + w_A * |A_curr - A|
    Goal: F(G) = 0.

    Mutation Operators:
        - Rewiring: remove 1 random edge, add 1 new edge
        - Bridging: break cycle to increase B (remove cycle edge,
                    reconnect through another vertex)

    Accepts worse states via Boltzmann probability:
        P = exp(-ΔF / T)

    Returns: Graph if F=0, None if not converged
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
        a = tarjan_analysis(G_try)
        if (a['nb_connected_components'] == C and
                a['nb_bridges'] == B and
                a['nb_articulation_points'] == A):
            return G_try

    # Initialize random graph
    edge_set = _random_initial_graph(V, E, C)
    if len(edge_set) != E:
        # Not enough edges, try adding cross-component edges
        all_verts = list(range(V))
        while len(edge_set) < E:
            u, v = random.sample(all_verts, 2)
            key = (min(u, v), max(u, v))
            if key not in edge_set:
                edge_set.add(key)

    def fitness(es):
        """Compute penalty fitness F(G)."""
        Gt = Graph(V)
        for u, v in es:
            Gt.AddEdge(u, v)
        info = tarjan_analysis(Gt)
        return (w_C * abs(info['nb_connected_components'] - C) +
                w_B * abs(info['nb_bridges'] - B) +
                w_A * abs(info['nb_articulation_points'] - A))

    current_fit = fitness(edge_set)
    if current_fit == 0:
        return _build_graph_from_edges(V, edge_set, shuffle=False)

    best_set = set(edge_set)
    best_fit = current_fit
    T = T_init

    for iteration in range(max_iterations):
        if current_fit == 0:
            break

        # --- Mutation operator: Rewiring ---
        new_set = set(edge_set)
        edges_list = list(new_set)

        # Select edge to remove
        old_edge = random.choice(edges_list)
        new_set.discard(old_edge)

        # Select new edge to add
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

        # Metropolis criterion
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


# ============================================================
# SECTION 6: UTILITY FUNCTIONS AND TESTING
# ============================================================

def verify_graph(G, nb_nodes, nb_edges, nb_connected_components,
                 nb_bridges, nb_articulation_points):
    """
    Verify whether graph G satisfies exactly (V, E, C, B, A).
    Uses Tarjan's algorithm.

    Returns: (is_valid: bool, details: dict)
    """
    V, E, C, B, A = (nb_nodes, nb_edges, nb_connected_components,
                      nb_bridges, nb_articulation_points)

    actual_V = G.n
    actual_E = sum(1 for e in G.edges if e.used == 1)

    analysis = tarjan_analysis(G)
    actual_C = analysis['nb_connected_components']
    actual_B = analysis['nb_bridges']
    actual_A = analysis['nb_articulation_points']

    is_valid = (actual_V == V and actual_E == E and actual_C == C and
                actual_B == B and actual_A == A)

    details = {
        'expected': {'V': V, 'E': E, 'C': C, 'B': B, 'A': A},
        'actual': {'V': actual_V, 'E': actual_E, 'C': actual_C,
                   'B': actual_B, 'A': actual_A},
        'is_valid': is_valid,
    }
    return is_valid, details


def generate_graph(nb_nodes, nb_edges, nb_connected_components,
                   nb_bridges, nb_articulation_points,
                   method='constructive', **kwargs):
    """
    Unified function: generate a graph using the selected method.

    method:
        'constructive' — Classical Constructive Algorithm (default)
        'mcmc'         — MCMC Edge Rewiring (random diversity)
        'sa'           — Simulated Annealing (heuristic)

    Returns: Graph or None
    """
    if method == 'constructive':
        return generate_constructive(nb_nodes, nb_edges,
                                     nb_connected_components,
                                     nb_bridges, nb_articulation_points)
    elif method == 'mcmc':
        iters = kwargs.get('nb_iterations', 1000)
        return generate_mcmc(nb_nodes, nb_edges, nb_connected_components,
                             nb_bridges, nb_articulation_points, iters)
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
        raise ValueError(f"Invalid method: {method}")


# ============================================================
# DEMO AND TESTING
# ============================================================

if __name__ == "__main__":
    test_cases = [
        # (V, E, C, B, A,  description)
        (1,  0, 1, 0, 0, "Single vertex"),
        (2,  1, 1, 1, 0, "Single edge (K2)"),
        (5,  4, 1, 4, 3, "Path graph (tree)"),
        (6,  4, 1, 4, 1, "Star graph (tree)"),
        (6,  6, 1, 0, 0, "Cycle C6 (biconnected)"),
        (7,  9, 1, 0, 2, "0 bridges, 2 APs"),
        (10, 12, 1, 2, 3, "2 bridges, 3 APs"),
        (8,  7, 1, 7, 6, "Tree with 8 vertices (full forest)"),
        (10, 15, 3, 0, 0, "3 components, 0 bridges"),
        (15, 20, 2, 3, 4, "2 components, 3 bridges, 4 APs"),
    ]

    print("=" * 70)
    print("TESTING CONSTRAINED STRUCTURAL GRAPH GENERATOR (V,E,C,B,A)")
    print("=" * 70)

    for V, E, C, B, A, desc in test_cases:
        print(f"\n--- Test: {desc} ---")
        print(f"    Target: V={V}, E={E}, C={C}, B={B}, A={A}")

        # Kiểm tra tính khả thi
        ok, msg = check_feasibility(V, E, C, B, A)
        if not ok:
            print(f"    Infeasible: {msg}")
            continue

        # Generate graph using Constructive
        G = generate_constructive(V, E, C, B, A)
        if G is None:
            print("    [Constructive] Failed -> trying SA...")
            G = generate_simulated_annealing(V, E, C, B, A,
                                             max_iterations=5000)
            if G is None:
                print("    [SA] Also failed to find a solution.")
                continue

        # Verify using Tarjan
        is_valid, details = verify_graph(G, V, E, C, B, A)
        act = details['actual']
        status = "✓ CORRECT" if is_valid else "✗ WRONG"
        print(f"    Result: V={act['V']}, E={act['E']}, C={act['C']}, "
              f"B={act['B']}, A={act['A']} -> {status}")

        if not is_valid:
            print(f"    Expected: {details['expected']}")
            print(f"    Actual:   {details['actual']}")

    # Test MCMC
    print(f"\n{'=' * 70}")
    print("TESTING MCMC EDGE REWIRING")
    print("=" * 70)
    V, E, C, B, A = 10, 15, 1, 2, 3
    print(f"\nGenerating 3 different graphs with V={V}, E={E}, C={C}, B={B}, A={A}:")
    for i in range(3):
        G = generate_mcmc(V, E, C, B, A, nb_iterations=500)
        if G is not None:
            is_valid, details = verify_graph(G, V, E, C, B, A)
            act = details['actual']
            status = "✓" if is_valid else "✗"
            edges_str = ", ".join(
                f"({e.fromNode}-{e.toNode})" for e in G.edges if e.used == 1)
            print(f"  [{i + 1}] {status} Edges: {edges_str}")
        else:
            print(f"  [{i + 1}] Failed to generate")
