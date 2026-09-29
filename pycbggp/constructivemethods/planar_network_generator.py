import sys, os, math, random
from typing import List, Optional, Tuple
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import DirectedGraph, Point2D
from DSU import DSU


def _build_points(nb_nodes: int, W: int, H: int):
    """ Node 0 = source at x=0, node n-1 = sink at x=W-1, both vertically centered.
        The n-2 remaining nodes are spread over up to W-2 middle columns, each
        column's points evenly spaced across the full height (spacious, no
        cramped/acute edges), and the columns themselves evenly spaced in x when
        there are fewer points than available column slots. Returns (points,
        columns) with columns = list of column id-lists ordered left to right
        (source and sink each their own 1-node "column"), or None if infeasible. """
    if W < 2:
        return None

    mid = nb_nodes - 2
    if mid == 0:
        return [Point2D(0, 0, H // 2), Point2D(1, W - 1, H // 2)], [[0], [1]]
    if W < 3:
        return None

    C = min(W - 2, mid)
    base, extra = divmod(mid, C)
    sizes = [base + (1 if c < extra else 0) for c in range(C)]
    if max(sizes) > H:
        return None

    xs = [(W - 1) // 2] if C == 1 else [1 + c * (W - 3) // (C - 1) for c in range(C)]

    points = [Point2D(0, 0, H // 2)]
    columns = [[0]]
    pid = 1
    for x, size in zip(xs, sizes):
        ys = [H // 2] if size == 1 else [k * (H - 1) // (size - 1) for k in range(size)]
        columns.append(list(range(pid, pid + size)))
        for y in ys:
            points.append(Point2D(pid, x, y))
            pid += 1
    points.append(Point2D(pid, W - 1, H // 2))
    columns.append([pid])
    return points, columns


def _strip_edges(colA: List[int], colB: List[int], points: List[Point2D]) -> List[Tuple[int, int]]:
    """ Non-crossing "ladder" between two adjacent columns (both already sorted by
        y): a two-pointer merge that always advances toward whichever side's next
        point is closer in y, so the connecting lines never cross. """
    i, j = 0, 0
    edges = [(colA[0], colB[0])]
    while i < len(colA) - 1 or j < len(colB) - 1:
        if i == len(colA) - 1:
            j += 1
        elif j == len(colB) - 1:
            i += 1
        elif points[colA[i + 1]].y <= points[colB[j + 1]].y:
            i += 1
        else:
            j += 1
        edges.append((colA[i], colB[j]))
    return edges


def _build_spanning_tree(edges: set, n: int) -> set:
    edges_list = list(edges)
    random.shuffle(edges_list)
    dsu = DSU(n)
    protected = set()
    for (u, v) in edges_list:
        if dsu.union(u, v):
            protected.add((u, v))
            if len(protected) == n - 1:
                break
    return protected


def _prune_to_target(edges: set, protected: set, k: int, points: List[Point2D], deg: List[int]) -> set:
    if k <= 0:
        return set(edges)
    candidates = list(edges - protected)
    avg_len = sum(math.hypot(points[u].x - points[v].x, points[u].y - points[v].y) for u, v in edges) / len(edges)

    def score(e):
        u, v = e
        length = math.hypot(points[u].x - points[v].x, points[u].y - points[v].y)
        return (deg[u] + deg[v]) + length / avg_len

    candidates.sort(key=score, reverse=True)
    result = set(edges)
    for (u, v) in candidates[:k]:
        result.discard((u, v))
    return result


def gen_planar_network(nb_nodes: int, nb_arcs: int, W: int, H: int) -> Optional[Tuple[DirectedGraph, List[Point2D]]]:
    if nb_nodes == 0:
        return (DirectedGraph(0), []) if nb_arcs == 0 else None
    if nb_nodes == 1:
        return (DirectedGraph(1), [Point2D(0, 0, H // 2)]) if nb_arcs == 0 else None

    built = _build_points(nb_nodes, W, H)
    if built is None:
        return None
    points, columns = built

    pool = set()
    for colA, colB in zip(columns, columns[1:]):
        for u, v in _strip_edges(colA, colB, points):
            pool.add((u, v) if u < v else (v, u))

    m_min, m_max = nb_nodes - 1, len(pool)
    if nb_arcs < m_min or nb_arcs > m_max:
        return None

    protected = _build_spanning_tree(pool, nb_nodes)
    deg = [0] * nb_nodes
    for u, v in pool:
        deg[u] += 1
        deg[v] += 1
    final_edges = _prune_to_target(pool, protected, len(pool) - nb_arcs, points, deg)

    G = DirectedGraph(nb_nodes)
    for u, v in final_edges:
        G.AddEdge(u, v)  # u < v always: edges only ever cross between adjacent columns, left to right
    return G, points


def plot_graph(G: DirectedGraph, points: List[Point2D]):
    if not G: return
    plt.figure(figsize=(8, 6))

    idx_to_pos = {idx: (p.x, p.y) for idx, p in enumerate(points)}
    for e in G.edges:
        plt.annotate('', xy=idx_to_pos[e.toNode], xytext=idx_to_pos[e.fromNode],
                      arrowprops=dict(arrowstyle='->', color='b', alpha=0.6, shrinkA=8, shrinkB=8))

    colors = ['green' if idx == 0 else 'orange' if idx == G.n - 1 else 'red' for idx in range(G.n)]
    plt.scatter([p.x for p in points], [p.y for p in points], c=colors, s=100, zorder=5)
    for p in points:
        tag = " (S)" if p.id == 0 else " (T)" if p.id == G.n - 1 else ""
        plt.annotate(f"{p.id}{tag}", (p.x, p.y), textcoords="offset points", xytext=(5, 5), ha='center', fontsize=12, fontweight='bold')

    plt.title(f"Planar Network (V = {G.n}, E = {G.m})")
    plt.xlabel("X"); plt.ylabel("Y")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


def _segments_cross(p1, p2, p3, p4) -> bool:
    """ Self-check only -- true if segments (p1,p2) and (p3,p4) overlap/cross
        anywhere other than at a shared endpoint. """
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def on_segment(p, a, b):
        return min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])

    if {p1, p2} & {p3, p4}:
        shared = ({p1, p2} & {p3, p4}).pop()
        others = [p for p in (p1, p2, p3, p4) if p != shared]
        return cross(shared, others[0], others[1]) == 0 and on_segment(others[0], shared, others[1])

    d1, d2 = cross(p3, p4, p1), cross(p3, p4, p2)
    d3, d4 = cross(p1, p2, p3), cross(p1, p2, p4)
    if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
        return True
    for p, a, b in ((p1, p3, p4), (p2, p3, p4), (p3, p1, p2), (p4, p1, p2)):
        if cross(a, b, p) == 0 and on_segment(p, a, b):
            return True
    return False


if __name__ == "__main__":
    result = gen_planar_network(12, 15, 6, 5)
    assert result is not None
    G, points = result
    assert G.n == 12 and G.m == 15

    in_deg = [0] * G.n
    for e in G.edges:
        assert e.fromNode < e.toNode, "every edge must point from a smaller id to a larger id"
        in_deg[e.toNode] += 1
    assert in_deg[0] == 0, "node 0 must be the unique source"
    assert len(G.Adj[G.n - 1]) == 0, "node n-1 must be the unique sink"

    pos = [(p.x, p.y) for p in points]
    segs = [(pos[e.fromNode], pos[e.toNode]) for e in G.edges]
    for i in range(len(segs)):
        for j in range(i + 1, len(segs)):
            assert not _segments_cross(*segs[i], *segs[j]), f"edges overlap: {segs[i]} vs {segs[j]}"

    parent = list(range(G.n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for e in G.edges:
        ru, rv = find(e.fromNode), find(e.toNode)
        if ru != rv:
            parent[ru] = rv
    assert len({find(i) for i in range(G.n)}) == 1, "must be connected"

    assert gen_planar_network(100, 10, 5, 5) is None
    assert gen_planar_network(0, 0, 3, 3) is not None
    assert gen_planar_network(1, 0, 3, 3) is not None

    print("OK!")
    plot_graph(G, points)
