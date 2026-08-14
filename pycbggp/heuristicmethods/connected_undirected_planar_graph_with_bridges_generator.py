import math, sys, os
import random
from matplotlib import pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph, Point2D


ESP = 1e-9 

def orientation(A: Point2D, B: Point2D, C: Point2D):
    orient = (B.x - A.x) * (C.y - B.y) - (B.y - A.y) * (C.x - B.x)
    if orient < -ESP: return -1  # cw
    elif orient > ESP: return 1  # ccw 
    else: return 0               # Collinear

def get_convex_hull(P: list[Point2D]) -> tuple[list[Point2D], int]: 
    S = sorted(P, key=lambda point: (point.x, point.y))
    V = len(S)
    if V <= 2:
        return S, V - 1
        
    def is_all_collinear(Temp): 
        for i in range(2, len(Temp)):
            if orientation(Temp[0], Temp[1], Temp[i]) != 0 :
                return False
        return True
        
    if is_all_collinear(S):
        return S, V - 1

    lower = []
    for p in S:
        while len(lower) >= 2 and orientation(lower[-2], lower[-1], p) < 0:
            lower.pop()
        lower.append(p)
        
    upper = []
    for p in reversed(S):
        while len(upper) >= 2 and orientation(upper[-2], upper[-1], p) < 0:
            upper.pop()
        upper.append(p)
        
    hull = lower[:-1] + upper[:-1]
    max_e = 3 * V - len(hull) - 3
    return hull, max_e # ccw order

def point_in_triangle(p: Point2D, a: Point2D, b: Point2D, c: Point2D) -> bool:
    o1 = orientation(a, b, p)
    o2 = orientation(b, c, p)
    o3 = orientation(c, a, p)
    
    # Check if point strictly sees opposite sides (means it is outside)
    has_neg = (o1 < 0) or (o2 < 0) or (o3 < 0)
    has_pos = (o1 > 0) or (o2 > 0) or (o3 > 0)
    
    # Valid inside or strictly on the edge
    return not (has_neg and has_pos)

class DSU:
    def __init__(self, n):
        self.parent = list(range(n))
    def find(self, i):
        if self.parent[i] == i:
            return i
        self.parent[i] = self.find(self.parent[i])
        return self.parent[i]
    def union(self, i, j):
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            self.parent[root_i] = root_j
            return True
        return False

def random_triangulation(points: list[Point2D]) -> list[tuple[int, int]]:
    n = len(points)
    if n < 3: return []
    
    hull, _ = get_convex_hull(points)
    
    hull_indices = [points.index(p) for p in hull]
    OnHull = set(hull_indices)
    rem_indices = [i for i in range(n) if i not in OnHull]
    
    # Step 1: Initialize triangle framework from convex hull
    triangles = []
    if not rem_indices:
        p_idx = hull_indices[0]
        for i in range(1, len(hull_indices) - 1):
            triangles.append((p_idx, hull_indices[i], hull_indices[i+1]))
    else:
        # Pick a random point as the center fan
        p_idx = random.choice(rem_indices)
        rem_indices.remove(p_idx)
        
        for i in range(len(hull_indices)):
            u = hull_indices[i]
            v = hull_indices[(i + 1) % len(hull_indices)]
            if orientation(points[u], points[v], points[p_idx]) != 0:
                triangles.append((u, v, p_idx))
                
    # Step 2: Random insertion
    random.shuffle(rem_indices)
    
    for i in rem_indices:
        p = points[i]
        bad_triangles = []
        
        for tri in triangles:
            if point_in_triangle(p, points[tri[0]], points[tri[1]], points[tri[2]]):
                bad_triangles.append(tri)
                
        if not bad_triangles:
            continue
            
        new_triangles = []
        for tri in bad_triangles:
            triangles.remove(tri)
            edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]
            
            # Form new triangles
            for u, v in edges:
                if orientation(points[u], points[v], p) != 0:
                    new_triangles.append((u, v, i))
                    
        triangles.extend(new_triangles)
        
    # Step 3: Extract unique edges
    final_edges = set()
    for t1, t2, t3 in triangles:
        final_edges.add((min(t1, t2), max(t1, t2)))
        final_edges.add((min(t2, t3), max(t2, t3)))
        final_edges.add((min(t3, t1), max(t3, t1)))
        
    return list(final_edges)

def heuristic_01(n: int, points: list[Point2D], cand_edges: list[tuple[int, int]], nb_edges: int, nb_bridges: int) -> list[tuple[int, int]]:
    """
    Heuristic 01: Tham lam giảm cầu dựa trên Cây khung nhỏ nhất (EMST) và DSU nén chu trình.
    Trả về danh sách cạnh hợp lệ nếu thành công, trả về None nếu thất bại.
    """
    # Sort candidate edges by Euclidean distance to build EMST 
    cand_edges.sort(key=lambda e: math.hypot(points[e[0]].x - points[e[1]].x, points[e[0]].y - points[e[1]].y))
    
    # 1. Extract Minimum Spanning Tree (MST)
    mst_dsu = DSU(n)
    mst_edges = []
    pool = []
    
    for u, v in cand_edges:
        if mst_dsu.union(u, v):
            mst_edges.append((min(u, v), max(u, v)))
        else:
            pool.append((min(u, v), max(u, v)))
            
    if len(mst_edges) < n - 1: return None
        
    # Adjacency list for path finding on the tree
    mst_adj = {i: [] for i in range(n)}
    for u, v in mst_edges:
        mst_adj[u].append(v)
        mst_adj[v].append(u)
        
    def get_path_in_tree(start, end):
        queue = [[start]]
        visited = {start}
        while queue:
            path = queue.pop(0)
            node = path[-1]
            if node == end: return path
            for neighbor in mst_adj[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(path + [neighbor])
        return []
        
    # Prioritize edges that form small cycles
    pool.sort(key=lambda e: len(get_path_in_tree(e[0], e[1])))
    
    # 2. Bridge Reduction
    bridge_dsu = DSU(n)
    bridges_to_destroy = (n - 1) - nb_bridges
    edges_needed = nb_edges - (n - 1)
    
    final_edges = list(mst_edges)
    final_edges_set = set(mst_edges)
    
    for u, v in pool:
        if bridges_to_destroy == 0 or edges_needed == 0: break
        
        path = get_path_in_tree(u, v)
        
        active_bridges = 0
        for i in range(len(path) - 1):
            if bridge_dsu.find(path[i]) != bridge_dsu.find(path[i+1]):
                active_bridges += 1
                
        if 0 < active_bridges <= bridges_to_destroy:
            final_edges.append((u, v))
            final_edges_set.add((min(u, v), max(u, v)))
            edges_needed -= 1
            bridges_to_destroy -= active_bridges
            
            for i in range(len(path) - 1):
                bridge_dsu.union(path[i], path[i+1])
                
    # 3. Edge Augmentation
    for u, v in pool:
        if edges_needed == 0: break
        if (min(u, v), max(u, v)) not in final_edges_set and bridge_dsu.find(u) == bridge_dsu.find(v):
            final_edges.append((u, v))
            final_edges_set.add((min(u, v), max(u, v)))
            edges_needed -= 1
            
    # 4. Success check
    if bridges_to_destroy == 0 and edges_needed == 0:
        return final_edges

    return None

def heuristic_02(n: int, points: list[Point2D], cand_edges: list[tuple[int, int]], nb_edges: int, nb_bridges: int) -> list[tuple[int, int]]:

    return None

def generate_connected_planar_graph_with_bridges(points: list[Point2D], nb_edges, nb_bridges):
    assert len({(p.x, p.y) for p in points}) == len(points), "!!!!Duplicate points detected in the input list!!!!"
    assert len({p.id for p in points}) == len(points), "!!!!Duplicate index detected in the input list!!!!"
    n = len(points)
    
    if n == 0: return None if nb_edges > 0 else Graph(0) 
    if n == 1: return Graph(1) if nb_edges == 0 else None
    if n == 2 and (nb_edges != 1 or nb_bridges != 1): return None
    
    # Constraints check
    hull, max_e = get_convex_hull(points)
    if nb_edges < n - 1 or nb_edges > max_e: return None
    if nb_bridges > n - 1 or nb_bridges < 0: return None

    # TODO by Nguyen Phuc Khanh 
    
    MAX_ATTEMPTS = 50 
    for attempt in range(MAX_ATTEMPTS):
        # Triangulation to create a maximum planar graph
        cand_edges = random_triangulation(points)
        if len(cand_edges) < nb_edges: continue
        
        final_edges = heuristic_01(n, points, cand_edges, nb_edges, nb_bridges)
        
        
        if final_edges is not None:
            G = Graph(n)
            for u, v in final_edges:
                G.AddEdge(u, v) 
            # print(attempt, file=sys.stderr)
            return G

    return None


def plot_graph(G, points):
    plt.figure(figsize = (8, 6))
    
    idx_to_pos = {idx: (p.x, p.y) for idx, p in enumerate(points)}
    
    if G != None:
        for e in G.edges:
            u_idx, v_idx = e.fromNode, e.toNode
            if u_idx in idx_to_pos and v_idx in idx_to_pos:
                u_pos, v_pos = idx_to_pos[u_idx], idx_to_pos[v_idx]
                plt.plot([u_pos[0], v_pos[0]], [u_pos[1], v_pos[1]], 'b-', alpha=0.6)
            
    x_coords, y_coords = [p.x for p in points], [p.y for p in points]
    plt.scatter(x_coords, y_coords, c = 'red', s = 100, zorder = 5)
    
    for p in points:
        plt.annotate(f"{p.id}", (p.x, p.y), textcoords = "offset points", xytext = (5,5), ha = 'center', fontsize = 12, fontweight = 'bold')
    stauts = ""
    if not G: Title = "Not complete"  
    else : Title = f"Connected Planar Graph (V = {G.n}, E = {G.m})"
    plt.title(Title)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True, linestyle='--', alpha = 0.5)
    plt.show()

if __name__ == "__main__":
    points = [
        Point2D(0, 5.5, 10.0), Point2D(1, 10.0, 0.0), 
        Point2D(2, 0.0, -10.0), Point2D(3, -10.0, 0.0),
        
        Point2D(4, 3.0, 3.0), Point2D(5, 3.0, -3.0), 
        Point2D(6, -3.0, -3.0), Point2D(7, -3.0, 3.0),
        
        Point2D(8, 0.0, 0.0)
    ]
    target_edges = 20
    for _ in range(5):
        G = generate_connected_planar_graph_with_bridges(points, 13, 3)
        plot_graph(G, points)
    