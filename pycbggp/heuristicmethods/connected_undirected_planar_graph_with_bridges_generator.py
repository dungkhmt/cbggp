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

def heuristic_02(V: int, points: list[Point2D],  E: int, B: int) -> list[tuple[int, int]]:
    # Unified Monotone Core-Branch
    '''
    Step1 : Divide points into 2 parts: Core & Leaves
            V - B points belong to Core & B remain points are leaves
    Step2:  Assign V - B edges for Core to create a planar graph without bridge
            Then connect B points and Core => we have V edges & B bridges
    Step3:  We need E edges so create more edges from Core 
    '''
    
    
    def is_cross(u1, v1, u2, v2):
        p1, p2, p3, p4 = points[u1], points[v1], points[u2], points[v2]
        def on_seg(p, a, b):
            return min(a.x, b.x) <= p.x <= max(a.x, b.x) and min(a.y, b.y) <= p.y <= max(a.y, b.y)

        if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
            shared = p1 if (u1 == u2 or u1 == v2) else p2
            o1 = p2 if shared == p1 else p1
            o2 = p4 if shared == p3 else p3
            if orientation(shared, o1, o2) == 0:
                if on_seg(o2, shared, o1) or on_seg(o1, shared, o2): return True
            return False
            
        o1, o2 = orientation(p1, p2, p3), orientation(p1, p2, p4)
        o3, o4 = orientation(p3, p4, p1), orientation(p3, p4, p2)
        
        if o1 != o2 and o3 != o4 and o1 != 0 and o2 != 0 and o3 != 0 and o4 != 0: return True
        if o1 == 0 and on_seg(p3, p1, p2): return True
        if o2 == 0 and on_seg(p4, p1, p2): return True
        if o3 == 0 and on_seg(p1, p3, p4): return True
        if o4 == 0 and on_seg(p2, p3, p4): return True
        return False

    if B > V - 3 and B != V - 1: return None
    
    # change the coordinate 
    # Apply a random projection matrix to change the sorting axis. 
    # This diversifies the X-monotone shapes generated across retries.
    cof = [[1, 0, 0, 1], [0, 1, 1, 0], [-1, 0, 0, 1], [0, -1, 1, 0], 
           [1, 1, 1, -1], [-1, -1, 1, -1], [-1, 1, 1, 1], [1, -1, -1, -1],
           [9, 2, 0, 7], [2, 0, 0, 7]]
    J = random.choice(cof)
    pts = sorted(range(V), key=lambda i: (points[i].x * J[0] + points[i].y * J[1], points[i].x * J[2] + points[i].y * J[3]))
        
    
    # Incase: Tree
    if B == V - 1:
        if E != V - 1: return None
        return [(min(pts[i], pts[i+1]), max(pts[i], pts[i+1])) for i in range(V - 1)]

    # Step 1:
    # Split the sorted vertices into Core (polygon) and Leaves (bridges)
    core_sz = V - B
    core = pts[:core_sz]
    leaf = pts[core_sz:]
    
    # Step 2:
    # Forming core 
    # Separate the core points into Upper and Lower bounds using orientation
    p_min, p_max = core[0], core[-1]
    up, dn = [], []
    for i in core[1:-1]:
        if orientation(points[p_min], points[p_max], points[i]) > 0:
            up.append(i)
        else:
            dn.append(i)
    cycle = [p_min] + up + [p_max] + dn[::-1] # Combine bounds to form a simple Monotone Polygon (0 bridges internally)
    
    edges = []
    for i in range(len(cycle)):
        u, v = cycle[i], cycle[(i + 1) % len(cycle)]
        edges.append((min(u, v), max(u, v)))
    # Connecting core with leaves    
    curr = p_max
    for L in leaf:
        edges.append((min(curr, L), max(curr, L)))
        curr = L
        
    rem = E - V
    if rem < 0: return None
    if rem == 0: return edges
    
    # Step 3: Ear Clipping Augmentation
    Upper = [p_min] + up + [p_max]
    Lower = [p_min] + dn + [p_max] 
    if rem == 0: return edges
    i = 1
    while i < len(Upper) - 1 and rem > 0:
        u, mid, v = Upper[i - 1], Upper[i], Upper[i + 1]
        if orientation(points[u], points[mid], points[v]) < 0:
            e = (min(u, v), max(u, v))
            conflict = False
            for eu, ev in edges:
                if is_cross(u, v, eu, ev):
                    conflict = True
                    break
            if not conflict and e not in edges:
                edges.append(e)
                rem -= 1
                Upper.pop(i) 
                i = max(1, i - 1)
                continue
        i += 1
    i = 1
    while i < len(Lower) - 1 and rem > 0:
        u, mid, v = Lower[i - 1], Lower[i], Lower[i + 1]
        if orientation(points[u], points[mid], points[v]) < 0:
            e = (min(u, v), max(u, v))
            conflict = False
            for eu, ev in edges:
                if is_cross(u, v, eu, ev):
                    conflict = True
                    break
            if not conflict and e not in edges:
                edges.append(e)
                rem -= 1
                Lower.pop(i) 
                i = max(1, i - 1) 
                continue
        i += 1
    if rem == 0:
        return edges
    return None
    

MAX_EDGES = 0
TYPE = 0
def generate_connected_planar_graph_with_bridges(points: list[Point2D], nb_edges, nb_bridges):
    assert len({(p.x, p.y) for p in points}) == len(points), "!!!!Duplicate points detected in the input list!!!!"
    assert len({p.id for p in points}) == len(points), "!!!!Duplicate index detected in the input list!!!!"
    n = len(points)
    
    if n == 0: return None if nb_edges > 0 else Graph(0) 
    if n == 1: return Graph(1) if nb_edges == 0 else None
    if n == 2 and (nb_edges != 1 or nb_bridges != 1): return None
    
    # Constraints check
    hull, max_e = get_convex_hull(points)
    global MAX_EDGES 
    global TYPE
    MAX_EDGES = max_e
    if nb_edges < n - 1 or nb_edges > max_e: return None
    if nb_bridges > n - 1 or nb_bridges < 0 or nb_bridges == n - 2: return None
    if nb_edges == n - 1 and nb_bridges != n - 1: return None
    if nb_edges != n - 1 and nb_bridges == n - 1: return None
    if nb_edges > n - 1 and nb_edges - nb_bridges < 3: return None
    if nb_edges == max_e and nb_bridges == 0: return None
    
    valid_edges = None
    for _ in range(30):
        if valid_edges: break
        cand_edges = random_triangulation(points)
        valid_edges = heuristic_01(n, points, cand_edges, nb_edges, nb_bridges)
        TYPE = 1
    for _ in range(30):
        if valid_edges: break
        valid_edges = heuristic_02(n, points, nb_edges, nb_bridges)
        TYPE = 2
    if valid_edges:
        G = Graph(n)
        for u, v in valid_edges:
            G.AddEdge(u, v) 
        return G

    return None


def plot_graph(G, points):
    plt.figure(figsize = (9, 6))
    
    idx_to_pos = {idx: (p.x, p.y) for idx, p in enumerate(points)}
    
    if G != None:
        for e in G.edges:
            u_idx, v_idx = e.fromNode, e.toNode
            if u_idx in idx_to_pos and v_idx in idx_to_pos:
                u_pos, v_pos = idx_to_pos[u_idx], idx_to_pos[v_idx]
                plt.plot([u_pos[0], v_pos[0]], [u_pos[1], v_pos[1]], 'b-', alpha=0.6)
            
    x_coords, y_coords = [p.x for p in points], [p.y for p in points]
    plt.scatter(x_coords, y_coords, c = 'red', s = 20, zorder = 5)
    
    for p in points:
        plt.annotate(f"{p.id}", (p.x, p.y), textcoords = "offset points", xytext = (5,5), ha = 'center', fontsize = 12, fontweight = 'bold')
    if not G: Title = f"Not complete (MAXE = {MAX_EDGES})"  
    else : Title = f"Connected Planar Graph (V = {G.n}, E = {G.m}, MAXE = {MAX_EDGES}) TYPE#{TYPE}"
    plt.title(Title)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True, linestyle='--', alpha = 0.5)
    plt.show()

# def gen_points(n: int, pattern: str = "random", bounds: tuple = (0, 100)) -> list[Point2D]:
#     """
#     Generate 2D points based on specific geometric patterns.
#     Supported patterns: "random", "dense", "sparse", "linear", "circular"
#     """
#     if n <= 0: return []
    
#     points = []
#     used_coords = set()
#     min_val, max_val = bounds
    
#     def add_point(id_val, x, y):
#         # Round to 3 decimal places to avoid floating point chaos and check duplicates
#         x, y = round(x, 3), round(y, 3)
#         if (x, y) not in used_coords:
#             used_coords.add((x, y))
#             points.append(Point2D(id_val, x, y))
#             return True
#         return False

#     if pattern == "random":
#         # Uniform random distribution
#         while len(points) < n:
#             add_point(len(points), random.uniform(min_val, max_val), random.uniform(min_val, max_val))
            
#     elif pattern == "dense":
#         # Gaussian distribution around a center (Clustered)
#         center_x = (min_val + max_val) / 2
#         center_y = (min_val + max_val) / 2
#         sigma = (max_val - min_val) / 6  # 99% of points within bounds
        
#         while len(points) < n:
#             x = random.gauss(center_x, sigma)
#             y = random.gauss(center_y, sigma)
#             add_point(len(points), x, y)
            
#     elif pattern == "sparse":
#         # Poisson-disk style (Minimum distance rejection sampling)
#         area = (max_val - min_val) ** 2
#         # Expected max radius to keep them sparse but fit n points
#         min_dist = math.sqrt(area / n) * 0.7 
        
#         attempts = 0
#         while len(points) < n:
#             x = random.uniform(min_val, max_val)
#             y = random.uniform(min_val, max_val)
            
#             # Check distance against all existing points
#             too_close = False
#             for p in points:
#                 if math.hypot(p.x - x, p.y - y) < min_dist:
#                     too_close = True
#                     break
                    
#             if not too_close or attempts > 100: # Force add if stuck
#                 if add_point(len(points), x, y):
#                     attempts = 0
#             else:
#                 attempts += 1
#                 if attempts > 50: min_dist *= 0.9 # Relax condition if getting stuck
                
#     elif pattern == "linear":
#         # Points along a line with Gaussian noise
#         slope = random.uniform(-2, 2)
#         intercept = random.uniform(min_val, max_val)
#         noise_level = (max_val - min_val) * 0.05
        
#         while len(points) < n:
#             # Distribute X evenly, add noise to Y
#             progress = len(points) / n
#             x_base = min_val + progress * (max_val - min_val)
#             x = x_base + random.gauss(0, noise_level)
#             y = (slope * x + intercept) + random.gauss(0, noise_level)
#             add_point(len(points), x, y)
            
#     elif pattern == "circular":
#         # Points along the perimeter of a circle with slight noise
#         center_x = (min_val + max_val) / 2
#         center_y = (min_val + max_val) / 2
#         radius = (max_val - min_val) * 0.4
#         noise_level = radius * 0.05
        
#         while len(points) < n:
#             angle = random.uniform(0, 2 * math.pi)
#             r = radius + random.gauss(0, noise_level)
#             x = center_x + r * math.cos(angle)
#             y = center_y + r * math.sin(angle)
#             add_point(len(points), x, y)
            
#     else:
#         raise ValueError("Unknown pattern. Choose from: random, dense, sparse, linear, circular")
        
#     return points

# if __name__ == "__main__":
#     patterns = ["random", "dense", "sparse", "linear", "circular"]
#     points = gen_points(n = 50, pattern = "sparse", bounds=(0, 100))
#     target_edges = 20
#     for _ in range(10):
#         G = generate_connected_planar_graph_with_bridges(points, 60, 10)
#         plot_graph(G, points)

    