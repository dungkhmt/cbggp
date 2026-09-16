import math, sys, os
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph, Point2D

EPS = 1e-9

def zigzag_connect(layer, add_edge):
    L = len(layer)
    if L <= 3: return
    
    left, right = 1, L - 1
    f = True
    while right - left > 1:
        add_edge(layer[left], layer[right])
        if f: left += 1
        else: right -= 1
        f = not f

def connect_layers(O, I, add_edge, cross):
    cx = sum(p.x for p in I) / len(I) + 1e-7
    cy = sum(p.y for p in I) / len(I) + 2.718e-7
    
    def ang(p): return math.atan2(p.y - cy, p.x - cx)
    
    I_angs = [ang(p) for p in I]
    O_angs = [ang(p) for p in O]
    
    def ang_dist(a, b):
        diff = (a - b) % (2 * math.pi)
        return min(diff, 2 * math.pi - diff)

    best_j = 0
    min_dist = float('inf')
    for j, oa in enumerate(O_angs):
        d = ang_dist(I_angs[0], oa)
        if d < min_dist:
            min_dist = d
            best_j = j

    O_rot = O[best_j:] + O[:best_j]
    O_angs_rot = O_angs[best_j:] + O_angs[:best_j]

    def unwrap(angles, base_angle=None):
        if not angles: return []
        unwrapped = [angles[0]]
        if base_angle is not None:
            diff = (unwrapped[0] - base_angle) % (2 * math.pi)
            if diff > math.pi: diff -= 2 * math.pi
            unwrapped[0] = base_angle + diff
            
        for k in range(1, len(angles)):
            diff = (angles[k] - unwrapped[-1]) % (2 * math.pi)
            if diff > math.pi:
                diff -= 2 * math.pi
            if diff <= 0:
                diff = 1e-6 
            unwrapped.append(unwrapped[-1] + diff)
        return unwrapped

    I_unwrapped = unwrap(I_angs)
    O_unwrapped = unwrap(O_angs_rot, I_unwrapped[0])

    I_rot = I[:] + [I[0]]
    I_unwrapped.append(I_unwrapped[0] + 2 * math.pi)
    
    O_rot = O_rot[:] + [O_rot[0]]
    O_unwrapped.append(O_unwrapped[0] + 2 * math.pi)

    i = j = 0
    add_edge(I_rot[0], O_rot[0])
    
    while i < len(I) or j < len(O):
        if i == len(I) and j == len(O):
            break
        elif i == len(I):
            j += 1
        elif j == len(O):
            i += 1
        else:
            c1 = cross(O_rot[j], I_rot[i + 1], I_rot[i])
            c2 = cross(I_rot[i], O_rot[j], O_rot[j + 1])
            
            valid_i = c1 > EPS
            valid_j = c2 > EPS
            
            if valid_i and not valid_j:
                i += 1
            elif valid_j and not valid_i:
                j += 1
            else:
                if I_unwrapped[i + 1] <= O_unwrapped[j + 1]:
                    i += 1
                else:
                    j += 1
        add_edge(I_rot[i], O_rot[j])


def _gen_connected_planar_graph(points, nb_edges):
    n = len(points)
    if n == 0: return None if nb_edges > 0 else Graph(0) 
    if n == 1: return Graph(1) if nb_edges == 0 else None

    id_to_idx = {p.id: idx for idx, p in enumerate(points)}

    def cross(A, B, C):
        return (B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y)

    def get_convex_hull(pts):
        pts = sorted(pts, key = lambda p: (p.x, p.y))
        if len(pts) <= 2: return pts[:]
        
        lower = []
        for p in pts:
            while len(lower) >= 2 and cross(lower[-2], lower[-1], p) < -EPS:
                lower.pop()
            lower.append(p)
            
        upper = []
        for p in reversed(pts):
            while len(upper) >= 2 and cross(upper[-2], upper[-1], p) < -EPS:
                upper.pop()
            upper.append(p)
            
        hull = lower[:-1] + upper[:-1]
        return hull if hull else pts[:]

    remaining = points[:]
    layers = []
    while remaining:
        hull = get_convex_hull(remaining)
        layers.append(hull)
        ids = {p.id for p in hull}
        remaining = [p for p in remaining if p.id not in ids]

    k = len(layers[0])
    max_edges = 3 * n - k - 3 if n >= 3 else (1 if n == 2 else 0)
    min_edges = n - 1
    
    if nb_edges < min_edges or nb_edges > max_edges:
        return None

    edge_pool = set()
    def add_edge(u, v):
        if u.id != v.id:
            edge_pool.add(tuple(sorted([u, v], key = lambda p: p.id)))

    for layer_idx, layer in enumerate(layers):
        L = len(layer)
        if L == 2:
            add_edge(layer[0], layer[1])
        elif L >= 3:
            is_coll = True
            for k_idx in range(2, L):
                if abs(cross(layer[0], layer[1], layer[k_idx])) > EPS:
                    is_coll = False
                    break
            
            if is_coll:
                for idx in range(L - 1): add_edge(layer[idx], layer[idx + 1])
            else:
                for idx in range(L): add_edge(layer[idx], layer[(idx + 1) % L])
                if layer_idx == len(layers) - 1:
                    zigzag_connect(layer, add_edge)

    for layer_idx in range(len(layers) - 1):
        O = layers[layer_idx]
        I = layers[layer_idx + 1]
        
        if len(I) == 1:
            for p in O: add_edge(I[0], p)
        elif len(O) == 1:
            for p in I: add_edge(p, O[0])
        else:
            connect_layers(O, I, add_edge, cross)

    pool_list = sorted(list(edge_pool), key = lambda x: (x[0].id, x[1].id))
    parent = {p.id: p.id for p in points}
    
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        root_i = find(i)
        root_j = find(j)
        if root_i != root_j:
            parent[root_i] = root_j
            return True
        return False

    final_edges = []
    unused_edges = []
    
    for u, v in pool_list:
        if union(u.id, v.id):
            final_edges.append((u, v))
        else:
            unused_edges.append((u, v))
            
    if len(final_edges) < min_edges or len(final_edges) + len(unused_edges) < nb_edges:
        return None
        
    needed = nb_edges - len(final_edges)
    final_edges.extend(unused_edges[:needed])

    G = Graph(n)
    for u, v in final_edges:
        G.AddEdge(id_to_idx[u.id], id_to_idx[v.id])
        
    return G


def plot_graph(G, points):
    if not G: return
    plt.figure(figsize = (8, 6))
    
    idx_to_pos = {idx: (p.x, p.y) for idx, p in enumerate(points)}
    
    for e in G.edges:
        u_idx, v_idx = e.fromNode, e.toNode
        if u_idx in idx_to_pos and v_idx in idx_to_pos:
            u_pos, v_pos = idx_to_pos[u_idx], idx_to_pos[v_idx]
            plt.plot([u_pos[0], v_pos[0]], [u_pos[1], v_pos[1]], 'b-', alpha=0.6)
            
    x_coords, y_coords = [p.x for p in points], [p.y for p in points]
    plt.scatter(x_coords, y_coords, c = 'red', s = 100, zorder = 5)
    
    for p in points:
        plt.annotate(f"{p.id}", (p.x, p.y), textcoords = "offset points", xytext = (5,5), ha = 'center', fontsize = 12, fontweight = 'bold')
        
    plt.title(f"Connected Planar Graph (V = {G.n}, E = {G.m})")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True, linestyle='--', alpha = 0.5)
    plt.show()

if __name__ == "__main__":
    # points = [
    #     Point2D(0, 0.0, 0.0), Point2D(1, 10.0, 0.0), Point2D(2, 5.0, 10.0),
    #     Point2D(3, 5.0, 0.00000001),
    #     Point2D(4, 5.0, 5.0)
    # ]
    # target_edges = 9

    # points = [
    #     Point2D(0, 0.0, 10.0), Point2D(1, 10.0, 0.0), 
    #     Point2D(2, 0.0, -10.0), Point2D(3, -10.0, 0.0),
        
    #     Point2D(4, 3.0, 3.0), Point2D(5, 3.0, -3.0), 
    #     Point2D(6, -3.0, -3.0), Point2D(7, -3.0, 3.0),
        
    #     Point2D(8, 0.0, 0.0)
    # ]
    # target_edges = 20

    points = [
        Point2D(10, 0.0, 0.0), Point2D(55, 10.0, 0.0), Point2D(99, 5.0, 10.0),
        Point2D(105, 5.0, 0.00000001),
        Point2D(202, 5.0, 5.0)
    ]
    target_edges = 9

    G = gen_connected_planar_graph(points, target_edges) 
    
    if G:
        print("OK!")
        G.Print()
        plot_graph(G, points)
    else:
        print("!OK")
