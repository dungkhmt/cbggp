import math, sys, os, random, heapq, collections
from matplotlib import pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph, Point2D

def get_part(V, k):
    if V < k: return None
    for _ in range(1000):
        p = [1] * k
        for _ in range(V - k): p[random.randint(0, k - 1)] += 1
        if all(x != 2 for x in p): return p
    return None

def split_rect(W, H, k):
    hp = [(-(W * H), 0, W, 0, H)]
    for _ in range(k - 1):
        _, x0, x1, y0, y1 = heapq.heappop(hp)
        w, h = x1 - x0, y1 - y0
        cd = 1 if w < 2 else (0 if h < 2 else random.choice([0, 1]))
        
        if cd == 0:
            xc = random.randint(x0 + 1, x1 - 1)
            a1, a2 = (xc - x0) * h, (x1 - (xc + 1)) * h
            heapq.heappush(hp, (-a1, x0, xc, y0, y1))
            heapq.heappush(hp, (-a2, xc + 1, x1, y0, y1))
        else:
            yc = random.randint(y0 + 1, y1 - 1)
            a1, a2 = w * (yc - y0), w * (y1 - (yc + 1))
            heapq.heappush(hp, (-a1, x0, x1, y0, yc))
            heapq.heappush(hp, (-a2, x0, x1, yc + 1, y1))
    return [(r[1], r[2], r[3], r[4]) for r in hp]

def gen_pts(V, B, W, H):
    k = B + 1
    p = get_part(V, k)
    if not p: raise ValueError(f"Cannot split {V} nodes into {k} parts.")
        
    rects = split_rect(W, H, k)
    used, pts = set(), []
    
    for i, n in enumerate(p):
        x0, x1, y0, y1 = rects[i]
        if (x1 - x0 + 1) * (y1 - y0 + 1) < n:
            raise ValueError(f"Rect {i} too small.")
            
        cur, att = 0, 0
        while cur < n:
            x, y = random.randint(x0, x1), random.randint(y0, y1)
            if (x, y) not in used:
                used.add((x, y))
                pts.append((x, y, i))
                cur += 1
            att += 1
            if att > n * 20: raise ValueError("Space constrained.")
    return pts

def cp(A, B, C):
    return (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)

def on_seg(p, q, r):
    return min(p.x, r.x) <= q.x <= max(p.x, r.x) and min(p.y, r.y) <= q.y <= max(p.y, r.y)

def isect(A, B, C, D):
    d1, d2, d3, d4 = cp(C, D, A), cp(C, D, B), cp(A, B, C), cp(A, B, D)
    if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)): return True
    if d1 == 0 and on_seg(C, A, D): return True
    if d2 == 0 and on_seg(C, B, D): return True
    if d3 == 0 and on_seg(A, C, B): return True
    if d4 == 0 and on_seg(A, D, B): return True
    return False

def check_isect(G, u, v):
    A, B = G.nodes[u], G.nodes[v]
    for e in G.edges:
        if u in (e.fromNode, e.toNode) or v in (e.fromNode, e.toNode): continue
        if isect(A, B, G.nodes[e.fromNode], G.nodes[e.toNode]): return True
    return False

def gen_graph(V, E, B, W, H):
    if E > 3 * V - 2 * B - 6: 
        print("Err: Max edges exceeded.")
        return None
    try:
        pts = gen_pts(V, B, W, H)
    except ValueError as e:
        print(e)
        return None

    G = Graph(V)
    grps = collections.defaultdict(list)
    G.nodes = []
    
    for i, (x, y, gid) in enumerate(pts):
        G.nodes.append(Point2D(i, x, y))
        grps[gid].append(i)
        
    # 1. Monotone cycles
    for gid, nodes in grps.items():
        if len(nodes) <= 2: continue 
        nodes.sort(key=lambda idx: (G.nodes[idx].x, G.nodes[idx].y))
        
        ps, pe = nodes[0], nodes[-1]
        up, dn = [], []
        x1, y1 = G.nodes[ps].x, G.nodes[ps].y
        x2, y2 = G.nodes[pe].x, G.nodes[pe].y
        
        for i in range(1, len(nodes) - 1):
            idx = nodes[i]
            if (x2 - x1) * (G.nodes[idx].y - y1) - (y2 - y1) * (G.nodes[idx].x - x1) > 0: up.append(idx)
            else: dn.append(idx)
                
        cyc = [ps] + up + [pe] + dn[::-1]
        for i in range(len(cyc)): G.AddEdge(cyc[i], cyc[(i + 1) % len(cyc)])
            
    # 2. MST bridges
    cents = {g: (sum(G.nodes[n].x for n in ns) / len(ns), sum(G.nodes[n].y for n in ns) / len(ns)) for g, ns in grps.items()}
    gids = list(grps.keys())
    ce = sorted((math.hypot(cents[gids[i]][0] - cents[gids[j]][0], cents[gids[i]][1] - cents[gids[j]][1]), gids[i], gids[j]) 
                for i in range(len(gids)) for j in range(i + 1, len(gids)))
    
    par = {g: g for g in gids}
    def find(i):
        if par[i] == i: return i
        par[i] = find(par[i])
        return par[i]
        
    for _, u_g, v_g in ce:
        ru, rv = find(u_g), find(v_g)
        if ru != rv:
            pairs = sorted((math.hypot(G.nodes[u].x - G.nodes[v].x, G.nodes[u].y - G.nodes[v].y), u, v) 
                           for u in grps[u_g] for v in grps[v_g])
            ba = False
            for _, u, v in pairs:
                if not check_isect(G, u, v):
                    G.AddEdge(u, v)
                    par[ru] = rv 
                    ba = True
                    break
            if not ba and pairs:
                G.AddEdge(pairs[0][1], pairs[0][2])
                par[ru] = rv

    # 3. Fill edges
    en = E - G.m
    if en > 0:
        for gid, nodes in grps.items():
            if en <= 0: break
            if len(nodes) < 4: continue 
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    if en <= 0: break
                    u, v = nodes[i], nodes[j]
                    if not any(((G.edges[eid].fromNode == u and G.edges[eid].toNode == v) or 
                                (G.edges[eid].fromNode == v and G.edges[eid].toNode == u)) for eid in G.Adj[u]):
                        if not check_isect(G, u, v):
                            G.AddEdge(u, v)
                            en -= 1
                            
    if G.m < E: print(f"Warn: Reached {G.m}/{E} edges.")
    return G

def plot_g(G, W, H):
    plt.figure(figsize=(12, 6), facecolor='white') 
    for e in G.edges:
        plt.plot([G.nodes[e.fromNode].x, G.nodes[e.toNode].x], [G.nodes[e.fromNode].y, G.nodes[e.toNode].y], color='slategray', ls='-', lw=1.5, alpha=0.8, zorder=1)
        
    plt.scatter([n.x for n in G.nodes], [n.y for n in G.nodes], s=120, c='red', ec='darkred', lw=1.5, zorder=2)
    for n in G.nodes: plt.text(n.x, n.y + 4, str(n.id), fontsize=10, ha='center', va='bottom', fontweight='bold', zorder=3)
        
    plt.title("Planar Graph with Bridges", fontsize=16, fontweight='bold', pad=20)
    plt.xlabel("X Coord", fontsize=12, fontweight='bold')
    plt.ylabel("Y Coord", fontsize=12, fontweight='bold')
    plt.xlim(-W * 0.05, W * 1.05)
    plt.ylim(-H * 0.05, H * 1.15)
    
    ax = plt.gca()
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']: ax.spines[sp].set_linewidth(1.5)
    plt.tight_layout() 
    plt.show()

if __name__ == "__main__":
    G = gen_graph(10, 13, 3, 30, 50)
    if G:
        print(f"Success: {G.n} nodes, {G.m} edges.")
        plot_g(G, 30, 50)
    else:
        print("Failed.")