import random
from CBSGG import Graph,DirectedGraph
from Geometry import Point
from Delaunay import Delaunay
from DSU import DSU

def generate_connected_planar_graph(nb_nodes,nb_edges):
    if nb_nodes <= 0:
        return None

    max_edges = 0 if nb_nodes == 1 else (1 if nb_nodes == 2 else 3 * nb_nodes - 6)
    if nb_edges < nb_nodes - 1 or nb_edges > max_edges:
        return None

    G = Graph(nb_nodes)
    if nb_nodes == 1:
        return G
    if nb_nodes == 2:
        G.AddEdge(0, 1)
        return G

    while True:
        bound = nb_nodes * nb_nodes * 100 + 1000
        x = [(0, 0), (bound, 0), (0, bound)]
        used = set(x)
        while len(x) < nb_nodes:
            px = random.randint(1, bound - 2)
            py = random.randint(1, bound - px - 1)
            if (px, py) not in used:
                used.add((px, py))
                x.append((px, py))

        points = [Point(a, b) for a, b in x]
        p_id = {(p.x, p.y): i for i, p in enumerate(points)}
        tri = Delaunay().triangulate(points)

        cand = set()
        for i in range(0, len(tri), 3):
            if i + 2 >= len(tri):
                break
            id = [
                p_id[(tri[i].x, tri[i].y)],
                p_id[(tri[i + 1].x, tri[i + 1].y)],
                p_id[(tri[i + 2].x, tri[i + 2].y)]
            ]
            for u, v in ((id[0], id[1]), (id[1], id[2]), (id[2], id[0])):
                if u != v:
                    cand.add((u, v) if u < v else (v, u))

        if len(cand) < nb_edges:
            continue

        edges = list(cand)
        random.shuffle(edges)
        dsu = DSU(nb_nodes)
        selected = []
        add = []
        for u, v in edges:
            if dsu.union(u, v):
                selected.append((u, v))
            else:
                add.append((u, v))

        if len(selected) != nb_nodes - 1:
            continue

        random.shuffle(add)
        selected.extend(add[:nb_edges - len(selected)])
        break

    for u, v in selected:
        G.AddEdge(u, v)

    return G 
