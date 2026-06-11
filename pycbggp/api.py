from CBSGG import Graph,DirectedGraph
from Geometry import Point
from Delaunay import Delaunay
from DSU import DSU
import random

from heuristicmethods.Undirected_Tree_diameter_between_P_and_Q_and_degree_at_most_D import UndirectedTreeDiameterBetweenPandQAndDegreeAtMostDGenerator

from heuristicmethods.undirectedconnectedgraphkbridges import generate_undirected_connected_graph_nb_bridges
from heuristicmethods.directedgraphkstronglyconnectedcomponent import generate_directed_graph_nb_strongly_connected_components
from heuristicmethods.constrained_structural_graph_generator import generate_graph



# undirected graphs
def gen_undirected_graph(nb_nodes, nb_edges, nb_connected_components,
                        nb_bridges, nb_articulation_points,
                        method='constructive', **kwargs):
    """
    Generate an undirected graph satisfying exact constraints (V, E, C, B, A).

    method:
        'constructive' — Direct construction, O(V+E), 100% exact.
        'mcmc'         — MCMC Edge Rewiring, generates diverse random graphs.
                         kwargs: nb_iterations=1000
        'sa'           — Simulated Annealing, heuristic for complex cases.
                         kwargs: max_iterations=10000, T_init=100.0, alpha=0.995

    Returns Graph or None. See constrained_structural_graph_generator.py for details.

    Examples:
        G = gen_undirected_graph(10, 15, 1, 2, 3)
        G = gen_undirected_graph(10, 15, 1, 2, 3, method='mcmc', nb_iterations=2000)
        G = gen_undirected_graph(8, 12, 1, 0, 0, method='sa', max_iterations=5000)
    """
    return generate_graph(nb_nodes, nb_edges, nb_connected_components,
                          nb_bridges, nb_articulation_points,
                          method=method, **kwargs)

def gen_undirected_connected_graph(nb_nodes, nb_edges):
 # TODO by Nguyen Ngoc Tuan Anh
 
 G = Graph(nb_nodes)
 
 return G  
 
def gen_undirected_connected_graph_no_bridge(nb_nodes, nb_edges):
 # TODO by Nguyen Ngoc Tuan Anh 
 G = Graph(nb_nodes)
 
 return G  

def gen_undirected_connected_graph_no_articulation_point(nb_nodes, nb_edges):
 # TODO by Nguyen Ngoc Tuan Anh
 
 G = Graph(nb_nodes)
 return G  

def gen_undirected_connected_graph_no_bridge_no_articulation_point(nb_nodes, nb_edges):
 # TODO by Nguyen Ngoc Tuan Anh
 
 G = Graph(nb_nodes)
 return G  
 
def gen_undirected_connected_graph_nb_bridges(nb_nodes, nb_edges, nb_bridges):
    return generate_undirected_connected_graph_nb_bridges(nb_nodes, nb_edges, nb_bridges)



# special undirected graphs 
def gen_undirected_complete_graph(nb_nodes):
 # TODO by QuyetLG
 G = Graph(nb_nodes)
 return G  

def gen_bipartite_graph(nb_left_nodes, nb_right_nodes,nb_edges):
 # TODO by QuyetLG
 G = Graph(nb_nodes)
 return G  

def gen_connected_bipartite_graph(nb_left_nodes, nb_right_nodes,nb_edges):
 # TODO by QuyetLG
 G = Graph(nb_nodes)
 return G  


# undirected trees
def gen_undirected_tree(nb_nodes):
 # TODO by Nguyen Ngoc Tuan Anh    
 G = Graph(nb_nodes)
 return G  

def gen_undirected_tree_bounded_diameter_degree(nb_nodes, ub_deg, ub_diameter):
 # TODO by Nguyen Ngoc Tuan Anh 
 G = Graph(nb_nodes)
 return G  

def gen_undirected_tree_diameter_between_P_and_Q_and_degree_at_most_D(n, p, q, d):
    gen = UndirectedTreeDiameterBetweenPandQAndDegreeAtMostDGenerator()
    G = gen.generate(n, p, q, d)
    return G
  
# directed graphs 
def gen_directed_graph(nb_nodes, nb_edges):
 # TODO by QuyetLG
 G = DirectedGraph(nb_nodes)
 return G  

def gen_directed_strongly_connected_graph(nb_nodes, nb_edges):
 # TODO by QuyetLG
 G = DirectedGraph(nb_nodes)
 return G 
 
def gen_directed_graph_nb_strongly_connected_components(nb_nodes, nb_edges, nb_strongly_connected_components):
    return generate_directed_graph_nb_strongly_connected_components(nb_nodes,nb_edges, nb_strongly_connected_components)

# planar graphs
def gen_connected_planar_graph(nb_nodes,nb_edges):
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
    
    
G = gen_connected_planar_graph(10, 15)
G.Print()    


# G = gen_undirected_connected_graph_nb_bridges(7, 8, 2)
# G.Print()

# G = gen_directed_graph_nb_strongly_connected_components(6, 9, 2)
# G.Print()

# G = gen_undirected_tree_diameter_between_P_and_Q_and_degree_at_most_D(6, 2, 3, 3)
# G.Print()
 
