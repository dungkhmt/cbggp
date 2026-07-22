from CBSGG import Graph,DirectedGraph
from Geometry import Point
from Delaunay import Delaunay
from DSU import DSU
import random

from heuristicmethods.Undirected_Tree_diameter_between_P_and_Q_and_degree_at_most_D import UndirectedTreeDiameterBetweenPandQAndDegreeAtMostDGenerator

from heuristicmethods.undirectedconnectedgraphkbridges import generate_undirected_connected_graph_nb_bridges
from heuristicmethods.directedgraphkstronglyconnectedcomponent import generate_directed_graph_nb_strongly_connected_components

from heuristicmethods.planargraph import generate_connected_planar_graph
from heuristicmethods.constrained_structural_graph_generator import generate_graph
from heuristicmethods.undirectedcompletegraph import generate_undirected_complete_graph
from heuristicmethods.bipartitegraph import generate_bipartite_graph
from heuristicmethods.connectedbipartitegraph import generate_connected_bipartite_graph
from heuristicmethods.directedgraph import generate_directed_graph
from heuristicmethods.directedstronglyconnectedgraph import generate_directed_strongly_connected_graph

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
 G = generate_undirected_complete_graph(nb_nodes)
 return G  

def gen_bipartite_graph(nb_left_nodes, nb_right_nodes,nb_edges):
 # TODO by QuyetLG
 G = generate_bipartite_graph(nb_left_nodes, nb_right_nodes,nb_edges)
 return G  

def gen_connected_bipartite_graph(nb_left_nodes, nb_right_nodes,nb_edges):
 # TODO by QuyetLG
 G = generate_connected_bipartite_graph(nb_left_nodes, nb_right_nodes, nb_edges)
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
 G = generate_directed_graph(nb_nodes, nb_edges)
 return G  

def gen_directed_strongly_connected_graph(nb_nodes, nb_edges):
 # TODO by QuyetLG
 G = generate_directed_strongly_connected_graph(nb_nodes, nb_edges)
 return G 
 
def gen_directed_graph_nb_strongly_connected_components(nb_nodes, nb_edges, nb_strongly_connected_components):
    return generate_directed_graph_nb_strongly_connected_components(nb_nodes,nb_edges, nb_strongly_connected_components)

# planar graphs
def gen_connected_planar_graph(nb_nodes,nb_edges):
    return generate_connected_planar_graph(nb_nodes, nb_edges)
    
   
def test1():    
    G = gen_connected_planar_graph(10, 15)
    G.Print()  
    G.PrintUnweighted()  
    G.gen_random_weight_edges(1,10)
    G1 = G.CopyGraphMapNewNodes(1,10)
    G1.SaveToFileWeighted('1.txt')

def test2():
    G1 = gen_connected_planar_graph(4,5)
    G1.Print()
    

    G2 = gen_connected_planar_graph(4,5)
    G2.Print()
    
    G3 = G2.CopyGraphMapNewNodes(4,7)
    G3.Print()
    G = G1.union_two_graphs(G1,G3)
    G.Print()
    
#test1()
test2()
G = generate_directed_strongly_connected_graph(4, 10)
G.Print()    


# G = gen_undirected_connected_graph_nb_bridges(7, 8, 2)
# G.Print()

# G = gen_directed_graph_nb_strongly_connected_components(6, 9, 2)
# G.Print()

# G = gen_undirected_tree_diameter_between_P_and_Q_and_degree_at_most_D(6, 2, 3, 3)
# G.Print()
 
