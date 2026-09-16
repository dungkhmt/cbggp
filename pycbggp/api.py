from CBSGG import Graph,DirectedGraph, Point2D
from Geometry import Point
from Delaunay import Delaunay
from DSU import DSU
import random

# Heuristic methods 
from heuristicmethods.Undirected_Tree_diameter_between_P_and_Q_and_degree_at_most_D import UndirectedTreeDiameterBetweenPandQAndDegreeAtMostDGenerator
from heuristicmethods.undirectedconnectedgraphkbridges import generate_undirected_connected_graph_nb_bridges
from heuristicmethods.directedgraphkstronglyconnectedcomponent import generate_directed_graph_nb_strongly_connected_components
from heuristicmethods.constrained_graph import generate_graph
from heuristicmethods.planargraph import generate_connected_planar_graph
from heuristicmethods.undirectedcompletegraph import generate_undirected_complete_graph
from heuristicmethods.bipartitegraph import generate_bipartite_graph
from heuristicmethods.connectedbipartitegraph import generate_connected_bipartite_graph
from heuristicmethods.directedgraph import generate_directed_graph
from heuristicmethods.directedstronglyconnectedgraph import generate_directed_strongly_connected_graph
from heuristicmethods.connected_undirected_planar_graph_with_bridges_generator import generate_connected_planar_graph_with_bridges

# Constructive methods 
from constructivemethods.undirected_connected_graph import generate_undirected_connected_graph
from constructivemethods.undirected_connected_graph_no_bridge import generate_undirected_connected_graph_no_bridge
from constructivemethods.biconnected_graph import generate_biconnected_graph
from constructivemethods.undirected_tree import generate_undirected_tree
from constructivemethods.undirected_tree_bounded_diameter_degree import generate_undirected_tree_bounded_diameter_degree
from constructivemethods.constrained_graph import (
    generate_constructive, tarjan_analysis, check_feasibility as check_feasibility_vecba,
    verify_graph,
)
from constructivemethods.connected_undirected_planar_generator import _gen_connected_planar_graph


def gen_connected_planar_graph_given_list_of_point(points: list[Point2D], nb_edges):
    """ Done by lmToT27!!!
        Step 1: Use Graham Scan to generate convex-hull layers, (allow collinear points on a convex hull).
            We can claim that theres at least V - 1 and at most 3V - k - 3 egdes of a connected planar graph (k is number of points on the outermost layer).
        Step 2: Add edges between 2 adjacent points on every layer.
        Step 3:
        Triangulate the space between two adjacent layers (Outer 'O' and Inner 'I').
            - Rotate the Outer layer to align its starting point with the Inner layer.
            - Sort points in both layers by polar angle with respect to the center point of the Inner layer.
            - We sweep using 2 pointers to determine the next edge to add, try to select closest pair of points from 2 layers.
            - Validate the chosen diagonal using the cross product to strictly prevent collinear or overlapping triangles.
        Step 4: In the innermost layer, triagulate by adding edges following a "zig-zag" pattern or something like so.
        Final step is easy we can use kruskal's algorithm to add edges until we reach the desired number of edges.
        Time complexity: O(|V|^2) for convex hull generation, O(|V|log|V|) for triangulation, O(|E|log|E|) for kruskal's algorithm.
        Space complexity: O(|V|) for convex hull generation, O(|V|) for triangulation, O(|E|) for kruskal's algorithm.
        It can be O(|V|log|V|) for convex hull generation but it doesn't seem like i can implement it TwT.
    """
    return _gen_connected_planar_graph(points, nb_edges)

def gen_connected_planar_graph_with_bridges(points: list[Point2D], nb_edges, nb_bridges):
    nb_nodes = len(points)
    # generate an undirected connected planar graph containing nb_nodes and nb_edges, nodes are located at points 
    # containing nb_bridges 
    # TODO by Nguyen Phuc Khanh 
    ''' Use Combination of heuristic algorithms 
        Heuristic01: greedy approach by generating a maximum planar graph, 
         finding MST, adding edges to get exactly expected number of bridges,
         finally adding edges for get exactly expected number of edges without changing number of bridges.
        Time complexity:  O(V^2)
        Works best for graphs with a moderate number of edges and a small number of bridges
        Heuristic02: greedy approach by using motone chain & spliting core-branch
        Works best for graphs with a small number of edges and a moderate number of bridges
    '''
    G = generate_connected_planar_graph_with_bridges(points, nb_edges, nb_bridges)
    return G 
    

    
# undirected graphs
def gen_undirected_graph(nb_nodes, nb_edges, nb_connected_components,
                        nb_bridges, nb_articulation_points,
                        method='constructive', **kwargs):
    # by Nguyen Ngoc Tuan Anh
    # Generate an undirected graph satisfying exact constraints (V, E, C, B, A).
    # method:
    #     'constructive' — 4-phase constructive algorithm, O(V+E), 100% exact (default).
    #     'mcmc'         — MCMC edge rewiring, produces diverse random graphs.
    #                      kwargs: nb_iterations=1000
    #     'sa'           — Simulated Annealing, heuristic fallback.
    #                      kwargs: max_iterations=10000, T_init=100.0, alpha=0.995

    # Returns Graph or None if the input is infeasible.
    return generate_graph(nb_nodes, nb_edges, nb_connected_components,
                          nb_bridges, nb_articulation_points,
                          method=method, **kwargs)

def gen_undirected_connected_graph(nb_nodes, nb_edges):
    # by Nguyen Ngoc Tuan Anh
    # Algorithm: Random Spanning Tree (DSU) + Random Edge Fill. O(V+E).
    return generate_undirected_connected_graph(nb_nodes, nb_edges)
 
def gen_undirected_connected_graph_no_bridge(nb_nodes, nb_edges):
    # by Nguyen Ngoc Tuan Anh
    # Algorithm: Hamiltonian Cycle + Random Extra Edges. O(V+E).
    return generate_undirected_connected_graph_no_bridge(nb_nodes, nb_edges)

def gen_undirected_connected_graph_no_articulation_point(nb_nodes, nb_edges):
    # by Nguyen Ngoc Tuan Anh
    # Algorithm: Open Ear Decomposition. O(V+E).
    # Biconnected graph ⟹ no articulation point AND no bridge.
    return generate_biconnected_graph(nb_nodes, nb_edges)

def gen_undirected_connected_graph_no_bridge_no_articulation_point(nb_nodes, nb_edges):
    # by Nguyen Ngoc Tuan Anh
    # Algorithm: Open Ear Decomposition. O(V+E).
    # Biconnected ⟺ no articulation point ⟹ no bridge (same as no-AP task).
    return generate_biconnected_graph(nb_nodes, nb_edges)
 
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
    # by Nguyen Ngoc Tuan Anh
    # Algorithm: Prüfer Sequence Decode — uniform distribution over V^(V-2) trees.
    # O(V log V).
    return generate_undirected_tree(nb_nodes)

def gen_undirected_tree_bounded_diameter_degree(nb_nodes, ub_deg, ub_diameter):
    # by Nguyen Ngoc Tuan Anh
    # Algorithm: Center-Rooted BFS Layer Growth. O(V).
    return generate_undirected_tree_bounded_diameter_degree(nb_nodes, ub_deg, ub_diameter)

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
 
