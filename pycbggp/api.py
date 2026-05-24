from CBSGG import Graph,DirectedGraph
import random
from Undirected_Tree_diameter_between_P_and_Q_and_degree_at_most_D import UndirectedTreeDiameterBetweenPandQAndDegreeAtMostDGenerator

from heuristicmethods.undirectedconnectedgraphkbridges import generate_undirected_connected_graph_nb_bridges
from heuristicmethods.directedgraphkstronglyconnectedcomponent import generate_directed_graph_nb_strongly_connected_components



# undirected graphs
def gen_undirected_graph(nb_nodes, nb_edges, nb_connected_components, nb_bridges, nb_articulation_points):
 # TODO by Nguyen Ngoc Tuan Anh
 G = Graph(nb_nodes)
 return G  

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
    #TODO by DungNT
    G = Graph(nb_nodes)
    return None 
    
    
    

G = gen_undirected_connected_graph_nb_bridges(7, 8, 2)
G.Print()

G = gen_directed_graph_nb_strongly_connected_components(6, 9, 2)
G.Print()

# G = gen_undirected_tree_diameter_between_P_and_Q_and_degree_at_most_D(6, 2, 3, 3)
# G.Print()
 