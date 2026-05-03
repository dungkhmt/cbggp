from CBSGG import Graph,DirectedGraph

# undirected graphs
def gen_undirected_graph(nb_nodes, nb_edges, nb_connected_components, nb_bridges, nb_articulation_points):
 G = Graph(nb_nodes)
 return G  

def gen_undirected_connected_graph(nb_nodes, nb_edges):
 G = Graph(nb_nodes)
 return G  
 
def gen_undirected_connected_graph_no_bridge(nb_nodes, nb_edges):
 G = Graph(nb_nodes)
 return G  

def gen_undirected_connected_graph_no_articulation_point(nb_nodes, nb_edges):
 G = Graph(nb_nodes)
 return G  

def gen_undirected_connected_graph_no_bridge_no_articulation_point(nb_nodes, nb_edges):
 G = Graph(nb_nodes)
 return G  

# special undirected graphs 
def gen_undirected_complete_graph(nb_nodes):
 G = Graph(nb_nodes)
 return G  

def gen_bipartite_graph(nb_left_nodes, nb_right_nodes,nb_edges):
 G = Graph(nb_nodes)
 return G  

def gen_connected_bipartite_graph(nb_left_nodes, nb_right_nodes,nb_edges):
 G = Graph(nb_nodes)
 return G  


# undirected trees
def gen_undirected_tree(nb_nodes):
 G = Graph(nb_nodes)
 return G  

def gen_undirected_tree_bounded_diameter_degree(nb_nodes, ub_deg, ub_diameter):
 G = Graph(nb_nodes)
 return G  
  
# directed graphs 
def gen_directed_graph(nb_nodes, nb_edges):
 G = DirectedGraph(nb_nodes)
 return G  

def gen_directed_strongly_connected_graph(nb_nodes, nb_edges):
 G = DirectedGraph(nb_nodes)
 return G 
 
def gen_directed_graph_nb_strongly_connected_components(nb_nodes, nb_edges, nb_strongly_connected_components):
 G = DirectedGraph(nb_nodes)
 return G
 
 
# main    
G = gen_undirected_connected_graph(5,10)
G.Print()

 
 