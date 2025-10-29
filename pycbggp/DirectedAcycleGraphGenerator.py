from CBSGG import Graph
import random

'''
    TODO by Nguyen Minh Giap
'''

class DirectedAcycleGraphGenerator:
  def generate(self, nbNodes, nbEdges):
    """
    Generate a random Directed Acycle Graph with nbNodes nodes, nbEdges edges
  
    Args:
      nbNodes: The number of nodes in the graph
      nbEdges: The number of edges in the graph
    
    Returns:
      A directed acycle graph if successful, None otherwise
    """
    
    if nbNodes < 0:
      print("Invalid number of nodes.")
      return None
    if nbEdges < 0 or nbEdges > nbNodes * (nbNodes - 1) // 2:
      print("Invalid number of edges.")
      return None
  
    G = Graph(nbNodes)
    topo_order = list(range(nbNodes))
    random.shuffle(topo_order)
    edges = set()
  
    while len(edges) < nbEdges:
        i, j = random.sample(range(nbNodes), 2)
        if i > j:
            i, j = j, i
        u, v = topo_order[i], topo_order[j]
        if (u, v) not in edges:
          edges.add((u, v))
          G.AddEdge(u, v)

    return G
  
  def generate_by_total_degree(self, nbNodes, nbEdges, max_degree):
    """
    Generate a random Directed Acycle Graph with nbNodes nodes, nbEdges edges,
    the maximum degree of each node
    
    Args:
      nbNodes: The number of nodes in the graph
      nbEdges: The number of edges in the graph
      max_degree: The maximum degree of each node
      
    Returns:
      A directed acycle graph if successful, None otherwise
    """
    
    if nbNodes < 0:
      print("Invalid number of nodes.")
      return None
    if nbEdges < 0 or nbEdges > nbNodes * (nbNodes - 1) // 2:
      print("Invalid number of edges.")
      return None
    
    G = Graph(nbNodes)
    topo_order = list(range(nbNodes))
    random.shuffle(topo_order)
    edges = set()
    in_degree = [0] * nbNodes
    out_degree = [0] * nbNodes
    
    while len(edges) < nbEdges:
      i,j = random.sample(range(nbNodes), 2)
      if i > j:
        i, j = j, i
      u, v = topo_order[i], topo_order[j]
      if (in_degree[u] + out_degree[u]) >= max_degree:
        continue
      if (in_degree[v] + out_degree[v]) >= max_degree:
        continue
      if (u, v) not in edges:
        edges.add((u, v))
        in_degree[v] += 1
        out_degree[u] += 1
        G.AddEdge(u, v)
        
    return G 

  def generate_by_io_degree(self, nbNodes, nbEdges, max_in_degree, max_out_degree):
    """
    Generate a random Directed Acycle Graph with nbNodes nodes, nbEdges edges, 
    the maximum in-degree of each node, the maximum out-degree of each node is
    
    Args:
      nbNodes: The number of nodes in the graph
      nbEdges: The number of edges in the graph
      max_in_degree: The maximum in-degree of each node
      max_out_degree: The maximum out-degree of each node
    
    Returns:
      A directed acycle graph if successful, None otherwise
    """
    if nbNodes < 0:
      print("Invalid number of nodes.")
      return None
    if nbEdges < 0 or nbEdges > nbNodes * (nbNodes - 1) // 2:
      print("Invalid number of edges.")
      return None
    
    G = Graph(nbNodes)
    topo_order = list(range(nbNodes))
    random.shuffle(topo_order)
    edges = set()
    in_degree = [0] * nbNodes
    out_degree = [0] * nbNodes
    
    while len(edges) < nbEdges:
      i,j = random.sample(range(nbNodes), 2)
      if i > j:
        i, j = j, i
      u, v = topo_order[i], topo_order[j]
      if (in_degree[v]) >= max_in_degree:
        continue
      if (out_degree[u]) >= max_out_degree:
        continue
      if (u, v) not in edges:
        edges.add((u, v))
        in_degree[v] += 1
        out_degree[u] += 1
        G.AddEdge(u, v)
        
    return G 
    
DAG = DirectedAcycleGraphGenerator()
G = DAG.generate(5, 4)
if G != None:
    G.Print()
    G.PrintUnweighted()
else:
  print("None")
    
G1 = DAG.generate_by_total_degree(10, 8, 2)
if G1 != None:
  G1.Print()
  G1.PrintUnweighted()
else:
  print("None")
  
G2 = DAG.generate_by_io_degree(100, 88, 2, 2)
if G2 != None:
  G2.Print()
  G2.PrintUnweighted()
else:
  print("None")
  