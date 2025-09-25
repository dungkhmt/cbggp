from CBSGG import Graph
import random

'''
    TODO by Nguyen Minh Giap
'''

class DirectedAcycleGraphGenerator:
 #generate randomly a directed acycle graph with n nodes and m edges
 def Generate(self, nbNodes, nbEdges):
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

DAG = DirectedAcycleGraphGenerator()
G = DAG.Generate(5, 4)
if G != None:
    G.Print()
    G.PrintUnweighted()