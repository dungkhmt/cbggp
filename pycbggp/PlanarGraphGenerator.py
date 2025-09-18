from CBSGG import Graph
from CBSGG import Node 
from CBSGG import Edge
'''
    TODO by Nguyen Tan Dung
    
''' 
class PlanarGraphGenerator:
 # generate planar graph with constraints
 
 def Generate(self,nbNodes,nbEdges, nbBridges, maxDegree):
  # return a random connected planar graph with nbNodes nodes and nbEdges edges and nbBridges bridges
  # degree of each node is <= maxDegree 
  
  nodes = []
  for i in range(nbNodes):
   nod = Node(i)
   nodes.append(nod)

  G = Graph(nodes)

  return G
  
 def SaveToFile(self, filename):
  return  

  
PG = PlanarGraphGenerator()
G = PG.Generate(5,9)
G.Print()  