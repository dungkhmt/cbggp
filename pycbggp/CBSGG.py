import sys
import random as rd
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Union

'''
Problems: generate a subgraph/subtree of a given graph with some properties:
  - connected
  - number of nodes is A and number of edges is B 
  - degree of each node is between L and U 
  - A specified node u has degree h 
  - Longest simple path is <= L 
  - contains A bridges (create A+1 subgraph without no bridge, and then reconnecting these subgraphs)
  - contains B articulation points 
  - contains NO bridges
  etc.
'''
# class Node:
#  def __init__(self,id):
#   self.id = id
#  def toStr(self):
#   return '' + str(self.id)


class Edge:
  def __init__(self, fromNode: int, toNode: int, id: int, weight: int = 0):
    self.fromNode = fromNode
    self.toNode = toNode
    self.weight = weight
    self.id = id

  def toStr(self):
    return '(' + str(self.fromNode) + ',' + str(self.toNode) + ')'


class Graph:
  '''
    TO BE DONE by Nguyen Tan Dung
    maintain: a list of edges E 
                A[v] is the list of indices of adjacent (outgoing) edges 
  '''
  # def __init__(self, nodes):
  #  self.nodes = nodes
  #  self.Adj = {} #adjacent edges
  #  for v in nodes:
  #   self.Adj[v] = []

  def __init__(self, n: int):
    self.n = n
    self.edges: list[Edge] = []
    self.Adj: list[list[int]] = [[] for _ in range(n)]
    self.m = 0

  def AddEdge(self, fromNode: int, toNode: int, weight: int = 0):
    e = Edge(fromNode, toNode, self.m, weight)
    self.Adj[fromNode].append(self.m)
    self.Adj[toNode].append(self.m)
    self.edges.append(e)
    self.m += 1

  def AddEdge1(self, e: Edge):
    self.Adj[e.fromNode].append(self.m)
    self.Adj[e.toNode].append(self.m)
    self.edges.append(e)
    self.m += 1
  
  def addFromAdj(self, adj : list[list[int]]):
    for u in range(self.n):
      for v in adj[u]:
        if u < v:
          self.AddEdge(u,v)

  def PopEdge(self):
    if self.m == 0:
      return None
    e = self.edges.pop()
    self.Adj[e.fromNode].pop()
    self.Adj[e.toNode].pop()
    self.m -= 1
    return e

  def CopyGraphMapNewNodes(self, nodeIds):
    # by DungNT
    # map id of the current graph to new nodeIds 
    return None 
  
  def CopyGraphMapNewNodes(self, fromId, toId):
    # by DungNT
    # map id of the current graph to new node ids in the range [fromId..toId]  
    return None 
    
  def LinkTwoGraphs(self, G1, node1, G2, node2):
    # by DungNT 
    # return a new created Graph by adding an edge connecting node1 of G1 and node2 of G2 
    return None 
    
    
    
  def Print(self):
    print('Nodes = ', end=' ')
    for v in range(self.n):
      print(v, end=' ')
    # for v in self.nodes:
    #  print(v.toStr(),end = ' ')
    print('')
    # for n in self.nodes:
    #  print('Adj[' + n.toStr() + ']: ',end = ' ')

    print('Edges = ', end=' ')
    for e in self.edges:
      print(e.toStr(), end=' ')
    print('')

    for n in range(self.n):
      print('Adj[' + str(n) + ']: ', end=' ')

      for e in self.Adj[n]:
        print(self.edges[e].toStr(), end=' ')
      print('')

  def PrintWeighted(self, offset = 0):
    print(self.n, self.m)
    for e in self.edges:
      print(e.fromNode + offset, e.toNode + offset, e.weight)

  def PrintUnweighted(self, offset = 0):
    print(self.n, self.m)
    for e in self.edges:
      print(e.fromNode + offset, e.toNode + offset)

  def PrintNodesEdgesFormat(self):
    V = '('
    for i in range(self.n):
      V = V + str(i)
      if i < self.n-1:
        V = V + ','
    V = V + ')'
    print(V)
    E = '{'
    for i in range(len(self.edges)):
      e = self.edges[i]
      E = E + '(' + str(e.fromNode) + ',' + str(e.toNode) + ')'
      if i < len(self.edges)-1:
        E = E + ','
    E = E + '}'
    print(E)    
  def SaveToFile(self, filename):
    with open(filename, 'w') as f:
      #print('Nodes = ', end=' ')
      line = ''
      for v in range(self.n):
        #print(v, end=' ')
        line = line + str(v) + ' '
      f.write(line + '\n')  
      # for v in self.nodes:
      #  print(v.toStr(),end = ' ')
      #print('')
      # for n in self.nodes:
      #  print('Adj[' + n.toStr() + ']: ',end = ' ')

      #print('Edges = ', end=' ')
      
      for e in self.edges:
       line = str(e.fromNode) + ' ' + str(e.toNode)
       #print(e.toStr(), end=' ')
       f.write(line + '\n')
    f.close()
    return

  def copy(self):
    G = Graph(self.n)
    for e in self.edges:
      G.AddEdge1(e)
    return G


class DirectedGraph(Graph):
  def __init__(self, n: int):
    super().__init__(n)

  def AddEdge(self, fromNode: int, toNode: int, weight: int = 0):
    e = Edge(fromNode, toNode, self.m, weight)
    self.Adj[fromNode].append(self.m)
    self.edges.append(e)
    self.m += 1

  def AddEdge1(self, e: Edge):
    self.Adj[e.fromNode].append(self.m)
    self.edges.append(e)
    self.m += 1

  def PopEdge(self):
    if self.m == 0:
      return None
    e = self.edges.pop()
    self.Adj[e.fromNode].pop()
    self.m -= 1
    return e

  def Print(self):
    print('Nodes = ', end=' ')
    for v in range(self.n):
      print(v, end=' ')
    print('')

    print('Edges = ', end=' ')
    for e in self.edges:
      print(e.toStr(), end=' ')
    print('')

    for n in range(self.n):
      print('Adj[' + str(n) + ']: ', end=' ')

      for e in self.Adj[n]:
        print(self.edges[e].toStr(), end=' ')
      print('')

  def PrintWeighted(self, offset = 0):
    print(self.n, self.m)
    for e in self.edges:
      print(e.fromNode + offset, e.toNode + offset, e.weight)

  def PrintUnweighted(self, offset = 0):
    print(self.n, self.m)
    for e in self.edges:
      print(e.fromNode + offset, e.toNode + offset)

  def addFromAdj(self, adj : list[list[int]]):
    for u in range(self.n):
      for v in adj[u]:
        self.AddEdge(u,v)

  def copy(self):
    G = DirectedGraph(self.n)
    for e in self.edges:
      G.AddEdge1(e)
    return G


class SubGraphGenerator:
  def __init__(self, G):
    self.G = G

  def GenGraph(self, nbNodes, nbEdges):
    return


class VarTree:
  '''
  object representing a subtree of a given graph
  provide CHANGE operators for generating neighboring trees: add, remove, replace edges 
  '''

  def __init__(self, G):
    self.G = G

  def AddEdge(self, e):
    return

  def RemoveEdge(self, e):
    return

  def ReplaceEdge(self, eIn, eOut):
    # replace edge eOut of current tree by another edge eIn
    return


class Weight:
  '''
  object representing the weight of a VarTree vt 
  '''

  def __init__(self, vt):
    self.vt = vt


class Diameter:
  '''
  object representing the longest path on the VarTree vt  
  '''

  def __init__(self, vt):
    self.vt = vt

  def Value(self):
    return 0

  def GetAddEdgeDelta(self, eIn):
    return 0

  def PropagateAddEdge(self, eIn):
    return

  def GetRemoveEdgeDelta(self, eOut):
    return 0

  def PropagateRemoveEdge(self, eOut):
    return

  def GetReplaceEdgeDelta(self, eIn, eOut):
    return 0

  def PropagateReplaceEdge(self, eIn, eOut):
    return


class VarRootedTree:
  '''
  object representing a subtree of G rooted at r 
  can be changed by applying AddEdge, RemoveEdge, ReplaceEdge  
  '''

  def __init__(self, G, r):
    self.G = G
    self.r = r


class VarPath:
  '''
  object representing a dynamic path from s to t on G 
  path can be change with CHANGE operators (replaceEdge) on VarRootedTree(G,t)
  '''

  def __init__(self, s, t, G):
    self.G = G
    self.s = s
    self.t = t

if __name__ == "__main__":
 connected = 'Y'
 nbNodes = 0
 nbEdges = 0
 nbBridges = 0
 maxDegree = 10
 gType = 'planar'
 directed = 'N' 
 filename = 'g.txt'
 for i in range(len(sys.argv)):
  #print('argv[',i,'] = ',sys.argv[i]) 
  if sys.argv[i] == '-connected':
   conneceted = sys.argv[i+1]
  elif sys.argv[i] == '-nbNodes':
   nbNodes = int(sys.argv[i+1])
  elif sys.argv[i] == '-nbEdges':
   nbEdges = int(sys.argv[i+1])
  elif sys.argv[i] == '-type':
   gType = sys.argv[i+1]
  elif sys.argv[i] == '-directed':
   directed = sys.argv[i+1] 
  elif sys.argv[i] == '-nbBridges':
   nbBridges = int(sys.argv[i+1])
  elif sys.argv[i] == '-maxDegree':
   maxDegree = int(sys.argv[i+1]) 
  elif sys.argv[i] == '-filename':
   filename = int(sys.argv[i+1])
   
 if gType == 'undirected':
  if nbBridges == 0 and connected == 'Y':
   CG = ConnectedGraphWithoutBridgeGenerator()
   G = CG.GenerateONE(nbNodes, nbEdges)
   G.SaveToFile(filename)

 elif gType == 'planar':
  PGG = PlanarGraphGenerator()
  G = PGG.Generate(nbNodes, nbEdges, nbBridges,maxDegree)
  G.Print()
  G.SaveToFile(filename)

# main test...
# n = 5
# # nodes = []
# # for i in range(n):
# #  nod = Node(i)
# #  nodes.append(nod)

# # G = Graph(nodes)

# G = Graph(n)
# # G.AddEdge(Edge(nodes[0],nodes[1]))
# # G.AddEdge(Edge(nodes[0],nodes[2]))
# # G.AddEdge(Edge(nodes[1],nodes[3]))
# # G.AddEdge(Edge(nodes[2],nodes[4]))

# G.AddEdge(0,1)
# G.AddEdge(0,2)
# G.AddEdge(1,3)
# G.AddEdge(2,4)

# G.Print()

# G1 = DirectedGraph(n)
# G1.AddEdge(0,1)
# G1.AddEdge(0,2)
# G1.AddEdge(1,3)
# G1.AddEdge(2,4)
# G1.Print()
