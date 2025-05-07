import sys 
import random as rd 

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
class Node:
 def __init__(self,id):
  self.id = id 
 def toStr(self):
  return '' + str(self.id)
  
class Edge:
 def __init__(self, u, v):
  self.u = u 
  self.v = v 
 
 def toStr(self):
  return '(' + self.fromNode.toStr() + ',' + self.toNode.toStr() + ')' 
  
class Graph:
 def __init__(self, nodes, sizeOfAdj):
  self.numberNodes = len(nodes)
  self.nodes = nodes 
  self.Adj = [[] for _ in range(sizeOfAdj)] #adjacent edges
  self.edges = [] #edges of the graph
  self.numberEdges = 0

 def AddEdge(self, newEdge):
  id = self.numberEdges
  self.Adj[newEdge.u].append(id)
  self.Adj[newEdge.v].append(id)
  self.edges.append(newEdge)
  self.numberEdges += 1
  
 def Print(self):
  print('Nodes = ',end = ' ')
  for v in self.nodes:
   print(str(v), end = ' ')  
  print('')
  for u in range(len(self.Adj)):
    if len(self.Adj[u]) == 0:
      continue
    print('Adj[' + str(u) + ']: ', end = ' ')
    for id in self.Adj[u]:
      v = self.edges[id].u + self.edges[id].v - u
      print(v, end = ' ')
    print('') 
   
   
class SubGraphGenerator:
 def __init__(self,G):
  self.G = G 
 def GenGraph(self, nbNodes, nbEdges):
  # generate all connected graph with nbNodes and nbEdges with m <= 20
  listGraph = []
  for mask in range(1, 1 << self.G.numberEdges):
    # check if the graph is connected
    # check if the number of nodes is nbNodes and number of edges is nbEdges
    if int.bit_count(mask) != nbEdges:
      continue
    visited = [False] * self.G.numberNodes
    listNodesCurrent = []
    for id in range(self.G.numberEdges):
      if mask & (1 << id):
        if not visited[self.G.edges[id].u]:
          listNodesCurrent.append(self.G.edges[id].u)
          visited[self.G.edges[id].u] = True
        if not visited[self.G.edges[id].v]:
          listNodesCurrent.append(self.G.edges[id].v)
          visited[self.G.edges[id].v] = True
    if len(listNodesCurrent) != nbNodes:
      continue
    # check if the graph is connected
    visited = [False] * self.G.numberNodes
    checkConnected = True
    for u in listNodesCurrent:
      if not visited[u]:
        ans = self.dfs(u, mask, visited)
        if ans != len(listNodesCurrent):
          checkConnected = False
          break
    if not checkConnected:
      continue
    # add the graph to the list
    GMask = Graph(listNodesCurrent, self.G.numberNodes)
    for id in range(self.G.numberEdges):
      if mask & (1 << id):
        GMask.AddEdge(self.G.edges[id])
    listGraph.append(GMask)
  return listGraph
 
 def dfs(self, u, mask, visited):
  visited[u] = True
  ans = 1
  for id in self.G.Adj[u]:
    if not (mask & (1 << id)):
      continue
    v = self.G.edges[id].u + self.G.edges[id].v - u
    if not visited[v]:
      ans += self.dfs(v, mask, visited)
  return ans
 
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
 def __init__(self,vt):
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
 
 def __init__(self,G,r):
  self.G = G
  self.r = r 
  
class VarPath:
 '''
 object representing a dynamic path from s to t on G 
 path can be change with CHANGE operators (replaceEdge) on VarRootedTree(G,t)
 ''' 
 def __init__(self,s,t,G):
  self.G = G 
  self.s = s 
  self.t = t 
  
# main test...   
n = 5   
nodes = []
for i in range(n):
 nodes.append(i)

G = Graph(nodes, n)
G.AddEdge(Edge(0, 1))
G.AddEdge(Edge(0, 2))
G.AddEdge(Edge(1, 3))
G.AddEdge(Edge(2, 4))

GraphGenerator = SubGraphGenerator(G)
listGraph = GraphGenerator.GenGraph(3, 2)
# print all subgraphs generated
for i in range(len(listGraph)):
  g = listGraph[i]
  g.Print()
  print('')