import random
from CBSGG import Graph, Edge

class DSUSave:
 def __init__(self, u: int, sizeU: int, v: int, sizeV: int):
  self.u = u
  self.sizeU = sizeU
  self.v = v
  self.sizeV = sizeV

class DSURollBack:
 def __init__(self, n: int):
  self.n = n
  self.root = [i for i in range(n)]
  self.size = [1] * n
  self.history = []
  self.components = n

 def getRoot(self, u: int):
  if self.root[u] != u:
   self.root[u] = self.getRoot(self.root[u])
  return self.root[u]
 
 def union(self, u: int, v: int):
  u = self.getRoot(u)
  v = self.getRoot(v)
  if u == v:
   return False
  if self.size[u] < self.size[v]:
   u, v = v, u
  self.history.append(DSUSave(u, self.size[u], v, self.size[v]))
  self.root[v] = u
  self.size[u] += self.size[v]
  self.components -= 1
  return True

 def rollBack(self):
  if not self.history:
   return
  last = self.history.pop()
  self.root[last.u] = last.u
  self.size[last.u] = last.sizeU
  self.root[last.v] = last.v
  self.size[last.v] = last.sizeV
  self.components += 1

class ConnectedGraphGenerator:
 '''
    TODO by Nguyen Ngoc Tuan Anh
 '''
 
 # generate a connected graph with n nodes and m edges
 
 def dfs(self, u: int, par: int, graph: Graph):
  graph.node += 1
  graph.num[u] = graph.node
  graph.low[u] = graph.node
  res = 0
  for id in graph.Adj[u]:
   e = graph.edges[id]
   v = e.fromNode + e.toNode - u
   if v == par:
    continue
   if graph.num[v] == 0:
    res += self.dfs(v, u, graph)
    graph.low[u] = min(graph.low[u], graph.low[v])
    if graph.low[v] == graph.num[v]:
     res += 1
   else:
    graph.low[u] = min(graph.low[u], graph.num[v])
  return res

 def countBridges(self, graph: Graph):
  graph.node = 0
  graph.num = [0] * graph.n
  graph.low = [0] * graph.n
  return self.dfs(0, -1, graph)

 def BnB(self, idx: int, nbEdges: int):
  if (len(self.listEdges) - idx < max(nbEdges, self.dsu.components - 1)):
   return
  if (idx == len(self.listEdges)):
   if (nbEdges == 0 and self.dsu.components == 1):
    numberBridges = self.countBridges(self.graph)
    if (numberBridges == self.nbBridges):
     self.listGraphs.append(self.graph.copy())
   return
  self.BnB(idx + 1, nbEdges)
  if (nbEdges > 0):
   self.graph.AddEdge(self.listEdges[idx][0], self.listEdges[idx][1])
   ok = self.dsu.union(self.listEdges[idx][0], self.listEdges[idx][1])
   self.BnB(idx + 1, nbEdges - 1)
   self.graph.PopEdge()
   if ok:
    self.dsu.rollBack()
  return

 def GenerateAll(self, nbNodes: int, nbEdges: int, nbBridges: int):
  # return a list of all undirected connected graphs with nbNodes nodes, nbEdges edges, and nbBridges bridges
  self.nbBridges = nbBridges
  self.graph = Graph(nbNodes)
  self.dsu = DSURollBack(nbNodes)
  self.listEdges = [[u, v] for u in range(nbNodes) for v in range(u + 1, nbNodes)]
  self.listGraphs = []
  self.BnB(0, nbEdges)
  print(len(self.listGraphs))
  return self.listGraphs

 def Generate(self, nbNodes, nbEdges, nbBridges):
  # return a undirected connected graph with nbNodes nodes, nbEdges edges, and nbBridges bridges
  listGraphs = self.GenerateAll(nbNodes, nbEdges, nbBridges)
  # Nếu tồn tại thì trả một đồ thị ngẫu nhiên trong danh sách này
  if listGraphs:
   return listGraphs[random.randint(0, len(listGraphs) - 1)]
  return None

 def SaveToFile(self, filename):
  return
 
gen = ConnectedGraphGenerator()
g = gen.Generate(7, 12, 2)
if g:
 g.Print()
