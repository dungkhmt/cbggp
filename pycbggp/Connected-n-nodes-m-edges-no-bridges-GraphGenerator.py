from CBSGG import Graph, Edge
from collections import deque

class DSU:
  def __init__(self, n):
    self.ncc = n
    self.parent = [i for i in range(n)]
    self.size = [1 for _ in range(n)]
    self.history = []

  def find(self, u):
    if self.parent[u] != u:
      return self.find(self.parent[u])
    return self.parent[u]

  def roll(self):
    if len(self.history) == 0:
      return
    u = self.history.pop()
    pu = self.parent[u]
    self.size[pu] -= self.size[u]
    self.parent[u] = u
    self.ncc += 1
  
  def union(self, u, v):
    pu = self.find(u)
    pv = self.find(v)
    if (pu == pv):
      return False
    if self.size[pu] < self.size[pv]:
      pu, pv = pv, pu
    self.parent[pv] = pu
    self.history.append(pv)
    self.size[pu] += self.size[pv]
    self.ncc -= 1
    return True

class ConnectedGraphWithoutBridgeGenerator:
  #generate randomly a connected graph with n nodes and m edges and no bridge
  # TODO by DungNT

  #generate all connected graphs with n nodes and m edges and no bridge
  def bfs(self, n):
    visited = 1
    queue = deque([0])
    while (len(queue) > 0):
      u = queue.popleft()
      mask = ~visited & self.adj[u]
      while mask > 0:
        v = (mask & -mask).bit_length() - 1
        visited |= (1 << v)
        queue.append(v)
        mask &= mask - 1

    return visited == (1 << n) - 1

  def dfs(self, u, p, tin, low, timer = 0):
    tin[u] = low[u] = timer
    timer += 1
    mask = self.adj[u]
    while (mask > 0):
      v = (mask & -mask).bit_length() - 1
      mask &= mask - 1
      if v == p:
        continue
      if tin[v] != -1:
        low[u] = min(low[u], tin[v])
      else:
        if not self.dfs(v, u, tin, low, timer):
          return False
        low[u] = min(low[u], low[v])
        if low[v] > tin[u]:
          return False


    return True

  def Backtrack(self, n, m, G, idx):
    if (self.dsu.ncc > 1 and m < self.dsu.ncc):
      return 
    if (n * (n - 1) // 2 - idx < m):
      return 
    if (self.bfs(n) == False):
      return 
    if (m == 0):
      tin = [-1 for _ in range(n)]
      low = [-1 for _ in range(n)]
      if self.dfs(0, -1, tin, low):
        self.graphs.append(G.copy())
      return 
    
    x = self.edges[idx].toNode
    y = self.edges[idx].fromNode
    self.adj[x] ^= 1 << y
    self.adj[y] ^= 1 << x
    self.Backtrack(n, m, G, idx + 1)
    self.adj[x] ^= 1 << y
    self.adj[y] ^= 1 << x

    G.AddEdge1(self.edges[idx])
    e = self.dsu.union(x, y)

    self.Backtrack(n, m - 1, G, idx + 1)

    if (e):
      self.dsu.roll()
    G.PopEdge()

    return None

  def Generate(self, n, m):
    if (m > n * (n - 1) // 2 or m < n):
      return []

    G = Graph(n)
    self.dsu = DSU(n)
    self.graphs = []
    self.adj = [(1 << n) - 1 for _ in range(n)]
    self.edges = []
    for i in range(n):
      for j in range(i + 1, n):
        e = Edge(i, j, len(self.edges))
        self.edges.append(e)

    self.Backtrack(n, m, G, 0)

    return self.graphs

n = 4
m = 5
CG = ConnectedGraphWithoutBridgeGenerator()
Gs = CG.Generate(n, m)
for G in Gs:
  G.Print()
  print('------------')