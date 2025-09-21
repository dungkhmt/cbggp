import random
import math
from CBSGG import Graph

# from CBSGG import Node 
from CBSGG import Edge

'''
  TODO by Nguyen Tan Dung
''' 

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


class PlanarGraphGenerator:
  # generate planar graph with constraints

  def intersect(self, p1, p2, p3, p4):
    # check if line segment p1p2 and p3p4 intersect
    # collinear points are considered intersecting

    def onSegment(p, q, r):
      if (q[0] <= max(p[0], r[0]) and q[0] >= min(p[0], r[0]) and
          q[1] <= max(p[1], r[1]) and q[1] >= min(p[1], r[1])):
        return True
      return False
    def orientation(p, q, r):
      # find the orientation of the ordered triplet (p, q, r)
      # 0 -> p, q and r are collinear
      # 1 -> Clockwise
      # 2 -> Counterclockwise
      val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
      if val == 0:
        return 0
      elif val > 0:
        return 1
      else:
        return 2
    
    o1 = orientation(p1, p2, p3)  
    o2 = orientation(p1, p2, p4)
    o3 = orientation(p3, p4, p1)
    o4 = orientation(p3, p4, p2)
    if o1 != o2 and o3 != o4:
      return True
    if o1 == 0 and onSegment(p1, p3, p2):
      return True
    if o2 == 0 and onSegment(p1, p4, p2):
      return True
    if o3 == 0 and onSegment(p3, p1, p4):
      return True
    if o4 == 0 and onSegment(p3, p2, p4):
      return True

    return False

  def Generate(self, n, m, k, d):
    #return a random connected planar graph with n nodes and m edges and k bridges
    #degree of each node is <= d 
    #TODO

    if m < n - 1 or m > 3 * n - 6:
      return None
    if k < 0 or k > n - 1:
      return None
    if d < 0 or n * d < 2 * m:
      return None
    if d == 0:
      if n == 1 and m == 0:
        return Graph(1)
      else:
        return None
    
    # print('valid')

    G = Graph(n)
    points = [() for _ in range(n)]
    for i in range(n):
      points[i] = (random.randint(0, 1000000000), random.randint(0, 1000000000))

    candidates = set()
    for i in range(n):
      candidates.add(i)

    degrees = [0 for _ in range(n)]
    mat = [[False for _ in range(n)] for _ in range(n)]
    check = [[False for _ in range(n)] for _ in range(n)]
    cnt = 0
    # for i in range(m):
    edges = []
    while True:
      # for _ in range(10000):
        if (len(candidates) < 2):
          # return None
          break
        u = random.choice(list(candidates))
        ok = False
        candidates_u = set()
        for v in candidates:
          if v != u and mat[u][v] == False and check[u][v] == False:
            ok = True
            candidates_u.add(v)
            break

        if ok == False:
          candidates.remove(u)
          continue

        v = random.choice(list(candidates_u))

        intersected = False
        for e in edges:
          x = e[0]
          y = e[1]
          if x == u or x == v or y == u or y == v:
            continue
          if self.intersect(points[u], points[v], points[x], points[y]):
            intersected = True
            break
        
        if intersected == True:
          check[u][v] = True
          check[v][u] = True
          continue

        mat[u][v] = True
        mat[v][u] = True
        # G.AddEdge(u, v)
        edges.append((u, v))
        degrees[u] += 1
        degrees[v] += 1
        if degrees[u] == d:
          candidates.remove(u)
        if degrees[v] == d:
          candidates.remove(v)
        cnt += 1
        # break

    if cnt < m:
      return None
    
    random.shuffle(edges)
    dsu = DSU(n)
    for i in range(n):
      for j in range(n):
        mat[i][j] = mat[j][i] = False

    for e in edges:
      if dsu.union(e[0], e[1]):
        G.AddEdge(e[0], e[1])
        mat[e[0]][e[1]] = True
        mat[e[1]][e[0]] = True
        if dsu.ncc == 1:
          print('connected')
          break

    cnt = m - n + 1
    random.shuffle(edges)
    for e in edges:
      if cnt == 0:
        break
      if mat[e[0]][e[1]] == False:
        G.AddEdge(e[0], e[1])
        cnt -= 1

    return G


'''
PG = PlanarGraphGenerator()
G = PG.Generate(5,9)
G.Print()  
'''

PG = PlanarGraphGenerator()
G = PG.Generate(50, 100, 1, 6)
if G != None:
  G.Print()
  G.PrintUnweighted()