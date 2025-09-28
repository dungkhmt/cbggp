import random
import math
from CBSGG import Graph
from functools import cmp_to_key

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
  
  def ccw(self, A, B, C):
    # val = (B[1] - A[1]) * (C[0] - B[0]) - (B[0] - A[0]) * (C[1] - B[1])
    val = (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0])
    if val == 0:
      return 0
    elif val > 0:
      return 1
    else:
      return -1
  
  def dfs(self, u, p, adj, tin, low, st, timer = 0):
    tin[u] = low[u] = timer
    timer += 1
    st.append(u)
    for e in adj[u]:
      v = e.toNode ^ e.fromNode ^ u
      if v == p:
        continue
      if tin[v] != -1:
        low[u] = min(low[u], tin[v])
      else:
        low[u] = min(low[u], low[v])
        if low[v] > tin[u]:
          None


  def Generate(self, n, m, k, d):
    # return a random connected planar graph with n nodes and m edges and k bridges
    # degree of each node is <= d 
    #TODO

    if m < n - 1 or m > 3 * n - 6:
      return None
    if k < 0 or k > n - 1:
      return None
    # if d < 0 or n * d < 2 * m:
    #   return None
    # if d == 0:
    #   if n == 1 and m == 0:
    #     return Graph(1)
    #   else:
    #     return None
    
    # print('valid')

    G = Graph(n)
    points = [() for _ in range(n)]
    marked = [False for _ in range(n)]
    degrees = [0 for _ in range(n)]
    mat = [[False for _ in range(n)] for _ in range(n)]
    adj = [[] for _ in range(n)]
    check = [[False for _ in range(n)] for _ in range(n)]
    cnt = 0
    # for i in range(m):
    edges = []
    dsu = DSU(n)
    for i in range(n):
      points[i] = (random.randint(0, 1000000000), random.randint(0, 1000000000))
    

    tot_hull_size = 0
    while True:
      p = -1
      idx = []
      for i in range(n):
        if marked[i] == False:
          idx.append(i)
          if p == -1 or points[i][1] < points[idx[p]][1] or (points[i][1] == points[idx[p]][1] and points[i][0] < points[idx[p]][0]):
            p = len(idx) - 1

      if p == -1:
        break
      # idx = [i for i in range(n)]
      (idx[0], idx[p]) = (idx[p], idx[0])
      def cmp(A, B):
        o = self.ccw(points[idx[0]], points[A], points[B])
        if o == 0:
          if math.dist(points[idx[0]], points[A]) < math.dist(points[idx[0]], points[B]):
            return -1
          else:
            return 1
        return -o

      # convex hull
      idx = idx[:1] + sorted(idx[1:], key=cmp_to_key(cmp))
      # for i in idx:
      #   print(points[i], math.atan2(points[i][1] - points[idx[0]][1], points[i][0] - points[idx[0]][0]))

      hull = []
      for i in idx:
        while len(hull) > 1 and self.ccw(points[hull[-2]], points[hull[-1]], points[i]) <= 0:
          hull.pop()
        hull.append(i)
      
      for i in hull:
        marked[i] = True
      #   print(points[i])


      # for (i, e) in enumerate(edges):
      #   if i < len(hull):
      #     break
      #   G.AddEdge(e[0], e[1])
      #   dsu.union(e[0], e[1])

      if len(hull) > 1:
        tot_hull_size += len(hull)
        m = len(hull)
        if m == 2:
          m = 1
        for i in range(m):
          u = hull[i]
          v = hull[(i + 1) % len(hull)]
          edges.append((u, v))
          mat[u][v] = True
          mat[v][u] = True
          degrees[u] += 1
          degrees[v] += 1
          adj[u].append(v)
          adj[v].append(u)
          cnt += 1
          dsu.union(u, v)
          G.AddEdge(u, v)
      else:
        break
    
    # print(edges)
    # print(degrees)

    candidates = set()
    for i in range(n):
      if degrees[i] < d:
        candidates.add(i)

    while True:
      # for _ in range(10000):
        # if cnt * 2 > n and random.randint(0, 10000) < 1000:
        #   for u in range(n):
        #     ok = False
        #     if d > 1 and degrees[u] == 1:
        #       for p in adj[u]:
        #         for v in adj[p]:
        #           if v != u and degrees[v] < d and mat[u][v] == False:
        #             mat[u][v] = mat[v][u] = True
        #             degrees[u] += 1
        #             degrees[v] += 1
        #             if degrees[u] == d:
        #               candidates.remove(u)
        #             if degrees[v] == d:
        #               candidates.remove(v)
        #             adj[u].append(v)
        #             adj[v].append(u)
        #             edges.append((u, v))
        #             cnt += 1
        #             ok = True
        #             break
        #     if ok == True:
        #       break

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
        adj[u].append(v)
        adj[v].append(u)
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

    # if cnt < m:
    #   return None
    

    edges = edges[tot_hull_size:]
    random.shuffle(edges)
    # for i in range(n):
    #   for j in range(n):
    #     mat[i][j] = mat[j][i] = False
    for (i, j) in edges:
      mat[i][j] = mat[j][i] = False

    # print('cnt =', len(G.edges))
    for e in edges:
      # if dsu.union(e[0], e[1]):
        G.AddEdge(e[0], e[1])
        mat[e[0]][e[1]] = True
        mat[e[1]][e[0]] = True
        # if dsu.ncc == 1:
        #   print('connected')
        #   break

    # if dsu.ncc > 1:
    #   return None

    # cnt = m - len(G.edges)
    # random.shuffle(edges)
    # for e in edges:
    #   if cnt <= 0:
    #     break
    #   if mat[e[0]][e[1]] == False:
    #     G.AddEdge(e[0], e[1])
    #     cnt -= 1
    
    # for u in range(n):
    #   if (degrees[u] == 1):
    #     for v in range(n):
    #       if (u != v and mat[u][v] == False and degrees[v] < d):
    #         intersected = False
    #         for e in G.edges:
    #           x = e.fromNode
    #           y = e.toNode
    #           if x == u or x == v or y == u or y == v:
    #             continue
    #           if self.intersect(points[u], points[v], points[x], points[y]):
    #             intersected = True
    #             break
    #         if intersected == True:
    #           continue
    #         G.AddEdge(u, v)
    #         degrees[u] += 1
    #         degrees[v] += 1
    #         mat[u][v] = True
    #         mat[v][u] = True
    #         break
    
    return G


'''
PG = PlanarGraphGenerator()
G = PG.Generate(5,9)
G.Print()  
'''

PG = PlanarGraphGenerator()
G = PG.Generate(30, 50, 1, 3)
if G != None:
  G.Print()
  G.PrintUnweighted(offset = 1)