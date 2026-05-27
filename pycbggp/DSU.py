class DSU:
  def __init__(self, n, is_rolled_back = False):
    self.n = n
    self.parent = list(range(n))
    self.size = [1] * n
    self.is_rolled_back = is_rolled_back 
    self.history = []
  
  def find(self, u):
    if self.parent[u] != u:
      if self.is_rolled_back:
        return self.find(self.parent[u])
      self.parent[u] = self.find(self.parent[u])
    return self.parent[u]
  
  def union(self, u, v):
    u = self.find(u)
    v = self.find(v)
    if u == v:
      return False
    if self.size[u] < self.size[v]:
      u, v = v, u
    self.parent[v] = u
    self.size[u] += self.size[v]
    if self.is_rolled_back:
      self.history.append(v)
    return True

  def roll(self):
    assert self.is_rolled_back
    if not self.history:
      return False
    v = self.history.pop()
    self.size[self.parent[v]] -= self.size[v]
    self.parent[v] = v
    return True
  