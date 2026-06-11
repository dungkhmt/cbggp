import sys
import random
from collections import namedtuple
from typing import List
from CBSGG import Graph

Constants = namedtuple('Constants', ['MAX_VALUE'])
constants = Constants(MAX_VALUE = int(1e9))

def power(a: int, b: int):
    if b == 0:
        return 1
    s = power(a, b // 2)
    s = s * s
    if b % 2 == 1:
        s = s * a
    if s > constants.MAX_VALUE:
        return constants.MAX_VALUE
    return s

class UndirectedTreeDiameterBetweenPandQAndDegreeAtMostDGenerator:
    '''
        TODO by Nguyen Ngoc Tuan Anh
    '''

    def check_constraint(self, n: int, d: int, D: int):
        if n <= 0 or d < 0 or D < 0:
            return False
        if D == 0:
            return n == 1 and d == 0
        if n - 1 < d:
            return False
        if D == 1:
            if n > 2 or d > 1:
                return False
            return True
        if D == 2:
            if d < n - 1:
                return False
            return True
        h = d // 2
        if d % 2 == 0 and n > 1 + D * (power(D - 1, h) - 1) // (D - 2):
            return False
        if d % 2 == 1 and n > 2 * (power(D - 1, h + 1) - 1) // (D - 2):
            return False
        return True

    def check(self, n: int, p: int, q: int, d: int):
        if n <= 0 or p > q or d < 0:
            return False
        lo = max(0, p)
        hi = min(q, n - 1)
        if lo > hi:
            return False
        for diameter in range(lo, hi + 1):
            if self.check_constraint(n, diameter, d):
                return True
        return False

    def initialize(self, n: int, p: int, q: int, d: int):
        G = Graph(n, 0)
        candidates_diameter = []
        lo = max(0, p)
        hi = min(q, n - 1)
        for diameter in range(lo, hi + 1):
            if self.check_constraint(n, diameter, d):
                candidates_diameter.append(diameter)
        if len(candidates_diameter) == 0:
            return None
        diameter = random.choice(candidates_diameter)
        for i in range (diameter):
            G.AddEdge(i, i + 1)
        max_distance = [0 for _ in range(n)]
        for i in range (diameter + 1):
            max_distance[i] = max(i, diameter - i)
        candidates_nodes = []
        for i in range (diameter + 1):
            if max_distance[i] < diameter and len(G.Adj[i]) < d:
                candidates_nodes.append(i)
        for i in range (diameter + 1, n):
            if (len(candidates_nodes) == 0):
                return None
            index = random.randint(0, len(candidates_nodes) - 1)
            node = candidates_nodes[index]
            G.AddEdge(node, i)
            if len(G.Adj[node]) == d:
                candidates_nodes[index] = candidates_nodes[-1]
                candidates_nodes.pop()
            max_distance[i] = max_distance[node] + 1
            if max_distance[i] < diameter and len(G.Adj[i]) < d:
                candidates_nodes.append(i)
        return G

    def dfs(self, G: Graph, s: int, d: int):
        cmp: List[int] = []
        vis = [0 for _ in range (G.n)]
        stack: List[int] = [s]
        while len(stack) > 0:
            u = stack.pop()
            if len(G.Adj[u]) < d:
                cmp.append(u)
            vis[u] = 1
            for idEdge in G.Adj[u]:
                if G.edges[idEdge].used == 0:
                    continue
                v = G.edges[idEdge].fromNode + G.edges[idEdge].toNode - u
                if vis[v] > 0:
                    continue;
                stack.append(v)
        return cmp
    
    def get_diameter(self, G: Graph):
        def dfs(s: int):
            dist = [-1 for _ in range (G.n)]
            dist[s] = 0
            stack: List[int] = [s]
            while len(stack) > 0:
                u = stack.pop()
                for idEdge in G.Adj[u]:
                    if G.edges[idEdge].used == 0:
                        continue
                    v = G.edges[idEdge].fromNode + G.edges[idEdge].toNode - u
                    if dist[v] >= 0:
                        continue;
                    dist[v] = dist[u] + 1
                    stack.append(v)
            return dist
        
        dist = dfs(0)
        s = dist.index(max(dist))
        
        dist = dfs(s)
        t = dist.index(max(dist))
        return dist[t]

    def generate(self, n: int, p: int, q: int, d: int):
        G_cur = self.initialize(n, p, q, d)
        if G_cur is None:
            return None
        # print('Initial graph: ', end = '')
        # G_cur.Print()
        for _ in range (100):
            G_new: Graph = G_cur.copy()
            if G_new is None:
                break
            idEdge = random.randint(0, G_new.m - 1);
            if G_new.edges[idEdge].used == 0:
                print('Error: edge with id ' + str(idEdge) + ' has used = 0')
                break

            e = G_new.DeleteEdge(idEdge)
            cmp1 = self.dfs(G_new, e.fromNode, d)
            cmp2 = self.dfs(G_new, e.toNode, d)
            if len(cmp1) == 0 or len(cmp2) == 0:
                continue
            node1 = random.choice(cmp1)
            node2 = random.choice(cmp2)
            G_new.AddEdge(node1, node2)
            diameter_new = self.get_diameter(G_new)
            if p <= diameter_new <= q:
                # print('-----------Iteration ' + str(_ + 1) + '-----------')
                # print('Current graph: ', end = '')
                # G_cur.Print()
                # print('Delete edge: ', e.toStr())
                # print('New graph: ', end = '')
                # G_new.Print()
                G_cur = G_new.copy()
        
        G_ans = Graph(n, 1)
        nodes = [i for i in range (n)]
        random.shuffle(nodes)
        for e in G_cur.edges:
            if e.used == 1:
                G_ans.AddEdge(nodes[e.fromNode], nodes[e.toNode])

        return G_ans

# def _token_generator():
#     for line in sys.stdin:
#         for tok in line.strip().split():
#             yield tok

# tokens = _token_generator()
# try:
#     n = int(next(tokens))
#     p = int(next(tokens))
#     q = int(next(tokens))
#     d = int(next(tokens))
# except (StopIteration, ValueError):
#     print("-1")
# else:
#     gen = UndirectedTreeDiameterBetweenPandQAndDegreeAtMostDGenerator()
#     G = gen.generate(n, p, q, d)
#     if G is None:
#         print("-1")
#     else:
#         G.Print()