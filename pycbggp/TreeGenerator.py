from cbggb import Graph, Edge
import random as rd
'''
    TODO by Quyet
    generate undirected tree with constraints
'''
class DSU:
    def __init__(self, n):
        self.parent = [i for i in range(n)]
    def find(self, u):
        if self.parent[u] != u:
            self.parent[u] = self.find(self.parent[u])
        return self.parent[u]
    def union(self, u, v):
        u = self.find(u)
        v = self.find(v)
        if (u == v):
            return 0
        self.parent[v] = u
        return 1

class TreeGenerator:
 # generate randomly an undirected tree with n nodes

    def Generate(self, nbNodes, maxDegree):
        #generate a random tree containsing nbNodes nodes and degree of each node is <= maxDegree 
        freeNode = [i for i in range(nbNodes)]
        treeNode = []
        Degree = [0 for i in range(nbNodes)]
        Tree = Graph(nbNodes)
        # Get root tree
        id = rd.randint(0, nbNodes - 1)
        freeNode[id], freeNode[nbNodes - 1] = freeNode[nbNodes - 1], freeNode[id]
        root = freeNode.pop()
        treeNode.append(root)
        # -------------------------------------------
        while (len(freeNode) > 0):
            n = len(treeNode)
            m = len(freeNode)
            # u from treeNode, v from freeNode
            u = rd.randint(0, n - 1)
            v = rd.randint(0, m - 1)
            valU = treeNode[u]
            valV = freeNode[v]
            Tree.AddEdge(valU, valV)
            # check Degree[u] == MaxDegree ??
            Degree[valU] += 1
            Degree[valV] += 1
            if (Degree[valU] == maxDegree) :
                treeNode[u], treeNode[n - 1] = treeNode[n - 1], treeNode[u]
                treeNode.pop()
            # remove v from freeNode && add v to treeNode
            freeNode[v], freeNode[m - 1] = freeNode[m - 1], freeNode[v]
            freeNode.pop()
            treeNode.append(valV)

        return Tree 
    def GenerateAll(self,n):
        # return all tree containing n nodes 
        edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
        m = len(edges)
        self.graphs = []
        def DSU_copy(d: DSU):
            new = DSU(n)
            new.parent = d.parent[:]
            return new
        def Backtrack(idx : int, NumEdge : int, G : Graph, dsu : DSU):
            if (NumEdge == n - 1):
                self.graphs.append(G)
                return
            if (idx == m) :
                return 
            # Don't select edge(u, v)
            Backtrack(idx + 1, NumEdge, G.copy(), DSU_copy(dsu))
            dsu2 = DSU_copy(dsu)
            u, v = edges[idx]
            # if u and v are not same connected component, add the edge (u, v)
            if (dsu2.union(u, v)) :
                G1 = G.copy()
                G1.AddEdge(u, v)
                Backtrack(idx + 1, NumEdge + 1, G1, dsu2)
        G = Graph(n)
        dsu0 = DSU(n)
        Backtrack(0, 0, G, dsu0)
        return self.graphs 
n = 4
m = 4
CT = TreeGenerator()
RandTree = CT.Generate(n, m)
RandTree.Print()
AllTree = CT.GenerateAll(n)
for Tree in AllTree :
    Tree.Print()
    print("-----------------------")
  
