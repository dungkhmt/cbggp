from CBSGG import DirectedGraph

'''
    TODO by Nguyen Phong 
    
'''

class DirectedGraphGenerator:
    """
    Generate strongly connected directed graphs.
    """

    def __init__(self):
        self.graph = None
        self.listEdges = []
        self.listGraphs = []

    def GenerateRandom(self, n: int, m: int) -> DirectedGraph:
    
        # Generate one random strongly connected directed graph with n nodes and m edges.
        # Idea:
        #   - Step 1: Create a Hamiltonian cycle with a random permutation of nodes.
        #   - Step 2: Randomly add the remaining edges until reaching m edges.
        
        if m < n:
            raise ValueError("Number of edges m must be >= n to guarantee strong connectivity")
        G = DirectedGraph(n)
        import random
        # Step 1: generate a cycle using a random permutation of nodes
        perm = list(range(n))
        random.shuffle(perm) 
        for i in range(n):
            u = perm[i]
            v = perm[(i + 1) % n]  # wrap around to form a cycle
            G.AddEdge(u, v)
        # Step 2: add random edges until having exactly m
        existing = {(e.fromNode, e.toNode) for e in G.edges}
        while G.m < m:
            u = random.randint(0, n - 1)
            v = random.randint(0, n - 1)
            if u != v and (u, v) not in existing:
                G.AddEdge(u, v)
                existing.add((u, v))
        return G

    def BnB(self, idx: int, nbEdges: int):
        """
        Recursive Branch & Bound search to build all directed graphs
        with exactly nbEdges edges.
        """
        # Pruning: not enough edges left to reach nbEdges
        if len(self.listEdges) - idx < nbEdges:
            return
        # End of edge list
        if idx == len(self.listEdges):
            if nbEdges == 0:
                # Check strong connectivity
                if self.is_strongly_connected(self.graph):
                    self.listGraphs.append(self.graph.copy())
            return
        # Case 1: do not take current edge
        self.BnB(idx + 1, nbEdges)
        # Case 2: take current edge (if still needed)
        if nbEdges > 0:
            u, v = self.listEdges[idx]
            self.graph.AddEdge(u, v)
            self.BnB(idx + 1, nbEdges - 1)
            self.graph.PopEdge()
    def GenerateAll(self, n: int, m: int):
        """
        Return all strongly connected directed graphs with n nodes and m edges.
        Warning: grows extremely fast, only usable for very small n (<= 6).
        """
        self.graph = DirectedGraph(n)
        self.listEdges = [(u, v) for u in range(n) for v in range(n) if u != v]
        self.listGraphs = []
        self.BnB(0, m)

        print("Total strongly connected digraphs:", len(self.listGraphs))
        return self.listGraphs

    def is_strongly_connected(self, G: DirectedGraph) -> bool:
        """
        Check strong connectivity using DFS twice (original + reversed graph).
        """
        def dfs(start, adj,c = 0):
            visited = [False] * G.n
            stack = [start]
            visited[start] = True
            while stack:
                u = stack.pop()
                for e in adj[u]:
                    v = G.edges[e].toNode
                    if c :
                        v = G.edges[e].fromNode
                    if not visited[v]:
                        visited[v] = True
                        stack.append(v)
            return visited
        
        # DFS in original graph
        visited1 = dfs(0, G.Adj)
        if not all(visited1):
            return False
        
        # Build reversed adjacency
        rev_adj = [[] for _ in range(G.n)]
        for e in G.edges:
            rev_adj[e.toNode].append(e.id)  # reversed edge
        visited2 = dfs(0, rev_adj,1)

        return all(visited2)
    
if __name__ == "__main__":
    gen = DirectedGraphGenerator()
    # all_graphs = gen.GenerateAll(5, 10)  # all strongly connected digraphs with 3 nodes and 5 edges
    # for g in all_graphs:
    #     g.Print()
    #     print("------")
    G = gen.GenerateRandom(3, 4)
    G.Print()
