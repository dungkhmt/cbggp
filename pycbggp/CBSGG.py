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

class Point2D:
  def __init__(self,id,x, y):
    self.id = id  
    self.x = x 
    self.y = y    
    
class Edge:
  def __init__(self, fromNode: int, toNode: int, id: int, weight: int = 0, offset: int = 0, used: int = 1):
    self.fromNode = fromNode
    self.toNode = toNode
    self.weight = weight
    self.id = id
    self.offset = offset
    self.used = used

  def copy(self):
    return Edge(self.fromNode, self.toNode, self.id, self.weight, self.offset, self.used)

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

  def __init__(self, n: int, offset: int = 0):
    self.n = n
    self.offset = offset
    self.edges: list[Edge] = []
    self.Adj: list[list[int]] = [[] for _ in range(n)]
    self.m = 0

  def AddEdge(self, fromNode: int, toNode: int, weight: int = 0):
    e = Edge(fromNode, toNode, self.m, weight, self.offset)
    self.Adj[fromNode].append(self.m)
    self.Adj[toNode].append(self.m)
    self.edges.append(e)
    self.m += 1
  
  def gen_random_weight_edges(self,min_w,max_w):
    for e in self.edges:
      e.weight = rd.randint(min_w,max_w)
      print('edge (',e.fromNode,',',e.toNode,') weight = ',e.weight)
      
  def AddEdge1(self, e: Edge):
    self.Adj[e.fromNode].append(self.m)
    self.Adj[e.toNode].append(self.m)
    self.edges.append(e)
    e.offset = self.offset
    self.m += 1

  def addFromAdj(self, adj: list[list[int]]):
    for u in range(self.n):
      for v in adj[u]:
        if u < v:
          self.AddEdge(u, v)

  def PopEdge(self):
    if self.m == 0:
      return None
    e = self.edges.pop()
    self.Adj[e.fromNode].pop()
    self.Adj[e.toNode].pop()
    self.m -= 1
    return e
  
  def DeleteEdge(self, idEdge: int):
    e: Edge = self.edges[idEdge]
    e.used = 0
    self.Adj[e.fromNode].remove(e.id)
    self.Adj[e.toNode].remove(e.id)
    return e

  def CopyGraphMapNewNodes(self, nodeIds):
    # by DungNT
    # map id of the current graph to new nodeIds

    max_id = -1
    min_id = int(1e9)
    for id in nodeIds:
      if id > max_id:
        max_id = id
      if id < min_id:
        min_id = id

    n = max_id - min_id + 1

    G = Graph(n, min_id)
    for e in self.edges:
      G.AddEdge(nodeIds[e.fromNode] - min_id,
                nodeIds[e.toNode] - min_id, e.weight)

    return G

  def CopyGraphMapNewNodes(self, fromId, toId):
    # by DungNT
    # map id of the current graph to new node ids in the range [fromId..toId]

    n = toId - fromId + 1
    assert n == self.n, 'Error: range_id of new graph must be equal to the original graph!'
    G = Graph(n, fromId)
    for e in self.edges:
      G.AddEdge(e.fromNode, e.toNode, e.weight)

    return G

  def union_two_graphs(self,g1, g2):
    # by PQD, to be checked carefully the ids consistence  
    n = g1.n + g2.n 
    G = Graph(n)
    for e in g1.edges:
      G.AddEdge(e.fromNode, e.toNode, e.weight)  
    for e in g2.edges:
      G.AddEdge(e.fromNode, e.toNode, e.weight)  

    return G
    
  def LinkTwoGraphMergeNode(self, G1, node1, G2, node2):
    # resulting graph G has nodes = (V(G1)\{nodes}) unify (V(G2)\{node2}) unify node (node1 fusing node2)
    # merge G2 into G1 and return G1, offset is the same as G1

    # TODO by DungNT 
    
    n = G1.n + G2.n - 1
    G = Graph(n, G1.offset)
    for e in G1.edges:
      G.AddEdge(e.fromNode, e.toNode, e.weight)
    delta = G1.n + G1.offset - G2.offset
    for e in G2.edges:
      o1 = 1 if e.fromNode > node2 else 0
      o2 = 1 if e.toNode > node2 else 0
      if e.fromNode != node2 and e.toNode != node2:
        G.AddEdge(e.fromNode + delta - o1, e.toNode + delta - o2, e.weight)
      elif e.fromNode == node2:
        G.AddEdge(node1, e.toNode + delta - o2, e.weight)
      else: 
        G.AddEdge(e.fromNode + delta - o1, node1, e.weight)

    return G
    
  def LinkTwoGraphs(self, G1, node1, G2, node2, is_offset = True):
    # by DungNT
    # return a new created Graph by adding an edge connecting node1 of G1 and node2 of G2

    if G1.offset > G2.offset:
      tmp = G1
      G1 = G2
      G2 = tmp

    delta = G2.offset - G1.offset
    n = max(G1.n, G2.n + delta)

    G = Graph(n, G1.offset)
    for e in G1.edges:
      G.AddEdge(e.fromNode, e.toNode, e.weight)
    
    for e in G2.edges:
      G.AddEdge(e.fromNode + delta, e.toNode + delta, e.weight)

    if is_offset:
      G.AddEdge(node1, node2 + delta)
    else:
      G.AddEdge(node1 - G1.offset, node2 - G1.offset)

    return G

  def Print(self):
    print('Nodes = ', end=' ')
    for v in range(self.n):
      print(v + self.offset, end=' ')
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

  def PrintWeighted(self, offset=0):
    print(self.n, self.m)
    for e in self.edges:
      print(e.fromNode + offset + self.offset,
            e.toNode + offset + self.offset, e.weight)

  def PrintUnweighted(self, offset=0):
    print(self.n, self.m)
    for e in self.edges:
      print(e.fromNode + offset + self.offset, e.toNode + offset + self.offset)

  def PrintNodesEdgesFormat(self):
    V = '('
    for i in range(self.n):
      V = V + str(i + self.offset)
      if i < self.n-1:
        V = V + ','
    V = V + ')'
    print(V)
    E = '{'
    for i in range(len(self.edges)):
      e = self.edges[i]
      E = E + '(' + str(e.fromNode + self.offset) + \
          ',' + str(e.toNode + self.offset) + ')'
      if i < len(self.edges)-1:
        E = E + ','
    E = E + '}'
    print(E)

  def SaveToFile(self, filename):
    with open(filename, 'w') as f:
      # print('Nodes = ', end=' ')
      line = ''
      for v in range(self.n):
        # print(v, end=' ')
        line = line + str(v + self.offset) + ' '
      f.write(line + '\n')
      # for v in self.nodes:
      #  print(v.toStr(),end = ' ')
      # print('')
      # for n in self.nodes:
      #  print('Adj[' + n.toStr() + ']: ',end = ' ')

      # print('Edges = ', end=' ')

      for e in self.edges:
        line = str(e.fromNode + self.offset) + \
            ' ' + str(e.toNode + self.offset)
        # print(e.toStr(), end=' ')
        f.write(line + '\n')
    f.close()
    return

  def SaveToFileWeighted(self, filename):
    with open(filename, 'w') as f:
      # print('Nodes = ', end=' ')
      line = ''
      for v in range(self.n):
        # print(v, end=' ')
        line = line + str(v + self.offset) + ' '
      f.write(line + '\n')
      # for v in self.nodes:
      #  print(v.toStr(),end = ' ')
      # print('')
      # for n in self.nodes:
      #  print('Adj[' + n.toStr() + ']: ',end = ' ')

      # print('Edges = ', end=' ')

      for e in self.edges:
        line = str(e.fromNode + self.offset) + ' ' + \
            str(e.toNode + self.offset) + ' ' + str(e.weight)
        # print(e.toStr(), end=' ')
        f.write(line + '\n')
    f.close()
    return

  def copy(self):
    G = Graph(self.n, self.offset)
    for e in self.edges:
      if e.used == 1:
        G.AddEdge(e.fromNode, e.toNode, e.weight)
    return G


class DirectedGraph(Graph):
  def __init__(self, n: int, offset: int = 0):
    super().__init__(n, offset)

  def AddEdge(self, fromNode: int, toNode: int, weight: int = 0):
    e = Edge(fromNode, toNode, self.m, weight, self.offset)
    self.Adj[fromNode].append(self.m)
    self.edges.append(e)
    self.m += 1

  def AddEdge1(self, e: Edge):
    self.Adj[e.fromNode].append(self.m)
    self.edges.append(e)
    e.offset = self.offset
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
      print(v + self.offset, end=' ')
    print('')

    print('Edges = ', end=' ')
    for e in self.edges:
      print(e.toStr(), end=' ')
    print('')

    for n in range(self.n):
      print('Adj[' + str(n + self.offset) + ']: ', end=' ')

      for e in self.Adj[n]:
        print(self.edges[e].toStr(), end=' ')
      print('')

  def PrintWeighted(self, offset=0):
    print(self.n, self.m)
    for e in self.edges:
      print(e.fromNode + offset + self.offset,
            e.toNode + offset + self.offset, e.weight)

  def PrintUnweighted(self, offset=0):
    print(self.n, self.m)
    for e in self.edges:
      print(e.fromNode + offset + self.offset, e.toNode + offset + self.offset)

  def addFromAdj(self, adj: list[list[int]]):
    for u in range(self.n):
      for v in adj[u]:
        self.AddEdge(u, v)

  def CopyGraphMapNewNodes(self, nodeIds):
    # by DungNT
    # map id of the current graph to new nodeIds

    max_id = -1
    min_id = int(1e9)
    for id in nodeIds:
      if id > max_id:
        max_id = id
      if id < min_id:
        min_id = id

    n = max_id - min_id + 1

    G = DirectedGraph(n, min_id)
    for e in self.edges:
      G.AddEdge(nodeIds[e.fromNode] - min_id,
                nodeIds[e.toNode] - min_id, e.weight)

    return G

  def CopyGraphMapNewNodes(self, fromId, toId):
    # by DungNT
    # map id of the current graph to new node ids in the range [fromId..toId]

    n = toId - fromId + 1
    assert n == self.n, 'Error: range_id of new graph must be equal to the original graph!'
    G = DirectedGraph(n, fromId)
    for e in self.edges:
      G.AddEdge(e.fromNode, e.toNode, e.weight)

    return G


  def LinkTwoGraphMergeNode(self, G1, node1, G2, node2):
    # resulting graph G has nodes = (V(G1)\{nodes}) unify (V(G2)\{node2}) unify node (node1 fusing node2)
    # TODO by DungNT 
    
    n = G1.n + G2.n - 1
    G = DirectedGraph(n, G1.offset)
    for e in G1.edges:
      G.AddEdge(e.fromNode, e.toNode, e.weight)
    delta = G1.n + G1.offset - G2.offset
    for e in G2.edges:
      o1 = 1 if e.fromNode > node2 else 0
      o2 = 1 if e.toNode > node2 else 0
      if e.fromNode != node2 and e.toNode != node2:
        G.AddEdge(e.fromNode + delta - o1, e.toNode + delta - o2, e.weight)
      elif e.fromNode == node2:
        G.AddEdge(node1, e.toNode + delta - o2, e.weight)
      else: 
        G.AddEdge(e.fromNode + delta - o1, node1, e.weight)

    return G
    
  def LinkTwoGraphs(self, G1, node1, G2, node2, is_offset = True):
    # by DungNT
    # return a new created Graph by adding an edge connecting node1 of G1 and node2 of G2

    if G1.offset > G2.offset:
      tmp = G1
      G1 = G2
      G2 = tmp

    delta = G2.offset - G1.offset
    n = max(G1.n, G2.n + delta)

    G = DirectedGraph(n, G1.offset)
    for e in G1.edges:
      G.AddEdge(e.fromNode, e.toNode, e.weight)
    
    for e in G2.edges:
      G.AddEdge(e.fromNode + delta, e.toNode + delta, e.weight)

    if is_offset:
      G.AddEdge(node1, node2 + delta)
    else:
      G.AddEdge(node1 - G1.offset, node2 - G1.offset)

    return G

  def SaveToFile(self, filename):
    with open(filename, 'w') as f:
      # print('Nodes = ', end=' ')
      line = ''
      for v in range(self.n):
        # print(v, end=' ')
        line = line + str(v + self.offset) + ' '
      f.write(line + '\n')
      # for v in self.nodes:
      #  print(v.toStr(),end = ' ')
      # print('')
      # for n in self.nodes:
      #  print('Adj[' + n.toStr() + ']: ',end = ' ')

      # print('Edges = ', end=' ')

      for e in self.edges:
        line = str(e.fromNode + self.offset) + \
            ' ' + str(e.toNode + self.offset)
        # print(e.toStr(), end=' ')
        f.write(line + '\n')
    f.close()
    return

  def SaveToFileWeighted(self, filename):
    with open(filename, 'w') as f:
      # print('Nodes = ', end=' ')
      line = ''
      for v in range(self.n):
        # print(v, end=' ')
        line = line + str(v + self.offset) + ' '
      f.write(line + '\n')
      # for v in self.nodes:
      #  print(v.toStr(),end = ' ')
      # print('')
      # for n in self.nodes:
      #  print('Adj[' + n.toStr() + ']: ',end = ' ')

      # print('Edges = ', end=' ')

      for e in self.edges:
        line = str(e.fromNode + self.offset) + ' ' + \
            str(e.toNode + self.offset) + ' ' + str(e.weight)
        # print(e.toStr(), end=' ')
        f.write(line + '\n')
    f.close()
    return

  def copy(self):
    G = DirectedGraph(self.n)
    for e in self.edges:
      G.AddEdge1(e)
    return G


class RootedTree:
    # todo by Chi Vi
    """A small weighted rooted tree.
    Edges are undirected and represented as ``(u, v, weight)``.  The root is
    used only to determine parent/child relationships; all distances use the
    supplied edge weights.
    """

    def __init__(self, root: int):
        self.root = root
        self.edges = []
        self._adj = {root: {}}

    @staticmethod
    def _edge_parts(edge, weight=None):
        if weight is not None:
            if len(edge) != 2:
                return None
            u, v = edge
            return u, v, weight
        if len(edge) == 2:
            return edge[0], edge[1], 1
        if len(edge) == 3:
            return edge[0], edge[1], edge[2]
        return None

    def _check_weight(self, weight):
        if not isinstance(weight, (int, float)):
            raise TypeError("edge weight must be a number")
        if weight < 0:
            return False
        return True

    def _has_edge(self, u, v):
        return u in self._adj and v in self._adj[u]

    def addEdge(self, edge, weight=None):
        """Add a leaf edge and return ``True``.
        Exactly one endpoint must already belong to the tree.  This keeps the
        object a rooted tree without needing a general graph connectivity
        check.
        """
        parts = self._edge_parts(edge, weight)
        if parts is None:
            return
        u, v, edge_weight = parts
        if not self._check_weight(edge_weight):
            return
        if u == v:
            return
        if self._has_edge(u, v):
            return
        u_exists = u in self._adj
        v_exists = v in self._adj
        if u_exists == v_exists:
            return

        self._adj.setdefault(u, {})[v] = edge_weight
        self._adj.setdefault(v, {})[u] = edge_weight
        self.edges.append((u, v, edge_weight))
        return True

    def add_child(self, parent_node: int, new_node: int, weight=1):
        return self.addEdge((parent_node, new_node), weight)

    def RemoveEdge(self, edge, weight=None):
        """Remove a leaf edge and return ``True``.
        Removing an internal edge would split the tree, so it is rejected.
        The weight is ignored when matching an existing edge.
        """
        parts = self._edge_parts(edge, weight)
        if parts is None:
            return
        u, v, _ = parts
        if not self._has_edge(u, v):
            return
        if len(self._adj[u]) > 1 and len(self._adj[v]) > 1:
            return

        del self._adj[u][v]
        del self._adj[v][u]
        if not self._adj[u] and u != self.root:
            del self._adj[u]
        if not self._adj[v] and v != self.root:
            del self._adj[v]
        self.edges.remove(next(e for e in self.edges
                               if {e[0], e[1]} == {u, v}))
        return True

    def ReplaceEdge(self, edge_in, edge_out):
        #Replace one leaf edge with another edge, preserving the tree.
        in_parts = self._edge_parts(edge_in)
        out_parts = self._edge_parts(edge_out)
        if in_parts is None or out_parts is None:
            return
        u, v, new_weight = in_parts
        old_u, old_v, _ = out_parts
        if not self._has_edge(old_u, old_v):
            return
        if len(self._adj[old_u]) > 1 and len(self._adj[old_v]) > 1:
            return
        if self.addEdge((u, v, new_weight)) is None:
            return
        self.RemoveEdge((old_u, old_v))
        return True

    def _rooted_distances(self):
        distances = {self.root: 0}
        parent = {self.root: None}
        stack = [self.root]
        while stack:
            node = stack.pop()
            for child, weight in self._adj[node].items():
                if child == parent[node]:
                    continue
                parent[child] = node
                distances[child] = distances[node] + weight
                stack.append(child)
        return distances, parent

    def _rooted_depths(self):
        depths = {self.root: 0}
        parent = {self.root: None}
        stack = [self.root]
        while stack:
            node = stack.pop()
            for child in self._adj[node]:
                if child == parent[node]:
                    continue
                parent[child] = node
                depths[child] = depths[node] + 1
                stack.append(child)
        return depths

    def height(self):
        #Return the maximum number of edges from root to any node.
        return max(self._rooted_depths().values())

    def weighted_height(self):
        #Return the maximum weighted distance from root to any node.
        return max(self._rooted_distances()[0].values())

    def diameter(self):
        #Return the maximum weighted distance between any two nodes.
        if len(self._adj) <= 1:
            return 0

        def farthest(start):
            distances = {start: 0}
            stack = [start]
            parent = {start: None}
            while stack:
                node = stack.pop()
                for neighbour, weight in self._adj[node].items():
                    if neighbour == parent[node]:
                        continue
                    parent[neighbour] = node
                    distances[neighbour] = distances[node] + weight
                    stack.append(neighbour)
            return max(distances.items(), key=lambda item: item[1])

        endpoint, _ = farthest(self.root)
        return farthest(endpoint)[1]

    def number_of_children(self, node):
        _, parent = self._rooted_distances()
        if node not in self._adj:
            return
        return sum(child != parent[node] for child in self._adj[node])

    def subtree_size(self, node):
        distances, parent = self._rooted_distances()
        if node not in distances:
            return
        return sum(
            1 for candidate in distances
            if candidate == node or self._is_descendant(candidate, node, parent)
        )

    @staticmethod
    def _is_descendant(node, ancestor, parent):
        while parent[node] is not None:
            node = parent[node]
            if node == ancestor:
                return True
        return False

    # Readable aliases for callers using snake_case.
    add_edge = addEdge
    remove_edge = RemoveEdge
    replace_edge = ReplaceEdge
    get_height = height
    get_weighted_height = weighted_height
    get_diameter = diameter
    count_children = number_of_children
    count_subtree_nodes = subtree_size

    def print(self):
        for edge in self.edges:
            print(edge)
        
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


if __name__ == "__main__" and False:
  connected = 'Y'
  nbNodes = 0
  nbEdges = 0
  nbBridges = 0
  maxDegree = 10
  gType = 'planar'
  directed = 'N'
  filename = 'g.txt'
  for i in range(len(sys.argv)):
    # print('argv[',i,'] = ',sys.argv[i])
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
    G = PGG.Generate(nbNodes, nbEdges, nbBridges, maxDegree)
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


# Chi Vi
def test_rooted_tree():
  tree = RootedTree(0)

  # Add edges and verify the basic tree structure.
  assert tree.addEdge((0, 1, 2)) is True
  assert tree.addEdge((0, 2, 3)) is True
  assert tree.addEdge((1, 3, 4)) is True
  assert tree.edges == [(0, 1, 2), (0, 2, 3), (1, 3, 4)]

  # Rooted metrics.
  assert tree.height() == 2
  assert tree.weighted_height() == 6
  assert tree.diameter() == 9
  assert tree.number_of_children(0) == 2
  assert tree.number_of_children(1) == 1
  assert tree.number_of_children(3) == 0
  assert tree.subtree_size(0) == 4
  assert tree.subtree_size(1) == 2
  assert tree.subtree_size(3) == 1

  # Snake_case aliases and add_child.
  assert tree.add_child(2, 4, 1) is True
  assert tree.get_height() == 2
  assert tree.get_weighted_height() == 6
  assert tree.count_children(2) == 1
  assert tree.count_subtree_nodes(2) == 2

  # Replace a leaf edge and verify that its new weight is used.
  assert tree.ReplaceEdge((1, 5, 5), (1, 3, 4)) is True
  assert (1, 3, 4) not in tree.edges
  assert (1, 5, 5) in tree.edges
  assert tree.get_height() == 2
  assert tree.get_weighted_height() == 7
  assert tree.get_diameter() == 11

  # Remove a leaf edge.
  assert tree.RemoveEdge((1, 5, 5)) is True
  assert tree.remove_edge((2, 4, 1)) is True
  assert tree.edges == [(0, 1, 2), (0, 2, 3)]

  print("RootedTree: all tests passed")


if __name__ == "__main__":
  test_rooted_tree()