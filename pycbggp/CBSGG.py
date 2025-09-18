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
 def __init__(self,fromNode, toNode):
  self.fromNode = fromNode 
  self.toNode = toNode 
  self.weight = 0
 
 def toStr(self):
  return '(' + self.fromNode.toStr() + ',' + self.toNode.toStr() + ')' 
  
class Graph:
 '''
    TO BE DONE by Nguyen Tan Dung
    maintain: a list of edges E 
               A[v] is the list of indices of adjacent (outgoing) edges 
 '''
 def __init__(self, nodes):
  self.nodes = nodes 
  self.Adj = {} #adjacent edges
  for v in nodes:
   self.Adj[v] = []
   
 def AddEdge(self, e):
  self.Adj[e.fromNode].append(e)
  self.Adj[e.toNode].append(e)
  
 def Print(self):
  print('Nodes = ',end = ' ')
  for v in self.nodes:
   print(v.toStr(),end = ' ')  
  print('')
  for n in self.nodes:
   print('Adj[' + n.toStr() + ']: ',end = ' ')
   
   for e in self.Adj[n]:
    print(e.toStr(), end = ' ')
   print('') 
   
   
class SubGraphGenerator:
 def __init__(self,G):
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
 nod = Node(i)
 nodes.append(nod)

G = Graph(nodes)
G.AddEdge(Edge(nodes[0],nodes[1]))
G.AddEdge(Edge(nodes[0],nodes[2]))
G.AddEdge(Edge(nodes[1],nodes[3]))
G.AddEdge(Edge(nodes[2],nodes[4]))

G.Print() 