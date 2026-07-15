import sys
import os
import random
import math

# Add parent directory to sys.path for CBSGG import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph, DirectedGraph

def generate_directed_strongly_connected_graph(nb_nodes, nb_edges):
  G = DirectedGraph(nb_nodes)
  
  if nb_nodes <= 1:
    if nb_edges > 0:
      return None
    return G
  max_edges = nb_nodes * (nb_nodes - 1)
  min_edges = nb_nodes 
  
  if nb_edges < min_edges:
    None
    
  if nb_edges > max_edges:
    None
    
  generated_edges = set()
  
  # Dựng chu trình 
  nodes = list(range(nb_nodes))
  random.shuffle(nodes)
  
  for i in range(nb_nodes):
    u = nodes[i]
    v = nodes[(i + 1) % nb_nodes] # Đỉnh cuối cùng sẽ vòng lại nối với đỉnh đầu tiên
    
    generated_edges.add((u, v))
    G.AddEdge(u, v)
    
  # Thêm các cạnh còn thiếu
  edges_to_add = nb_edges - len(generated_edges)
  
  if edges_to_add > 0:
    if max_edges <= 10**6:
      all_possible_edges = [
        (u, v) 
        for u in range(nb_nodes) 
        for v in range(nb_nodes)
        if u != v
      ]
      
      available_edges = [e for e in all_possible_edges if e not in generated_edges]
      sampled_edges = random.sample(available_edges, edges_to_add)
      
      for u, v in sampled_edges:
        generated_edges.add((u, v))
        G.AddEdge(u, v)
            
    else:
      while len(generated_edges) < nb_edges:
        u = random.randint(0, nb_nodes - 1)
        v = random.randint(0, nb_nodes - 1)
        
        if u != v:
          edge = (u, v)
          if edge not in generated_edges:
            generated_edges.add(edge)
            G.AddEdge(u, v)
            
  return G
