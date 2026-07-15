import sys
import os
import random
import math

# Add parent directory to sys.path for CBSGG import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph, DirectedGraph

def generate_directed_graph(nb_nodes, nb_edges):
  G = DirectedGraph(nb_nodes)
  
  max_edges = nb_nodes * (nb_nodes - 1)
  
  if nb_edges > max_edges:
    return None
    
  # Trường hợp 1: Đối với đồ thị dày
  if max_edges <= 10**6:
    all_possible_edges = [
      (u, v) 
      for u in range(nb_nodes) 
      for v in range(nb_nodes)
      if u != v # Đảm bảo không tạo khuyên (self-loop)
    ]
    
    sampled_edges = random.sample(all_possible_edges, nb_edges)
    
    for u, v in sampled_edges:
      G.AddEdge(u, v)
          
  # Trường hợp 2: Đồ thị lớn và thưa 
  else:
    generated_edges = set()
    
    while len(generated_edges) < nb_edges:
      u = random.randint(0, nb_nodes - 1)
      v = random.randint(0, nb_nodes - 1)
      
      if u != v: # Bỏ qua nếu tạo ra khuyên
        edge = (u, v)
        
        if edge not in generated_edges:
          generated_edges.add(edge)
          G.AddEdge(u, v)

  return G
