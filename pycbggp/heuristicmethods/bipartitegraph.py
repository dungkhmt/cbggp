import sys
import os
import random
import math

# Add parent directory to sys.path for CBSGG import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph

def generate_bipartite_graph(nb_left_nodes, nb_right_nodes, nb_edges):
  nb_nodes = nb_left_nodes + nb_right_nodes
  G = Graph(nb_nodes)
  
  max_edges = nb_left_nodes * nb_right_nodes
  
  if nb_edges > max_edges:
    return None
  if max_edges <= 10**6:
    # Tạo danh sách tất cả các cạnh có thể có
    all_possible_edges = [
      (u, v) 
      for u in range(nb_left_nodes) 
      for v in range(nb_left_nodes, nb_nodes)
    ]
    
    # Dùng random.sample để nhặt ra đúng nb_edges cạnh mà không bị trùng lặp
    sampled_edges = random.sample(all_possible_edges, nb_edges)
    
    for u, v in sampled_edges:
      G.AddEdge(u, v)
  else:
    generated_edges = set()
    
    while len(generated_edges) < nb_edges:
      u = random.randint(0, nb_left_nodes - 1)
      v = random.randint(nb_left_nodes, nb_nodes - 1)
      
      edge = (u, v)
      
      if edge not in generated_edges:
          generated_edges.add(edge)
          G.AddEdge(u, v)
              
  return G
