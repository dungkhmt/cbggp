import sys
import os
import random
import math

# Add parent directory to sys.path for CBSGG import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph

def generate_connected_bipartite_graph(nb_left_nodes, nb_right_nodes, nb_edges):
  nb_nodes = nb_left_nodes + nb_right_nodes
  G = Graph(nb_nodes)
  
  max_edges = nb_left_nodes * nb_right_nodes
  min_edges = nb_nodes - 1
  
  if nb_left_nodes == 0 or nb_right_nodes == 0:
    if nb_nodes > 1:
      None
          
  if nb_edges < min_edges:
    None
      
  if nb_edges > max_edges:
    None

  left_nodes = list(range(nb_left_nodes))
  right_nodes = list(range(nb_left_nodes, nb_nodes))
  
  random.shuffle(left_nodes)
  random.shuffle(right_nodes)
  
  generated_edges = set()
  
  # BƯỚC 1: Dựng cây khung 
  
  # Rút 1 đỉnh đầu tiên từ mỗi tập để làm gốc của cây
  first_left = left_nodes.pop()
  first_right = right_nodes.pop()
  
  generated_edges.add((first_left, first_right))
  G.AddEdge(first_left, first_right)
  
  visited_left = [first_left]
  visited_right = [first_right]
  
  # Gộp tất cả các đỉnh chưa thăm còn lại và trộn ngẫu nhiên
  remaining_nodes = left_nodes + right_nodes
  random.shuffle(remaining_nodes)
  
  # Duyệt qua từng đỉnh chưa thăm để nối vào cây
  for node in remaining_nodes:
    # Nếu node thuộc tập bên trái 
    if node < nb_left_nodes:
      v = random.choice(visited_right)
      generated_edges.add((node, v))
      G.AddEdge(node, v)
      visited_left.append(node)
    # Nếu node thuộc tập bên phải
    else:
      u = random.choice(visited_left)
      generated_edges.add((u, node))
      G.AddEdge(u, node)
      visited_right.append(node)

  # Thêm các cạnh còn thiếu
  
  edges_to_add = nb_edges - len(generated_edges)
  
  if edges_to_add > 0:
    if max_edges <= 10**6:
      all_possible_edges = [
        (u, v) 
        for u in range(nb_left_nodes) 
        for v in range(nb_left_nodes, nb_nodes)
      ]
      
      available_edges = [e for e in all_possible_edges if e not in generated_edges]
      sampled_edges = random.sample(available_edges, edges_to_add)
      
      for u, v in sampled_edges:
        generated_edges.add((u, v))
        G.AddEdge(u, v)
              
    else:
      while len(generated_edges) < nb_edges:
        u = random.randint(0, nb_left_nodes - 1)
        v = random.randint(nb_left_nodes, nb_nodes - 1)
        edge = (u, v)
        
        if edge not in generated_edges:
          generated_edges.add(edge)
          G.AddEdge(u, v)
                  
  return G
