import sys
import os
import random
import math

# Add parent directory to sys.path for CBSGG import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from CBSGG import Graph

def generate_undirected_complete_graph(nb_nodes) :
  G = Graph(nb_nodes)

  for i in range(nb_nodes):
    for j in range(i + 1, nb_nodes):
      G.AddEdge(i, j)
  
  return G
