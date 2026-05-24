from PlanarGraphGenerator import PlanarGraphGenerator
from CBSGG import Graph


PG = PlanarGraphGenerator()
G = PG.Generate(30, 50, 1, 3)
if G != None:
  G.Print()
  G.PrintUnweighted(offset = 1)