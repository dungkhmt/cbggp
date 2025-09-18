import sys
'''
    Generate graphs under constraints
    -connected undirected graph with n nodes and m edges and k bridges
    -strongly connected graph with n nodes and m edges
    -tree with bounded diameter, degree
    -directed acyclic graph (DAG) with n nodes and m edges 
    -connected graph with n nodes and m edges and k bi-connected components
    -planar graph with constraint
    ...
'''

def generateConnectedSubgraph():
 return None 
 
if __name__ == "__main__":
     
   
 connected = 'Y'
 nbNodes = 0
 nbEdges = 0
 gType = 'planar'
 directed = 'N' 
 for i in range(len(sys.argv)):
  #print('argv[',i,'] = ',sys.argv[i]) 
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
 