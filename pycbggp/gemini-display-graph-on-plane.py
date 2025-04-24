import networkx as nx
import matplotlib.pyplot as plt

# Define the nodes and edges
nodes = [1, 2, 3, 4, 5]
edges = [(1, 2), (1, 3), (1, 5), (2, 4), (3, 5)]

# Create the graph
graph = nx.Graph()
graph.add_nodes_from(nodes)
graph.add_edges_from(edges)

# Choose a layout for the graph.  Spring layout generally looks nice.
pos = nx.spring_layout(graph)  # You can experiment with other layouts like:
# pos = nx.circular_layout(graph)
# pos = nx.random_layout(graph)
# pos = nx.shell_layout(graph)
# pos = nx.spectral_layout(graph)

# Draw the graph
plt.figure(figsize=(8, 6))  # Adjust figure size for better visualization
nx.draw(graph, pos, with_labels=True, node_color='skyblue', edge_color='gray', node_size=500)
plt.title("Graph with Nodes {1, 2, 3, 4, 5} and Edges {(1, 2), (1, 3), (1, 5), (2, 4), (3, 5)}")
plt.show()
