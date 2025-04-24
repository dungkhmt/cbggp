import networkx as nx
import random
import matplotlib.pyplot as plt

def generate_connected_graph(num_nodes, num_edges):
    """
    Generates a connected graph with the specified number of nodes and edges.

    Args:
        num_nodes: The number of nodes in the graph.
        num_edges: The number of edges in the graph.

    Returns:
        A networkx Graph object if successful, None otherwise.
    """
    if num_edges < num_nodes - 1 or num_edges > num_nodes * (num_nodes - 1) / 2:
        print("Invalid number of edges.  Must be at least num_nodes - 1 and at most n*(n-1)/2 for simple graph")
        return None

    # Start with an empty graph
    graph = nx.Graph()
    # Add the nodes
    graph.add_nodes_from(range(num_nodes))

    # 1. Ensure connectivity: Create a spanning tree first.
    #    This guarantees the graph is connected.
    for i in range(num_nodes - 1):
        graph.add_edge(i, i + 1)  # Creates a simple path (0-1-2-...-9)

    # 2. Add the remaining edges randomly.
    remaining_edges = num_edges - (num_nodes - 1)
    possible_edges = [(u, v) for u in range(num_nodes) for v in range(u + 1, num_nodes)]
    # Remove the edges already added in the spanning tree part
    possible_edges = [e for e in possible_edges if not graph.has_edge(e[0],e[1])]

    if remaining_edges > len(possible_edges):
        print(f"Cannot create {num_edges} edges.  Maximum possible is {len(possible_edges) + (num_nodes-1)}")
        return None
    
    edges_to_add = random.sample(possible_edges, remaining_edges) # Faster and avoids duplicates
    graph.add_edges_from(edges_to_add)
    
    # Double check connectivity (optional, but good for testing)
    if not nx.is_connected(graph):
        print("Error: Graph is not connected (this should not happen).  Retrying...")
        return None # Indicate failure, so we can try again.

    return graph

if __name__ == "__main__":
    num_nodes = 10
    num_edges = 15
    graph = generate_connected_graph(num_nodes, num_edges)

    if graph is not None:
        print(f"Generated a connected graph with {num_nodes} nodes and {num_edges} edges.")
        # Visualize the graph
        plt.figure(figsize=(8, 6))  # Adjust figure size for better visualization
        nx.draw(graph, with_labels=True, node_color='skyblue', edge_color='gray', node_size=500)
        plt.title(f"Connected Graph ({num_nodes} Nodes, {num_edges} Edges)")
        plt.show()
    else:
        print("Failed to generate a connected graph.")
