import networkx as nx
import random
import matplotlib.pyplot as plt

def generate_graph_with_bridges_and_component_size(num_nodes, num_edges, num_bridges, min_component_size):
    """
    Generates a graph with the specified number of nodes, edges, and bridges,
    with the added constraint that each 2-edge-connected component has at
    least the specified minimum size.

    Args:
        num_nodes: The number of nodes in the graph.
        num_edges: The number of edges in the graph.
        num_bridges: The number of bridges the graph should contain.
        min_component_size: The minimum number of nodes in each 2-edge-connected component.

    Returns:
        A networkx Graph object if successful, None otherwise.
    """
    if num_edges < num_nodes - 1 or num_edges > num_nodes * (num_nodes - 1) / 2:
        print("Invalid number of edges.")
        return None

    if num_bridges > num_nodes - 1 or num_bridges > num_edges:
        print("Invalid number of bridges. Cannot have more bridges than nodes - 1 or edges.")
        return None

    if min_component_size < 1:
        print("Invalid minimum component size.")
        return None

    for attempt in range(1000): # Try multiple times to generate a suitable graph
        # 1. Start with a path graph
        graph = nx.path_graph(num_nodes)
        edges_added = num_nodes - 1
        bridge_edges = []

        # 2. Add bridges, ensuring they connect to end nodes of paths
        for i in range(num_bridges):
            n1 = i
            n2 = num_nodes - 1 - i
            graph.add_edge(n1, n2)
            bridge_edges.append((n1,n2))
            edges_added += 1

        # Remove path edges to make space for bridges
        path_edges_to_remove = list(graph.edges())[:num_bridges]
        graph.remove_edges_from(path_edges_to_remove)


        # 3. Add remaining non-bridge edges randomly
        remaining_edges = num_edges - edges_added
        possible_edges = [(u, v) for u in range(num_nodes) for v in range(u + 1, num_nodes) if not graph.has_edge(u, v)]
        if remaining_edges > len(possible_edges):
            print("Not enough unique edges to fulfill request")
            return None

        edges_to_add = random.sample(possible_edges, remaining_edges)
        graph.add_edges_from(edges_to_add)



        # 4. Check 2-edge-connected component sizes
        connected_components = list(nx.biconnected_components(graph)) # Gets the 2-edge connected components
        valid = True
        for component in connected_components:
            if len(component) < min_component_size:
                valid = False
                break # Exit loop, reject this graph

        if valid:
            bridges = list(nx.bridges(graph)) # Check the number of bridges
            if len(bridges) == num_bridges:
                return graph # Found a valid graph
            else:
                print(f"Incorrect number of bridges: Expected {num_bridges}, found {len(bridges)}. Retrying...")
        else:
            print(f"A component is smaller than {min_component_size}. Retrying...")


    print(f"Failed to generate a graph meeting criteria after {attempt} attempts.")
    return None # Indicate failure after max attempts



if __name__ == "__main__":
    num_nodes = 40
    num_edges = 80
    num_bridges = 2
    min_component_size = 3
    graph = generate_graph_with_bridges_and_component_size(num_nodes, num_edges, num_bridges, min_component_size)

    if graph is not None:
        print(f"Generated a graph with {num_nodes} nodes, {num_edges} edges, {num_bridges} bridges, and minimum 2-edge-connected component size of {min_component_size}.")

        # Visualize
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(graph)
        nx.draw(graph, pos, with_labels=True, node_color='skyblue', edge_color='gray', node_size=500)
        bridges = list(nx.bridges(graph))
        nx.draw_networkx_edges(graph, pos, edgelist=bridges, edge_color='red', width=2)
        plt.title(f"Graph with {num_nodes} Nodes, {num_edges} Edges, {num_bridges} Bridges, Min Component Size {min_component_size} (Red=Bridges)")
        plt.show()
    else:
        print("Failed to generate a graph meeting the criteria.")
