import matplotlib.pyplot as plt
import networkx as nx

# Define graph
G = nx.Graph()
G.add_nodes_from([1, 2, 3, 4, 5, 6])
G.add_edges_from([(1,2), (1,4), (2,3), (2,5), (3,4)])

# Define positions manually for a nice layout
pos = {
    1: (-1, 1),
    2: (0, 1.5),
    3: (1, 1),
    4: (0, 0),
    5: (0, -1),
    6: (2, 0),  # isolated node
}

# Draw
nx.draw(
    G, pos,
    with_labels=True,
    node_size=800,
    node_color="skyblue",
    font_size=14,
    font_weight="bold",
    edgecolors="black"
)
plt.title("Nice Layout of the Graph")
plt.show()
