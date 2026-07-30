import json
import random
from pathlib import Path
import networkx as nx

NUM_GRAPHS = 500
MIN_SIZE = 3
MAX_SIZE = 8
EDGE_REMOVAL_PROB = 0.0
OUTPUT_PATH = Path(__file__).resolve().parent / "synthetic_grid_graphs.json"


def generate_connected_grid(rows, cols):

    G = nx.Graph()

    node_id = lambda r, c: r * cols + c + 1
    w1  = random.random()
    w2 = random.random()
    # add nodes
    for r in range(rows):
        for c in range(cols):
            G.add_node(node_id(r, c))

    # add horizontal edges
    for r in range(rows):
        for c in range(cols - 1):
            u = node_id(r, c)
            v = node_id(r, c + 1)
            G.add_edge(u, v, weight=w1)

    # add vertical edges
    for r in range(rows - 1):
        for c in range(cols):
            u = node_id(r, c)
            v = node_id(r + 1, c)
            G.add_edge(u, v, weight=w2)

    # attempt random edge removals while preserving connectivity
    edges = list(G.edges())

    for u, v in edges:
        if random.random() < EDGE_REMOVAL_PROB:

            weight = G[u][v]["weight"]
            G.remove_edge(u, v)

            # if graph becomes disconnected, undo removal
            if not nx.is_connected(G):
                G.add_edge(u, v, weight=weight)

    return G


def graph_to_json(G, name):

    nodes = sorted(list(G.nodes()))

    edges = []
    for u, v, data in G.edges(data=True):

        # ensure unique ordering (avoid duplicates)
        if u > v:
            u, v = v, u

        edges.append([int(u), int(v), float(data["weight"])])

    return {
        "name": name,
        "nodes": nodes,
        "edges": edges
    }


def main():

    graphs_json = {"graphs": []}

    for i in range(NUM_GRAPHS):

        rows = random.randint(MIN_SIZE, MAX_SIZE)
        cols = random.randint(MIN_SIZE, MAX_SIZE)

        G = generate_connected_grid(rows, cols)

        name = f"G{i+1}_grid_{rows}x{cols}"

        graphs_json["graphs"].append(graph_to_json(G, name))

    with OUTPUT_PATH.open("w") as f:
        json.dump(graphs_json, f, indent=2)


if __name__ == "__main__":
    main()
