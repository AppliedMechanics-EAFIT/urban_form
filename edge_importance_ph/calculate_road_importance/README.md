# mde_dijsktra

Workspace for street-network analysis with OSMnx and NetworkX, focused on Medellin and comparison cities. It stores GraphML road networks, cached OSMnx data, betweenness-centrality outputs, and figures used to study shortest paths and important road links.

## Main files

- `dijkstra_medellin.ipynb`: samples Dijkstra shortest paths on the Medellin road graph.
- `mde_edge_betweenness_dir.ipynb`: computes directed edge betweenness using travel time and exports/plots the result of proxy 'road importance' measure
- `mde_example/`, `osmnx_graphs_betweeness/`, `experiment_bounding_city/`: saved graphs and experiment outputs.
