# Urban Form Research Projects

This repository collects research prototypes for studying urban form with
topological data analysis, persistent homology, graph kernels, and street-network
analysis. Most folders are self-contained case studies with their own notebooks,
data, figures, and local README files.

The common theme is to represent cities or street networks as geometric,
topological, or graph objects, then compare how urban connectivity,
accessibility, morphology, and service coverage change under different
filtrations or graph-similarity methods.

## Project folders

- `15_min_city_ph/`: proof of concept for studying 15-minute-city accessibility
  in Medellin with persistent homology. The main workflow builds a walking-time
  accessibility surface to shopping malls and retail buildings, converts it into
  a cubical-complex filtration, and studies connected accessibility islands and
  persistent service voids. The folder includes the main mall-accessibility
  notebook, exploratory service-query notebooks, cached walking-time data, a
  Medellin GraphML network, and exported figures.

- `edge_importance_ph/`: persistent-homology experiments based on road
  importance. This project computes edge-importance proxies such as edge
  betweenness centrality, uses those scores to define filtrations of road
  networks, and studies how connectivity and cycles evolve as less important
  roads are added. It also includes Corcoran-and-Jones-style road-hierarchy
  replication work, executable filtration pipelines, comparison-city graph
  exports, experiment outputs, and Medellin road-importance figures.

- `elevation_filtration_ph/`: experiments on how city connectedness changes
  under elevation and street-slope filtrations. The notebooks download or load
  OSMnx street networks, enrich them with elevation data, and build sublevel or
  superlevel filtrations by node elevation or edge slope. The folder contains
  GraphML inputs for cities such as Medellin, Manizales, La Paz, Chicago,
  Lisbon, and Funchal, plus generated figures for slope and elevation
  experiments.

- `graph_kernels/`: graph-comparison experiments around Metric Graph Kernels via
  the Tropical Torelli Map, following Cao and Monod. The folder contains the
  kernel implementation, graph dataset preparation utilities, toy graph
  notebooks, MUTAG benchmark tests, city-classification replication notebooks,
  Boston-specific urban-pattern experiments, GraKeL tutorials, and synthetic
  grid/non-grid classification tests.

- `levelset_urban_ph/`: level-set persistent-homology experiments for urban
  street-network morphology. The pipeline rasterizes street networks or polygon
  samples, initializes a signed-distance field, propagates a level-set front
  through open space, computes cubical persistent homology, and exports
  persistence diagrams, barcodes, and landscapes. The folder also contains
  clustering and classification notebooks using Wasserstein distances and
  persistence landscapes, processed outputs for Medellin and U.S. city samples,
  and labeled input geometries.

## Repository structure

Each top-level project folder is organized around some combination of:

- notebooks for exploratory analysis and reproducible case studies;
- `src/` or `pipeline/` folders for reusable implementation code;
- `data/`, `labeled_data/`, or `processed_data/` folders for inputs and derived
  outputs;
- `figs/` or `figures/` folders for exported plots, maps, persistence diagrams,
  and filtration snapshots;
- `cache/`, `osmnx_cache/`, and `__pycache__/` folders for generated
  intermediate files.

The folder-level README files provide more specific notes about methods,
notebooks, and outputs inside each project.


# License

All code is licensed under the MIT License, and the content is licensed under
the Creative Commons Attribution License.

The content in this repository is licensed under the
[Creative Commons Attribution 4.0 International License](http://choosealicense.com/licenses/cc-by-4.0/),
and the accompanying source code is licensed under the [MIT License](https://opensource.org/licenses/mit-license.php).