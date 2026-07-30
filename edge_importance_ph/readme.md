# Edge Importance Filtration

This folder contains the edge importance filtration work for persistent homology.
It was developed as a case study for measuring road-network connectivity with
persistent homology by filtering roads according to their importance.

The general workflow is to compute an importance score for road edges, use that
score to define a filtration, and then study how connected components, cycles,
and persistent-density summaries evolve as increasingly less important roads are
added back into the network.

## Folder overview

- `calculate_road_importance/`: notebooks, GraphML files, and figures used to
  compute road-importance proxies, mainly edge betweenness centrality and
  shortest-path examples with OSMnx and NetworkX. It includes Medellin examples,
  comparison-city graph exports, and a bounding-city experiment workspace.
- `calculate_road_importance/mde_example/`: saved Medellin graph variants used
  as local inputs for road-importance calculations.

- `filtration_with_importance/`: implementations of persistent homology
  filtrations driven by the road-importance scores computed in
  `calculate_road_importance/`. The notebooks and scripts compare log/binning,
  decile, and Gaussian-smoothed betweenness filtrations.
- `filtration_with_importance/pipeline/`: Python pipeline scripts for running
  the road-importance filtrations and exporting persistent-density outputs.
- `filtration_with_importance/notebooks/`: exploratory Medellin notebooks that
  match the pipeline variants.

- `corcoran_replication/`: notebooks and scripts that replicate or closely
  follow the Corcoran and Jones road-hierarchy filtration methodology. These
  experiments provide a baseline for comparing road-class filtrations with
  edge-importance filtrations.
- `corcoran_replication/corcoran_replication_pipeline/`: executable pipeline
  version of the Corcoran-style node-first road hierarchy filtration.
