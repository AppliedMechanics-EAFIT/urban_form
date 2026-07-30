# Urban Level-Set Persistent Homology

This folder contains experiments that use level sets to define filtrations for
persistent homology on urban street-network data. The project studies how
topological summaries from persistence diagrams and persistence landscapes can
support clustering and classification of urban forms.

The general workflow is to rasterize a street network or polygon-based urban
sample, initialize a signed distance function around the street geometry,
propagate a level-set front through open space, compute cubical persistent
homology from the resulting arrival-time map, and compare the resulting
persistence diagrams or landscapes with clustering and classification methods.

## Folder overview

- `src/`: executable pipeline code for the level-set filtration. It downloads or
  loads street networks, rasterizes streets, evolves the front, computes cubical
  persistent homology, and exports diagrams, barcodes, landscapes, and figures.
- `src/utils/`: reusable utilities for network handling, level-set evolution,
  homology computation, and visualization.

- `notebooks_clustering/`: notebooks for running the urban level-set persistent
  homology workflow and experimenting with downstream analysis. These include
  Wasserstein-distance clustering, Wasserstein-based classification,
  persistence-landscape classification, and PCA-space clustering of persistence
  landscapes.
- `notebooks_clustering/clustering_medellin/`: Medellin-specific clustering
  experiments, distance matrices, OSMnx street data, and exported cluster maps.
- `notebooks_clustering/clustering_medellin/osmnx_cache/`: OSMnx cache for the
  Medellin clustering notebooks.

- `labeled_data/`: input labels and geometries used to build classification and
  clustering datasets.
