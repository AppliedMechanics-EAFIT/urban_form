# Metric Graph Kernels

This folder contains experiments that replicate and extend the metric graph
kernel methodology from Cao and Monod's work on Metric Graph Kernels via the
Tropical Torelli Map. The goal is to compare graphs by converting their metric
structure into kernel features and then using those features in examples such as
synthetic graph tests, urban street-network comparisons, and classification
experiments.

The general workflow is to build or load a graph dataset, compute graph features
through the tropical Torelli map kernel implementation, and evaluate whether the
resulting similarities separate graph families or urban street-network patterns.

## Folder overview

- `src/`: core implementation of the metric graph kernel described by Cao and
  Monod, together with utilities for preparing graph datasets and computing
  graph features.
- `notebooks/`: exploratory notebooks for toy graphs, MUTAG tests, Medellin
  street-network tests, appendix examples, and normalized or unweighted graph
  comparisons.
- `city_replication/`: notebooks and outputs for replicating the city
  classification experiment from the paper. The local notes indicate that the
  result is similar to the paper, but the SVM appears to overfit.
- `boston_example/`: Boston-specific experiments using urban pattern clusters,
  including tests with clustered tracts and street-network samples.
- `boston_example/utils/`: helper code for Boston subgraph sampling and local
  data preparation.
- `grakel_tutorials/`: small notebooks for learning and testing the GraKeL graph
  kernel library.
- `test_grid_no_grid/`: synthetic grid and non-grid graph experiment. It
  generates many grid-like synthetic graphs and trains a one-class SVM to test
  whether the kernel can distinguish grid structure.
