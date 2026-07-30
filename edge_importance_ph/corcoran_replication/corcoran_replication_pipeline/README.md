# Corcoran Replication Pipeline

This pipeline reproduces the notebook's node-first road hierarchy filtration:

- all graph nodes are present from filtration time `0`
- edges are added cumulatively by road class
- road order is `motorway`, `trunk`, `primary`, `secondary`, `tertiary`, `unclassified`, `residential`
- H0 and H1 density CSVs always contain seven rows
- infinite deaths are replaced by `7`, matching the notebook convention
- H1 densities are reversed as in the notebook plotting cell
- `average_H0_density.csv` and `average_H1_density.csv` are saved for the processed set of cities

## Usage

Process one GraphML file:

```bash
python filtration_with_importance/corcoran_replication/corcoran_node_filtration_pipeline.py \
  filtration_with_importance/medellin_osmnx_dir.graphml \
  --out filtration_with_importance/corcoran_replication/results
```

Process every GraphML file in a folder:

```bash
python filtration_with_importance/corcoran_replication/corcoran_node_filtration_pipeline.py \
  filtration_with_importance \
  --out filtration_with_importance/corcoran_replication/results
```

Process only selected cities from a folder:

```bash
python filtration_with_importance/corcoran_replication/corcoran_node_filtration_pipeline.py \
  filtration_with_importance \
  --cities filtration_with_importance/corcoran_replication/cities.txt \
  --out filtration_with_importance/corcoran_replication/results
```

Example `cities.txt`:

```text
medellin
liverpool
vienna
```

The city filter matches normalized GraphML stems, so `medellin` matches
`medellin_osmnx_dir.graphml`.
