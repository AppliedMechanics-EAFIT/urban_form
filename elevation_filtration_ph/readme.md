This folder contains different experiments done concerning how the connectedness of a city changes as its elevation or its street's slope changes.

## Notebooks

- `get_geojson_elevation.ipynb`: loads an urban boundary from `data/query.geojson`, downloads the road network inside that polygon with OSMnx, adds node elevations from OpenTopoData, checks sample node elevations, and saves the enriched graph as GraphML.
- `get_osmnx_elevation.ipynb`: downloads a city boundary directly with OSMnx/Nominatim, builds the road network for that boundary, adds node elevations from OpenTopoData, and saves the resulting elevation-enriched graph.
- `elevation_sublevel_filtration.ipynb`: loads `data/medellin_drive_with_elevation.graphml` and builds a 10-step sublevel filtration by node elevation, adding lower-elevation parts of the network first. It plots each step and summarizes connected components and cycles.
- `slope_percentiles_sublevel_filtration.ipynb`: computes edge slope as `abs(elevation_u - elevation_v) / length_m`, builds a 10-step sublevel filtration from flatter streets to steeper streets, and analyzes how components and cycles appear as slope thresholds increase.
- `slope_percentiles_superlevel_filtration.ipynb`: computes the same edge-slope measure but builds a superlevel filtration, starting with the steepest streets and progressively adding less steep streets. It is used to compare how high-slope street segments organize connectivity and cycle formation.
