# Urban Level-Set Pipeline

This folder contains scripts for converting OpenStreetMap street networks into
raster fields, evolving a level-set front through the rasterized urban space,
and extracting persistent homology features from the resulting arrival-time
map.

## What the Pipeline Does

1. Fetches a drivable street network with OSMnx, either around a centroid or
   inside a polygon.
2. Projects the network to metric coordinates and rasterizes street polylines
   into a binary image.
3. Builds a signed distance function where street pixels are inside the initial
   front and open-space pixels are outside it.
4. Propagates the level set with a Godunov upwind update, producing snapshots
   and an arrival-time map.
5. Runs cubical persistent homology with GUDHI on the arrival-time map.
6. Saves H0/H1 diagrams, barcodes, H1 persistence landscapes, and related data
   files.

## Main Entry Points

- `main.py`: runs the pipeline for one configured center point.
- `run_batch.py`: runs the same centroid-radius pipeline for every entry in
  `labeled_data/with_centroids/sampled_centroids.json`.
- `run_batch_polygon.py`: runs the pipeline for each polygon geometry in
  `labeled_data/with_centroids/sample_everything.json`.

## How to Use

Run commands from the project root unless you intentionally want outputs
relative to another working directory.

```bash
python src/main.py
```

To change the single-location run, edit the parameters at the top of
`main.py`, especially:

- `CENTER_POINT`
- `DIST_M`
- `IMG_SIZE`
- `STREET_WIDTH`
- `DT`, `V`, and `N_STEPS`
- `SAVE_DIR`

For batch runs, edit the input and output constants at the top of the relevant
script, then run:

```bash
python src/run_batch.py
python src/run_batch_polygon.py
```

The centroid batch input expects entries with an `id` and `coordinates`
`[latitude, longitude]`. The polygon batch input expects entries with an `id`
and GeoJSON-style `geometry`.

## Outputs

Each run writes figures and data to `SAVE_DIR` or to one directory per input
entry under `OUTPUT_BASE`:

- `street_binary.pdf`
- `sdf.pdf`
- `levelset_propagation_steps.pdf`
- `arrival_time_map.pdf`
- `persistence_diagrams.pdf`
- `persistence_barcodes.pdf`
- `persistence_landscape.pdf`
- `h1_pairs.csv`
- `persistence_landscape_h1.npz`

## Dependencies

The scripts require an environment with `numpy`, `scipy`, `scikit-image`,
`matplotlib`, `osmnx`, `gudhi`, and `shapely` installed. Network access is also
needed when OSMnx downloads street data from OpenStreetMap.
