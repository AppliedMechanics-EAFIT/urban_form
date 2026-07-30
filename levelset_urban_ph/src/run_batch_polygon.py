"""
run_batch_polygon.py — Run the persistent homology pipeline using tract polygon geometries.

Instead of fetching a street network by radius around a centroid, this script
uses the MultiPolygon geometry from each entry to query osmnx directly.

Usage:
    python src/run_batch_polygon.py
"""
import sys
import json
import warnings
import gc
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')

from shapely.geometry import shape

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = Path(__file__).resolve().parent
CENTROIDS_PATH = PROJECT_ROOT / "labeled_data" / "with_centroids" / "sample_everything.json"
OUTPUT_BASE = PROJECT_ROOT / "figs" / "test_us_cities" / "everything"

IMG_SIZE     = 1000
STREET_WIDTH = 2
V            = 1.0
DT           = 0.1
N_STEPS      = 2000
NUM_LANDSCAPES = 50
RESOLUTION     = 500
# ---------------------------------------------------------------------------


def fetch_network_polygon(polygon):
    import osmnx as ox
    G = ox.graph_from_polygon(polygon, network_type='drive', retain_all=True)
    G_proj = ox.project_graph(G)
    nodes, edges = ox.graph_to_gdfs(G_proj)
    x_min = nodes.geometry.x.min()
    y_min = nodes.geometry.y.min()
    x_max = nodes.geometry.x.max()
    y_max = nodes.geometry.y.max()
    print(f'  Nodes: {len(nodes)}, Edges: {len(edges)}')
    print(f'  Bounding box (m): Δx={x_max-x_min:.0f}, Δy={y_max-y_min:.0f}')
    return nodes, edges, x_min, x_max, y_min, y_max


def run_point(entry):
    point_id = entry["id"]
    polygon  = shape(entry["geometry"])  # shapely geometry from MultiPolygon
    save_dir = OUTPUT_BASE / point_id
    save_dir.mkdir(parents=True, exist_ok=True)

    def _path(fname):
        return str(save_dir / fname)

    from utils import (
        rasterize_streets,
        build_sdf, propagate,
        compute_persistence, compute_landscape,
        plot_street_binary, plot_sdf, plot_propagation_steps,
        plot_arrival_time_map, plot_persistence_diagrams,
        plot_barcodes, plot_landscape,
    )

    print(f"\n[{point_id}] Fetching network from polygon geometry ...")
    nodes, edges, x_min, x_max, y_min, y_max = fetch_network_polygon(polygon)

    street_binary = rasterize_streets(
        edges, nodes, IMG_SIZE, x_min, x_max, y_min, y_max, line_width=STREET_WIDTH
    )
    plot_street_binary(street_binary, save_files=True, save_path=_path('street_binary.pdf'))

    phi0 = build_sdf(street_binary)
    plot_sdf(phi0, save_files=True, save_path=_path('sdf.pdf'))

    arrival_time_map, snapshots = propagate(phi0, IMG_SIZE, DT, V, N_STEPS)
    plot_propagation_steps(snapshots, phi0, street_binary, DT,
                           save_files=True, save_path=_path('levelset_propagation_steps.pdf'))
    plot_arrival_time_map(arrival_time_map, save_files=True, save_path=_path('arrival_time_map.pdf'))

    pairs_h0, pairs_h1 = compute_persistence(arrival_time_map, dt=DT, min_h1_steps=10)
    inf_death_cap = arrival_time_map.max() * 1.05 if arrival_time_map.max() > 0 else 1.0

    plot_persistence_diagrams(pairs_h0, pairs_h1, inf_death_cap=inf_death_cap,
                              save_files=True, save_path=_path('persistence_diagrams.pdf'))
    plot_barcodes(pairs_h0, pairs_h1, inf_death_cap=inf_death_cap,
                  save_files=True, save_path=_path('persistence_barcodes.pdf'))

    if pairs_h1:
        landscapes, x_grid, flat_vec, L_obj = compute_landscape(pairs_h1, NUM_LANDSCAPES, RESOLUTION)
        plot_landscape(landscapes, x_grid, save_files=True, save_path=_path('persistence_landscape.pdf'))
        np.savetxt(_path('h1_pairs.csv'), np.array(pairs_h1), delimiter=',')
        np.savez(
            _path('persistence_landscape_h1.npz'),
            landscape_values=landscapes,
            x_grid=x_grid,
            flat_vector=flat_vec,
            num_landscapes=L_obj.num_landscapes,
            resolution=L_obj.resolution,
        )
    else:
        print(f'  [{point_id}] No H1 features found — skipping landscape.')

    del street_binary, phi0, arrival_time_map, snapshots, pairs_h0, pairs_h1
    gc.collect()

    print(f'  [{point_id}] Done ✓  →  {save_dir}')


def main():
    if str(SRC_DIR) not in sys.path:
        sys.path.insert(0, str(SRC_DIR))

    with CENTROIDS_PATH.open() as f:
        entries = json.load(f)

    print(f"Processing {len(entries)} polygons → {OUTPUT_BASE}")

    for entry in entries:
        try:
            run_point(entry)
        except Exception as e:
            print(f"  ERROR [{entry['id']}]: {e}")


if __name__ == '__main__':
    main()
