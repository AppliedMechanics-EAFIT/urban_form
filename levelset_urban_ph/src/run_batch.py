"""
run_batch.py — Run the persistent homology pipeline for all sampled centroids.

Usage:
    python src/run_batch.py
"""
import sys
import json
import warnings
import gc
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, no GUI overhead

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = Path(__file__).resolve().parent
CENTROIDS_PATH = PROJECT_ROOT / "labeled_data" / "with_centroids" / "sampled_centroids.json"
OUTPUT_BASE = PROJECT_ROOT / "figs" / "test_us_cities" / "test2"

DIST_M       = 500
IMG_SIZE     = 500
STREET_WIDTH = 2
V            = 1.0
DT           = 0.1
N_STEPS      = 1000
NUM_LANDSCAPES = 50
RESOLUTION     = 500
# ---------------------------------------------------------------------------

def run_point(entry):
    point_id   = entry["id"]
    lat, lon   = entry["coordinates"]
    center     = (lat, lon)
    save_dir   = OUTPUT_BASE / point_id
    save_dir.mkdir(parents=True, exist_ok=True)

    def _path(fname):
        return str(save_dir / fname)

    from utils import (
        fetch_network, rasterize_streets,
        build_sdf, propagate,
        compute_persistence, compute_landscape,
        plot_street_binary, plot_sdf, plot_propagation_steps,
        plot_arrival_time_map, plot_persistence_diagrams,
        plot_barcodes, plot_landscape,
    )

    print(f"\n[{point_id}] Fetching network around {center} ...")
    nodes, edges, x_min, x_max, y_min, y_max = fetch_network(center, DIST_M)

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

    # free large arrays explicitly
    del street_binary, phi0, arrival_time_map, snapshots, pairs_h0, pairs_h1
    gc.collect()

    print(f'  [{point_id}] Done ✓  →  {save_dir}')


def main():
    if str(SRC_DIR) not in sys.path:
        sys.path.insert(0, str(SRC_DIR))

    with CENTROIDS_PATH.open() as f:
        centroids = json.load(f)

    print(f"Processing {len(centroids)} points → {OUTPUT_BASE}")

    for entry in centroids:
        try:
            run_point(entry)
        except Exception as e:
            print(f"  ERROR [{entry['id']}]: {e}")


if __name__ == '__main__':
    main()
