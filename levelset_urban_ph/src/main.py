"""
main.py — Urban Level-Set Persistent Homology Pipeline

Edit the parameters below and run:
    python src/main.py
"""
import sys
import warnings
from pathlib import Path
import numpy as np

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# PARAMETERS — edit here
# ---------------------------------------------------------------------------

# Override with Boston / Roxbury config

L1 = 42.33978754488314
L2 = -71.06703730407725

CENTER_POINT = (L1, L2)
DIST_M = 1000
IMG_SIZE = 1000
STREET_WIDTH = 2

# --- Propagation ---
V = 1.0
DT = 0.1
N_STEPS = 1000

# --- Persistent homology ---
NUM_LANDSCAPES = 30
RESOLUTION = 500

# --- Output ---
SAVE_FILES = True          # set True to save all figures / data
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = Path(__file__).resolve().parent
SAVE_DIR = PROJECT_ROOT / "figs" / "boston_test"         # directory for saved files (created if needed)

# ---------------------------------------------------------------------------
# END OF PARAMETERS
# ---------------------------------------------------------------------------


def _path(filename):
    """Return full save path if SAVE_FILES is True, else None."""
    if not SAVE_FILES:
        return None
    return str(SAVE_DIR / filename)


def main():
    # Ensure src/ is on the path so `utils` package is importable
    if str(SRC_DIR) not in sys.path:
        sys.path.insert(0, str(SRC_DIR))

    from utils import (
        fetch_network, rasterize_streets,
        build_sdf, propagate,
        compute_persistence, compute_landscape,
        plot_street_binary, plot_sdf, plot_propagation_steps,
        plot_arrival_time_map, plot_persistence_diagrams,
        plot_barcodes, plot_landscape,
    )

    if SAVE_FILES:
        SAVE_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Fetch street network
    print(f'Fetching street network around {CENTER_POINT} ...')
    nodes, edges, x_min, x_max, y_min, y_max = fetch_network(CENTER_POINT, DIST_M)

    # 2. Rasterize
    street_binary = rasterize_streets(
        edges, nodes, IMG_SIZE, x_min, x_max, y_min, y_max, line_width=STREET_WIDTH
    )
    plot_street_binary(street_binary, save_files=SAVE_FILES, save_path=_path('street_binary.pdf'))

    # 3. Build signed distance function
    phi0 = build_sdf(street_binary)
    print(f'SDF range: [{phi0.min():.1f}, {phi0.max():.1f}]')
    plot_sdf(phi0, save_files=SAVE_FILES, save_path=_path('sdf.pdf'))

    # 4. Propagate level set
    arrival_time_map, snapshots = propagate(phi0, IMG_SIZE, DT, V, N_STEPS)
    plot_propagation_steps(
        snapshots, phi0, street_binary, DT,
        save_files=SAVE_FILES, save_path=_path('levelset_propagation_steps.pdf')
    )
    plot_arrival_time_map(
        arrival_time_map,
        save_files=SAVE_FILES, save_path=_path('arrival_time_map.pdf')
    )

    # 5. Persistent homology
    pairs_h0, pairs_h1 = compute_persistence(arrival_time_map, dt=DT, min_h1_steps=10)
    inf_death_cap = arrival_time_map.max() * 1.05 if arrival_time_map.max() > 0 else 1.0

    plot_persistence_diagrams(
        pairs_h0, pairs_h1, inf_death_cap=inf_death_cap,
        save_files=SAVE_FILES, save_path=_path('persistence_diagrams.pdf')
    )
    plot_barcodes(
        pairs_h0, pairs_h1, inf_death_cap=inf_death_cap,
        save_files=SAVE_FILES, save_path=_path('persistence_barcodes.pdf')
    )

    # 6. Persistence landscape (H1 only)
    if pairs_h1:
        landscapes, x_grid, flat_vec, L_obj = compute_landscape(pairs_h1, NUM_LANDSCAPES, RESOLUTION)
        plot_landscape(landscapes, x_grid, save_files=SAVE_FILES, save_path=_path('persistence_landscape.pdf'))

        if SAVE_FILES:
            np.savetxt(SAVE_DIR / 'h1_pairs.csv', np.array(pairs_h1), delimiter=',')
            np.savez(
                SAVE_DIR / 'persistence_landscape_h1.npz',
                landscape_values=landscapes,
                x_grid=x_grid,
                flat_vector=flat_vec,
                num_landscapes=L_obj.num_landscapes,
                resolution=L_obj.resolution,
            )
            print(f'Data saved to {SAVE_DIR}/')
    else:
        print('No H1 features found — skipping landscape computation.')

    print('Pipeline complete ✓')


if __name__ == '__main__':
    main()
