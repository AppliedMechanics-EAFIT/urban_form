"""Persistent homology computation and persistence landscapes."""
import numpy as np
import gudhi
from gudhi.representations import Landscape


def compute_persistence(arrival_time_map, dt=1.0, min_h1_steps=0):
    cc = gudhi.CubicalComplex(
        dimensions=arrival_time_map.shape,
        top_dimensional_cells=arrival_time_map.flatten()
    )
    cc.persistence()
    min_h1_lifetime = min_h1_steps * dt

    pairs_h0, pairs_h1 = [], []
    for dim, (b, d) in cc.persistence():
        if dim == 0:
            if d == float('inf'):
                pairs_h0.append((b, d))
                continue
            if d - b > 0:
                pairs_h0.append((b, d))
        elif dim == 1:
            if d == float('inf'):
                continue
            lifetime = d - b
            if lifetime > 0 and lifetime > min_h1_lifetime:
                pairs_h1.append((b, d))

    print(f'H0 features: {len(pairs_h0)}, H1 features: {len(pairs_h1)}')
    return pairs_h0, pairs_h1


def compute_landscape(pairs_h1, num_landscapes, resolution):
    pairs_arr = np.array(pairs_h1, dtype=float)
    L = Landscape(num_landscapes=num_landscapes, resolution=resolution, keep_endpoints=True)
    L.fit([pairs_arr])
    vec = L.transform([pairs_arr])[0]
    landscapes = vec.reshape(L.num_landscapes, L.resolution)
    x_grid = L.grid_
    return landscapes, x_grid, vec, L
