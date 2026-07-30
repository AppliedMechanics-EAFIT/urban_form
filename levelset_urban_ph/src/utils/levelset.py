"""Signed distance function and level-set propagation."""
import numpy as np
from scipy.ndimage import distance_transform_edt


def build_sdf(binary_img):
    street_mask = (binary_img == 0)
    open_mask = (binary_img == 1)
    dist_from_street = distance_transform_edt(open_mask)
    dist_from_open = distance_transform_edt(street_mask)
    phi = dist_from_street.astype(float) - dist_from_open.astype(float)
    return phi


def godunov_grad_mag(phi, v_sign):
    phi_p_x = np.zeros_like(phi)
    phi_n_x = np.zeros_like(phi)
    phi_p_y = np.zeros_like(phi)
    phi_n_y = np.zeros_like(phi)

    phi_p_x[:, :-1] = phi[:, 1:] - phi[:, :-1]
    phi_n_x[:, 1:] = phi[:, 1:] - phi[:, :-1]
    phi_p_y[:-1, :] = phi[1:, :] - phi[:-1, :]
    phi_n_y[1:, :] = phi[1:, :] - phi[:-1, :]

    if v_sign > 0:
        Dx = np.maximum(np.maximum(phi_n_x, 0) ** 2, np.minimum(phi_p_x, 0) ** 2)
        Dy = np.maximum(np.maximum(phi_n_y, 0) ** 2, np.minimum(phi_p_y, 0) ** 2)
    else:
        Dx = np.maximum(np.minimum(phi_n_x, 0) ** 2, np.maximum(phi_p_x, 0) ** 2)
        Dy = np.maximum(np.minimum(phi_n_y, 0) ** 2, np.maximum(phi_p_y, 0) ** 2)

    return np.sqrt(Dx + Dy)


def advance_level_set(phi, dt=0.4, v=1.0):
    grad_mag = godunov_grad_mag(phi, v_sign=np.sign(v))
    return phi - dt * v * grad_mag


def propagate(phi0, img_size, dt, v, n_steps):
    arrival_time_map = np.full((img_size, img_size), np.inf, dtype=np.float64)
    arrival_time_map[phi0 <= 0] = 0.0

    phi = phi0.copy()
    num_snapshots = min(n_steps, max(10, 12))
    snapshot_steps = np.unique(np.linspace(1, n_steps, num=num_snapshots, dtype=int)).tolist()
    snapshot_lookup = set(snapshot_steps)
    snapshots = {}

    for step in range(1, n_steps + 1):
        phi = advance_level_set(phi, dt=dt, v=v)
        newly_reached = (phi <= 0) & (arrival_time_map == np.inf)
        arrival_time_map[newly_reached] = step * dt
        if step in snapshot_lookup:
            snapshots[step] = phi.copy()

    max_time = n_steps * dt
    arrival_time_map[arrival_time_map == np.inf] = max_time

    print(f'Propagation complete ✓')
    print(f'Arrival time range: [{arrival_time_map.min():.2f}, {arrival_time_map.max():.2f}]')
    return arrival_time_map, snapshots
