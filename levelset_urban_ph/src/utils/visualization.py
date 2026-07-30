"""Visualization utilities for level-set and persistence results."""
import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'axes.facecolor': 'white',
    'figure.facecolor': 'white',
    'text.color': '#111111',
    'axes.labelcolor': '#111111',
    'xtick.color': '#111111',
    'ytick.color': '#111111',
})


def _maybe_save(fig, path, save_files):
    if save_files and path:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path, bbox_inches='tight', facecolor='white')
        print(f'Saved → {path}')
    plt.close(fig)


def plot_street_binary(street_binary, save_files=False, save_path=None):
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.imshow(street_binary, cmap='bone_r', interpolation='nearest')
    ax.set_title('Initial Street Network (binary)', color='#111111', fontsize=12)
    ax.axis('off')
    plt.tight_layout()
    _maybe_save(fig, save_path, save_files)


def plot_sdf(phi0, save_files=False, save_path=None):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    im0 = axes[0].imshow(phi0, cmap='RdBu_r', vmin=-20, vmax=20)
    axes[0].set_title('Initial SDF  φ₀', color='#111111')
    axes[0].axis('off')
    plt.colorbar(im0, ax=axes[0], fraction=0.046)
    axes[1].imshow(phi0 > 0, cmap='inferno')
    axes[1].set_title('Open space  (φ > 0)', color='#111111')
    axes[1].axis('off')
    plt.tight_layout()
    _maybe_save(fig, save_path, save_files)


def plot_propagation_steps(snapshots, phi0, street_binary, dt, save_files=False, save_path=None):
    snapshot_steps = sorted(snapshots.keys())
    n_snap = len(snapshot_steps)
    ncols = 4
    nrows = int(np.ceil(n_snap / ncols))
    max_phi = phi0.max()

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 4.2), facecolor='white')
    axes = np.atleast_1d(axes).ravel()

    for idx, step in enumerate(snapshot_steps):
        ax = axes[idx]
        phi_s = snapshots[step]
        T = step * dt
        reached = phi_s <= 0
        ax.imshow(phi0, cmap='gray', alpha=0.15, vmin=-30, vmax=80)
        fill_display = np.where(reached, -phi_s, np.nan)
        ax.imshow(fill_display, cmap='inferno', alpha=0.85, vmin=0, vmax=max(max_phi * 0.9, 1e-6))
        ax.contour(phi_s, levels=[0], colors=['#e63946'], linewidths=[0.8], alpha=0.9)
        ax.contour(street_binary, levels=[0.5], colors=['#222222'], linewidths=[0.4], alpha=0.5)
        coverage = 100.0 * reached.sum() / reached.size
        ax.set_title(f't = {T:.1f}   |   {coverage:.1f}% covered', color='#111111', fontsize=9, pad=4)
        ax.set_facecolor('white')
        ax.axis('off')

    for j in range(n_snap, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle('Level Set Propagation on Urban Street Network', color='#111111', fontsize=13, y=1.01)
    plt.tight_layout()
    _maybe_save(fig, save_path, save_files)


def plot_arrival_time_map(arrival_time_map, save_files=False, save_path=None):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    cmap = plt.cm.magma.copy()
    cmap.set_under('black')
    im = ax.imshow(arrival_time_map, cmap=cmap, vmin=0.001, interpolation='bilinear')
    plt.colorbar(im, ax=ax, label='Arrival time tau(x)', fraction=0.046)
    ax.set_title('Arrival-Time Map tau(x)', color='#111111', fontsize=13)
    ax.axis('off')
    plt.tight_layout()
    _maybe_save(fig, save_path, save_files)


def plot_persistence_diagrams(pairs_h0, pairs_h1, inf_death_cap=None, save_files=False, save_path=None):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, pairs, dim, color in [
        (axes[0], pairs_h0, 0, '#e63946'),
        (axes[1], pairs_h1, 1, '#457b9d'),
    ]:
        finite_pairs = [(b, d) for b, d in pairs if np.isfinite(d)]
        essential_pairs = [(b, d) for b, d in pairs if not np.isfinite(d)]
        finite_births = [b for b, _ in finite_pairs]
        finite_deaths = [d for _, d in finite_pairs]
        essential_births = [b for b, _ in essential_pairs]
        reference_vals = finite_births + finite_deaths + essential_births
        lo, hi = (min(reference_vals), max(reference_vals)) if reference_vals else (0.0, 1.0)
        local_cap = hi + max(0.5, 0.05 * (hi - lo if hi > lo else 1.0))
        cap = inf_death_cap if inf_death_cap is not None else local_cap
        if finite_pairs:
            ax.scatter(finite_births, finite_deaths, s=15, alpha=0.7, color=color, edgecolors='none')
        if essential_pairs:
            ax.scatter(essential_births, [cap] * len(essential_pairs), s=28, marker='^', alpha=0.9, color=color, edgecolors='none')
            ax.text(0.03, 0.97, '▲ death = ∞ (capped)', transform=ax.transAxes, va='top', fontsize=8, color='#333333')
        all_vals = reference_vals + ([cap] if essential_pairs else [])
        if not all_vals:
            all_vals = [0.0, 1.0]
        lo_plot, hi_plot = min(all_vals), max(all_vals)
        margin = (hi_plot - lo_plot) * 0.05 if hi_plot > lo_plot else 0.5
        ax.plot([lo_plot - margin, hi_plot + margin], [lo_plot - margin, hi_plot + margin], 'k--', linewidth=0.8, alpha=0.5)
        ax.set_xlim(lo_plot - margin, hi_plot + margin)
        ax.set_ylim(lo_plot - margin, hi_plot + margin)
        ax.set_xlabel('Birth time')
        ax.set_ylabel('Death time')
        ax.set_title(f'Persistence Diagram — H{dim}')
        ax.set_aspect('equal')
    plt.tight_layout()
    _maybe_save(fig, save_path, save_files)


def plot_barcodes(pairs_h0, pairs_h1, inf_death_cap=None, save_files=False, save_path=None):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, pairs, dim, color in [
        (axes[0], pairs_h0, 0, '#e63946'),
        (axes[1], pairs_h1, 1, '#457b9d'),
    ]:
        if pairs:
            sorted_pairs = sorted(pairs, key=lambda p: p[0])
            finite_vals = [v for pair in sorted_pairs for v in pair if np.isfinite(v)]
            local_cap = (max(finite_vals) + max(0.5, 0.05 * (max(finite_vals) - min(finite_vals) if max(finite_vals) > min(finite_vals) else 1.0))) if finite_vals else 1.0
            cap = inf_death_cap if inf_death_cap is not None else local_cap
            x_vals = []
            has_essential = False
            for i, (b, d) in enumerate(sorted_pairs):
                d_plot = d if np.isfinite(d) else cap
                has_essential = has_essential or (not np.isfinite(d))
                x_vals.extend([b, d_plot])
                ax.plot([b, d_plot], [i, i], color=color, linewidth=1.5, solid_capstyle='butt')
                if not np.isfinite(d):
                    ax.plot(d_plot, i, marker='>', markersize=4, color=color)
            ax.set_ylim(-1, len(sorted_pairs))
            x_lo, x_hi = min(x_vals), max(x_vals)
            margin = (x_hi - x_lo) * 0.05 if x_hi > x_lo else 0.5
            ax.set_xlim(max(0.0, x_lo - margin), x_hi + margin)
            if has_essential:
                ax.text(0.03, 0.97, '> death = ∞ (capped)', transform=ax.transAxes, va='top', fontsize=8, color='#333333')
        ax.set_xlabel('Filtration value (arrival time)')
        ax.set_ylabel('Feature index')
        ax.set_title(f'Persistence Barcode — H{dim}')
    plt.tight_layout()
    _maybe_save(fig, save_path, save_files)


def plot_landscape(landscapes, x_grid, save_files=False, save_path=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    for k in range(landscapes.shape[0]):
        ax.plot(x_grid, landscapes[k])
    ax.set_xlabel('t')
    ax.set_ylabel('Landscape value')
    ax.set_title('Persistence Landscape — H1')
    plt.tight_layout()
    _maybe_save(fig, save_path, save_files)
