"""Figures 13 (composite barcode) and 14 (H1 void attribution map).

Expects in the calling namespace: PH, SURF, C, KEY, LBL, mats, grid, ROOT,
SENTINEL, N_ROWS, N_COLS, city_cells, META_GREY, np, pd, plt.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import label as cc_label

CAP = 70.0


def fig13(PH, C, ROOT, META_GREY):
    """Barcode for the composite (all-three) filtration, H0 and H1."""
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.4), sharex=True)
    for ax, dim, lbl in [(axes[0], 'h0', r'$H_0$ — accessibility islands'),
                         (axes[1], 'h1', r'$H_1$ — enclosed service voids')]:
        pairs = PH['t_max'][dim]
        # longest-lived first: these are the structurally meaningful features
        ps = sorted(pairs, key=lambda x: -((min(x[1], CAP)) - x[0]))[:60]
        for i, (b, d) in enumerate(ps):
            dp = min(d, CAP)
            ax.plot([b, dp], [i, i], color=C['composite'], lw=1.5, alpha=0.85,
                    solid_capstyle='butt')
            if d >= SENTINEL:
                ax.annotate('', xy=(CAP + 1.6, i), xytext=(CAP - 0.4, i),
                            arrowprops=dict(arrowstyle='->', color=C['composite'], lw=0.9))
        ax.axvline(15, color=META_GREY, lw=1.0, ls=':')
        ax.set_ylabel(f'{lbl}\n(60 longest-lived)', fontsize=7)
        ax.set_xlim(-1, CAP + 4)
        ax.invert_yaxis()
        ax.set_yticks([])
    axes[0].annotate('15 min', xy=(15, 0), xycoords=('data', 'axes fraction'),
                     fontsize=6.5, color=META_GREY, ha='center', va='bottom',
                     xytext=(0, 3), textcoords='offset points',
                     bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.5))
    axes[1].set_xlabel('Walking-time threshold (min)')
    axes[0].set_title('Composite (all-three) barcode — bars crossing 15 min are '
                      'features the 15-minute standard does not resolve')
    fig.tight_layout()
    fig.savefig(ROOT / 'figs/fig13_barcodes_composite.png', dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig14(PH, mats, grid, C, ROOT, META_GREY, city_cells):
    """Locate the composite surface's long-lived H1 voids and attribute each
    to the service that binds it (supplies the max walking time)."""
    M = mats['t_max']
    voids = sorted([(b, d) for b, d in PH['t_max']['h1'] if d > 15 and d < SENTINEL],
                   key=lambda x: -(x[1] - x[0]))[:10]

    # rim-proxy: cells whose filtration value lies inside a void's birth-death window
    vm = np.zeros(M.shape)
    for b, d in voids:
        vm += ((M >= b) & (M <= d)).astype(float)
    vm = np.clip(vm, 0, len(voids)) / max(len(voids), 1)

    # binding-service raster
    bind = np.full(M.shape, np.nan)
    code = {'mall': 0, 'health': 1, 'school': 2}
    bind[grid.row_idx.values, grid.col_idx.values] = grid.binding.map(code).values

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.0), sharex=True, sharey=True)

    ax = axes[0]
    ax.imshow(np.where(M < SENTINEL, M, np.nan), origin='lower', cmap='magma_r',
              vmin=0, vmax=45, interpolation='nearest')
    ax.set_title('Composite accessibility\n(all three services)')

    ax = axes[1]
    ax.imshow(np.where(M <= 15, 1.0, np.nan), origin='lower', cmap='Greys',
              vmin=0, vmax=3, alpha=0.35, interpolation='nearest')
    ax.imshow(np.where(vm > 0, vm, np.nan), origin='lower', cmap='Reds',
              vmin=0, vmax=1, interpolation='nearest')
    ax.set_title(f'Top-{len(voids)} persistent $H_1$ voids\n(red = void interior / rim)')

    ax = axes[2]
    from matplotlib.colors import ListedColormap
    cmap3 = ListedColormap([C['mall'], C['health'], C['school']])
    ax.imshow(bind, origin='lower', cmap=cmap3, vmin=-0.5, vmax=2.5,
              interpolation='nearest')
    ax.set_title('Binding service\n(which service is furthest away)')
    import matplotlib.patches as mpatches
    shares = grid.binding.value_counts(normalize=True) * 100
    handles = [mpatches.Patch(color=C[k],
                              label=f'{n} ({shares.get(k, 0):.1f}% of cells)')
               for k, n in [('mall', 'Malls / retail'),
                            ('health', 'Clinics & hospitals'),
                            ('school', 'Schools')]]
    ax.legend(handles=handles, frameon=False, fontsize=6.5, loc='lower left')

    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

    fig.suptitle('Where the 15-minute city fails, and which service is responsible', y=1.0)
    fig.tight_layout()
    fig.savefig(ROOT / 'figs/fig14_void_attribution_map.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # per-void attribution table
    rows = []
    for i, (b, d) in enumerate(voids, 1):
        mask = (M >= b) & (M <= d) & (M < SENTINEL)
        rr, cc_ = np.where(mask)
        sub = grid[grid.set_index(['row_idx', 'col_idx']).index.isin(list(zip(rr, cc_)))]
        share = sub.binding.value_counts(normalize=True) * 100 if len(sub) else None
        rows.append({
            'void': i, 'birth_min': round(b, 2), 'death_min': round(d, 2),
            'persistence_min': round(d - b, 2), 'n_cells_in_window': int(mask.sum()),
            'binding_mall_pct': round(share.get('mall', 0), 1) if share is not None else np.nan,
            'binding_health_pct': round(share.get('health', 0), 1) if share is not None else np.nan,
            'binding_school_pct': round(share.get('school', 0), 1) if share is not None else np.nan,
        })
    return voids, rows
