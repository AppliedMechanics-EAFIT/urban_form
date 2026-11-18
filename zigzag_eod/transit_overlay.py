"""Do the notebook's longest-lived H1 cycles trace Medellín's rail transit corridors?

Pulls the operator-published Metro/Metrocable/Tranvía line geometry from Medellín's
municipal GIS server (servidormapas), keeps only lines that had opened by the 2017
survey year, extracts the notebook's top-5 longest H1 bars via its own
`indicative_h1_cycles` logic, and tests cycle-to-transit proximity against a
degree-agnostic spatial null: random same-size tuples of zones drawn from the same
"ever active" pool the cycles are drawn from.

Produces:
  data/transit/transit_lines_2017.shp   -- the 7 lines open by end of 2017
  fig_transit_overlay_metro.png         -- Metro-only overlay
  fig_transit_overlay_full.png          -- Metro + Metrocable + Tranvía overlay
  transit_overlay_stats.json            -- observed vs. null proximity statistics
  transit_cycle_distances.csv           -- per-cycle distance table

Run from this directory: python transit_overlay.py
"""
import itertools
import json
import os

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import nbformat as nbf
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString

SEED = 20170101
TOP_N_BARS = 20     # bars to draw cycles from for the null comparison
TOP_K_SHOWN = 5      # distinct cycles drawn on the map
N_NULL_DRAWS = 3000

# -- lines confirmed open by the 2017 survey year (checked against operator history) --
# Excludes: Línea M (Metrocable, opened 2019), Línea P (2021), "La Aldea" cable (2021),
# Línea O (under construction in the service's own `estado` field), and Líneas 1/2
# (tipo=3, BRT bus pretroncales -- not rail).
LINES_2017 = ["Línea A", "Línea B", "Línea K", "Línea J", "Línea L", "Línea H", "Línea T"]
MODE_MAP = {1: "Metro", 2: "Metrocable", 4: "Tranvía"}


def fetch_transit_lines_2017(out_path):
    url = (
        "https://www.medellin.gov.co/servidormapas/rest/services/"
        "transporte/VC_Infraest_Gestion_Transporte/MapServer/12/query"
        "?where=1%3D1&outFields=*&returnGeometry=true&outSR=4326&f=geojson"
    )
    all_lines = gpd.read_file(url)
    gdf = all_lines[all_lines.linea.isin(LINES_2017)].copy()
    gdf["mode"] = gdf["tipo"].map(MODE_MAP)
    gdf = gdf.rename(columns={"fecha_actualizacion": "fecha_act", "st_length(shape)": "len_m"})
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    gdf.to_file(out_path)
    return gdf


def exec_notebook_through_cycles(nb_path="eod_zigzag_persistence.ipynb"):
    """Run the built notebook's cells up through H1-cycle extraction; return its namespace."""
    nb = nbf.read(nb_path, as_version=4)
    ns = {"__name__": "__main__"}
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        exec(cell.source, ns)
        if "indicative_h1_cycles(top_h1.birth_min" in cell.source:
            break
    return ns


def indicative_cycles_for_bar(row, starts, WINDOW_MIN, snapshot_graph, max_len=6, topk=2):
    """Generalizes the notebook's indicative_h1_cycles to any H1 bar row."""
    spanned = [s for s in starts if row.birth_min <= s + WINDOW_MIN / 2 <= row.death_min]
    if not spanned:
        return []
    graphs = [snapshot_graph(s) for s in spanned]
    common = set.intersection(*[{tuple(sorted(e)) for e in g.edges} for g in graphs])
    Gc = nx.Graph(list(common))
    out = [c for c in nx.chordless_cycles(Gc, length_bound=max_len) if len(c) >= 4][:topk]
    if not out:
        Gb = snapshot_graph(int(row.birth_min - WINDOW_MIN / 2))
        out = list(itertools.islice(
            (c for c in nx.chordless_cycles(Gb, length_bound=max_len) if len(c) >= 4), topk))
    return [[int(x) for x in c] for c in out]


def cycle_geom(cyc, cent):
    pts = [cent[z] for z in cyc if z in cent.index]
    if len(pts) < 3:
        return None
    xs = [p.x for p in pts] + [pts[0].x]
    ys = [p.y for p in pts] + [pts[0].y]
    return LineString(list(zip(xs, ys)))


def main():
    transit_path = "data/transit/transit_lines_2017.shp"
    if os.path.exists(transit_path):
        transit_all = gpd.read_file(transit_path)
    else:
        transit_all = fetch_transit_lines_2017(transit_path)
    print(f"transit lines (2017): {len(transit_all)}")

    ns = exec_notebook_through_cycles()
    zones, bars = ns["zones"], ns["bars"]
    starts, WINDOW_MIN, snapshot_graph = ns["starts"], ns["WINDOW_MIN"], ns["snapshot_graph"]
    cent = zones.set_index("SIT_2017").geometry.centroid
    transit_all = transit_all.to_crs(zones.crs)
    metro_only = transit_all[transit_all["mode"] == "Metro"]
    metrocable = transit_all[transit_all["mode"] == "Metrocable"]
    tranvia = transit_all[transit_all["mode"] == "Tranvía"]
    metro_union = metro_only.geometry.union_all()
    transit_union = transit_all.geometry.union_all()

    h1 = bars[bars.dim == 1].sort_values("persistence_h", ascending=False).reset_index(drop=True)
    top_bars = h1.head(TOP_N_BARS)
    all_cycles = []
    for idx, row in top_bars.iterrows():
        for c in indicative_cycles_for_bar(row, starts, WINDOW_MIN, snapshot_graph):
            all_cycles.append({"bar_rank": idx, "persistence_h": row.persistence_h, "cycle": c})
    print(f"{len(all_cycles)} indicative cycles from top {TOP_N_BARS} H1 bars")

    def zone_dists(cyc, target_union):
        ds = [cent[z].distance(target_union) for z in cyc if z in cent.index]
        return np.mean(ds) if ds else np.nan

    for r in all_cycles:
        r["n_zones"] = len(r["cycle"])
        r["mean_zone_dist_metro"] = zone_dists(r["cycle"], metro_union)
        r["mean_zone_dist_all"] = zone_dists(r["cycle"], transit_union)
    res = pd.DataFrame(all_cycles)

    # -- spatial null: random same-size zone tuples from the same "ever active" pool --
    all_active_zones = sorted(set().union(*[set(snapshot_graph(s).nodes()) for s in starts]))
    active_with_geom = [z for z in all_active_zones if z in cent.index]
    rng = np.random.default_rng(SEED)

    def random_zone_tuple_dist(k, n_draws=N_NULL_DRAWS):
        dm = np.empty(n_draws)
        da = np.empty(n_draws)
        for i in range(n_draws):
            zs = rng.choice(active_with_geom, size=k, replace=False)
            dm[i] = np.mean([cent[z].distance(metro_union) for z in zs])
            da[i] = np.mean([cent[z].distance(transit_union) for z in zs])
        return dm, da

    null_by_k = {int(k): random_zone_tuple_dist(int(k)) for k in sorted(res.n_zones.unique())}
    res["pct_metro"] = [float(np.mean(null_by_k[row.n_zones][0] <= row.mean_zone_dist_metro))
                         for _, row in res.iterrows()]
    res["pct_all"] = [float(np.mean(null_by_k[row.n_zones][1] <= row.mean_zone_dist_all))
                       for _, row in res.iterrows()]
    res.to_csv("transit_cycle_distances.csv", index=False)

    stats = {
        "n_cycles_examined": len(res),
        "observed_median_mean_zone_dist_metro_m": float(res.mean_zone_dist_metro.median()),
        "observed_median_mean_zone_dist_all_transit_m": float(res.mean_zone_dist_all.median()),
        "null_by_k_median_dist_metro_m": {k: float(np.median(v[0])) for k, v in null_by_k.items()},
        "null_by_k_median_dist_all_m": {k: float(np.median(v[1])) for k, v in null_by_k.items()},
        "frac_cycles_closer_to_metro_than_null_median": float((res.pct_metro < 0.5).mean()),
        "frac_cycles_closer_to_any_transit_than_null_median": float((res.pct_all < 0.5).mean()),
        "mean_percentile_rank_metro": float(res.pct_metro.mean()),
        "mean_percentile_rank_all_transit": float(res.pct_all.mean()),
        "note": "percentile rank = P(random same-size zone tuple has dist <= observed); "
                "lower rank = cycle sits closer to transit than a random zone tuple of the same size",
    }
    with open("transit_overlay_stats.json", "w") as f:
        json.dump(stats, f, indent=1, ensure_ascii=False)
    print(json.dumps(stats, indent=1))

    # -- figures --
    seen_ranks, picked = set(), []
    for r in all_cycles:
        if r["bar_rank"] in seen_ranks:
            continue
        seen_ranks.add(r["bar_rank"])
        picked.append(r)
        if len(picked) >= TOP_K_SHOWN:
            break
    cycle_colors = ["#C0392B", "#D68910", "#2E86C1", "#7D3C98", "#229954"]
    urban = zones[zones.geometry.area < zones.geometry.area.quantile(.75)]
    ux0, uy0, ux1, uy1 = urban.total_bounds
    upad = 0.08 * max(ux1 - ux0, uy1 - uy0)
    zone_bg = "#F7F7F7"

    def build_overlay_fig(transit_layers, title_suffix, fname):
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.6))
        ax = axes[0]
        zones.plot(ax=ax, color=zone_bg, edgecolor="white", linewidth=.15)
        for gdf_, _, color in transit_layers:
            gdf_.plot(ax=ax, color=color, linewidth=2.0, zorder=5)
        for p, c in zip(picked, cycle_colors):
            g = cycle_geom(p["cycle"], cent)
            xs, ys = g.xy
            ax.plot(xs, ys, "-", color=c, lw=1.3, alpha=.9, zorder=6)
            ax.scatter(xs[:-1], ys[:-1], s=10, color=c, zorder=7)
        ax.set_xlim(ux0 - upad, ux1 + upad)
        ax.set_ylim(uy0 - upad, uy1 + upad)
        ax.set_axis_off()
        ax.set_title(f"Urban core{title_suffix}")

        ax = axes[1]
        zone_ids = sorted({z for p in picked for z in p["cycle"]})
        sub = zones[zones.SIT_2017.isin(zone_ids)]
        bx0, by0, bx1, by1 = sub.total_bounds
        pad = 0.25 * max(bx1 - bx0, by1 - by0, 1)
        zones.plot(ax=ax, color=zone_bg, edgecolor="white", linewidth=.2)
        for gdf_, label, color in transit_layers:
            gdf_.plot(ax=ax, color=color, linewidth=2.2, zorder=5, label=label)
        for p, c in zip(picked, cycle_colors):
            g = cycle_geom(p["cycle"], cent)
            xs, ys = g.xy
            ax.plot(xs, ys, "-", color=c, lw=1.8, marker="o", ms=5, alpha=.95, zorder=6)
        ax.set_xlim(bx0 - pad, bx1 + pad)
        ax.set_ylim(by0 - pad, by1 + pad)
        ax.set_axis_off()
        ax.set_title(f"Top-{TOP_K_SHOWN} longest H$_1$ bars, zoomed{title_suffix}")

        transit_handles = [mpl.lines.Line2D([0], [0], color=color, lw=2.2, label=label)
                            for _, label, color in transit_layers]
        cycle_handles = [mpl.lines.Line2D([0], [0], color=c, lw=1.8, marker="o", ms=5,
                                           label=f"H$_1$ bar, {p['persistence_h']:.1f} h")
                          for p, c in zip(picked, cycle_colors)]
        ax.legend(handles=transit_handles + cycle_handles, loc="upper left", frameon=False, fontsize=6.3)

        fig.suptitle("Do long-lived H$_1$ cycles trace transit corridors?", fontsize=10, y=1.00)
        fig.tight_layout()
        fig.savefig(fname, dpi=200, bbox_inches="tight")
        plt.close(fig)

    build_overlay_fig([(metro_only, "Metro (Líneas A, B)", "#1A5276")],
                       " — Metro only", "fig_transit_overlay_metro.png")
    build_overlay_fig(
        [(metro_only, "Metro (A, B)", "#1A5276"),
         (metrocable, "Metrocable (H, J, K, L)", "#8E44AD"),
         (tranvia, "Tranvía (T)", "#B7950B")],
        " — full 2017 network", "fig_transit_overlay_full.png")
    print("figures saved: fig_transit_overlay_metro.png, fig_transit_overlay_full.png")


if __name__ == "__main__":
    main()
