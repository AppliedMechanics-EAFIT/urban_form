#!/usr/bin/env python3
"""Batch Vietoris--Rips persistent homology for Colombian cities.

The implementation follows ``persistent_homology.ipynb``.  For every city in
an input text file it runs all six equipment categories and writes results to

    results/vr_1/<city>/<equipment>/

City names are matched without regard to accents, capitalization, spaces, or
underscores.  Put one city name on each non-empty line of the input file;
lines whose first non-space character is ``#`` are ignored.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

# The pipeline is intended for unattended batch execution.
matplotlib.use("Agg")

import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.spatial import cKDTree


LOGGER = logging.getLogger("persistent_homology_pipeline")

METRIC_CRS = "EPSG:32618"
DEFAULT_MAX_FILTRATION_MIN = 60.0
N_FILTRATION_STAGES = 6
MAX_DISPLAY_EDGES = 20_000
MAX_DISPLAY_TRIANGLES = 10_000

EQUIPMENT_CODES = {
    "health": ["021001", "021002", "021003", "021004", "021011"],
    "sport": ["020601", "020602", "020603", "020604", "020605", "020606"],
    "education": ["020903", "020905"],
    "financial": ["020501", "020502"],
    "cultural": ["020701", "020702", "020706", "020709", "020710"],
    "parks": ["021101", "021102", "021104", "021105"],
}

EQUIPMENT_LABELS = {
    "health": {
        "title": "Health",
        "singular": "health center",
        "plural": "health centers",
    },
    "sport": {
        "title": "Sport",
        "singular": "sports facility",
        "plural": "sports facilities",
    },
    "education": {
        "title": "Education",
        "singular": "educational facility",
        "plural": "educational facilities",
    },
    "financial": {
        "title": "Financial",
        "singular": "financial facility",
        "plural": "financial facilities",
    },
    "cultural": {
        "title": "Cultural",
        "singular": "cultural facility",
        "plural": "cultural facilities",
    },
    "parks": {"title": "Parks", "singular": "park", "plural": "parks"},
}


@dataclass(frozen=True)
class CityNetwork:
    """Input paths and canonical names for a city pedestrian network."""

    key: str
    slug: str
    network_id: str
    edge_path: Path
    node_path: Path


@dataclass
class CityData:
    """Objects shared by all equipment categories for one city."""

    network: CityNetwork
    label: str
    edges: gpd.GeoDataFrame
    nodes: gpd.GeoDataFrame
    graph: nx.MultiGraph
    walking_graph: nx.MultiGraph
    network_nodes: gpd.GeoDataFrame
    metro_perimeter: gpd.GeoDataFrame
    blocks: gpd.GeoDataFrame


def normalize_city_name(value: object) -> str:
    """Return the city key used by the network shapefile convention."""

    ascii_value = (
        unicodedata.normalize("NFKD", str(value))
        .encode("ascii", "ignore")
        .decode("ascii")
        .casefold()
    )
    return re.sub(r"[^a-z0-9]+", "_", ascii_value).strip("_")


def find_repo_root(script_path: Path) -> Path:
    """Locate the repository from this script's directory or its parents."""

    for candidate in (script_path.parent, *script_path.parents):
        if (candidate / "data/input/shape").is_dir():
            return candidate
    raise FileNotFoundError(
        "Could not find data/input/shape in the script directory or its parents."
    )


def discover_city_networks(shape_dir: Path) -> dict[str, CityNetwork]:
    """Discover available city networks from ``*_edges.shp`` filenames."""

    networks: dict[str, CityNetwork] = {}
    for edge_path in sorted((shape_dir / "edges").glob("*_edges.shp")):
        network_stem = edge_path.stem.removesuffix("_edges")
        try:
            city_slug, network_id = network_stem.rsplit("-", maxsplit=1)
        except ValueError as exc:
            raise ValueError(
                f"Unexpected edge filename (expected city-ID_edges.shp): {edge_path}"
            ) from exc

        city_key = normalize_city_name(city_slug)
        if city_key in networks:
            raise ValueError(f"Duplicate normalized city name: {city_key}")

        node_path = shape_dir / "nodes" / f"{network_stem}_nodes.shp"
        if not node_path.exists():
            raise FileNotFoundError(f"Node file not found for {city_slug}: {node_path}")

        networks[city_key] = CityNetwork(
            key=city_key,
            slug=city_slug,
            network_id=network_id,
            edge_path=edge_path,
            node_path=node_path,
        )

    if not networks:
        raise FileNotFoundError(f"No city edge shapefiles found under {shape_dir}")
    return networks


def read_city_list(path: Path) -> list[str]:
    """Read unique city names, preserving their input order."""

    if not path.is_file():
        raise FileNotFoundError(f"City list not found: {path}")

    cities: list[str] = []
    seen: set[str] = set()
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        city = line.strip()
        if not city or city.startswith("#"):
            continue
        key = normalize_city_name(city)
        if key and key not in seen:
            cities.append(city)
            seen.add(key)

    if not cities:
        raise ValueError(f"The city list contains no city names: {path}")
    return cities


def resolve_cities(
    requested_names: Iterable[str], networks: dict[str, CityNetwork]
) -> list[CityNetwork]:
    """Resolve requested names before beginning expensive calculations."""

    resolved: list[CityNetwork] = []
    missing: list[str] = []
    for name in requested_names:
        network = networks.get(normalize_city_name(name))
        if network is None:
            missing.append(name)
        else:
            resolved.append(network)

    if missing:
        available = ", ".join(network.slug for network in networks.values())
        raise ValueError(
            f"Unknown cities: {', '.join(missing)}. Available cities: {available}"
        )
    return resolved


def apply_walking_impedance(
    edges: gpd.GeoDataFrame, speed_table: pd.DataFrame
) -> gpd.GeoDataFrame:
    """Add notebook-equivalent Naismith--Langmuir walking times to edges."""

    required_edge_columns = {"from", "to", "key", "length", "grade_abs", "geometry"}
    missing_edge_columns = required_edge_columns.difference(edges.columns)
    if missing_edge_columns:
        raise ValueError(
            "Edge shapefile is missing columns: "
            + ", ".join(sorted(missing_edge_columns))
        )

    required_speed_columns = {"s1", "s2", "speed"}
    missing_speed_columns = required_speed_columns.difference(speed_table.columns)
    if missing_speed_columns:
        raise ValueError(
            "Slope speed table is missing columns: "
            + ", ".join(sorted(missing_speed_columns))
        )

    edges = edges.copy()
    edges["from"] = edges["from"].astype(str)
    edges["to"] = edges["to"].astype(str)
    edges["speed_kmh"] = 0.0

    for interval in speed_table.itertuples(index=False):
        in_interval = (edges["grade_abs"] > interval.s1) & (
            edges["grade_abs"] <= interval.s2
        )
        edges.loc[in_interval, "speed_kmh"] = interval.speed

    edges.loc[edges["grade_abs"] == 0, "speed_kmh"] = 5.0
    edges.loc[edges["grade_abs"] > speed_table["s2"].max(), "speed_kmh"] = (
        speed_table.iloc[0]["speed"]
    )

    invalid_speed = ~np.isfinite(edges["speed_kmh"]) | (edges["speed_kmh"] <= 0)
    if invalid_speed.any():
        raise ValueError(
            f"{int(invalid_speed.sum()):,} edges have no valid walking speed."
        )

    edges["impedance_min"] = edges["length"] / (
        edges["speed_kmh"] * 1_000 / 60
    )
    return edges


def prepare_city_data(
    network: CityNetwork,
    all_perimeters: gpd.GeoDataFrame,
    all_blocks: gpd.GeoDataFrame,
    block_stats: pd.DataFrame,
    speed_table: pd.DataFrame,
) -> CityData:
    """Load and prepare all city-level geographic and network data."""

    city_label = network.slug.replace("_", " ").title()
    LOGGER.info("Preparing city %s (network %s)", city_label, network.network_id)

    edges = gpd.read_file(network.edge_path).to_crs(METRIC_CRS)
    nodes = gpd.read_file(network.node_path).to_crs(METRIC_CRS)
    nodes["osmid"] = nodes["osmid"].astype(str)
    edges = apply_walking_impedance(edges, speed_table)

    perimeter_matches = gpd.sjoin(
        edges[["geometry"]],
        all_perimeters[["ID_mpio", "Name", "geometry"]],
        how="inner",
        predicate="intersects",
    )
    if perimeter_matches.empty:
        raise ValueError(f"No urban perimeter intersects the {city_label} network.")

    metro_units = all_perimeters[
        all_perimeters["ID_mpio"].isin(perimeter_matches["ID_mpio"].unique())
    ].copy()
    metro_perimeter = gpd.GeoDataFrame(
        {"geometry": [metro_units.geometry.union_all()]}, crs=METRIC_CRS
    )

    blocks = (
        gpd.sjoin(all_blocks, metro_perimeter, how="inner", predicate="intersects")
        .drop(columns="index_right")
        .drop_duplicates(subset="MANZ_CCNCT")
        .reset_index(drop=True)
    )
    blocks["MANZ_CCNCT"] = blocks["MANZ_CCNCT"].astype("string")
    blocks = blocks.merge(
        block_stats, on="MANZ_CCNCT", how="left", validate="one_to_one"
    )

    graph = nx.from_pandas_edgelist(
        edges,
        source="from",
        target="to",
        edge_attr=[
            "length",
            "grade_abs",
            "speed_kmh",
            "impedance_min",
            "geometry",
        ],
        edge_key="key",
        create_using=nx.MultiGraph(),
    )
    node_attributes = (
        nodes.set_index("osmid")[["elevation", "geometry"]].to_dict(orient="index")
    )
    nx.set_node_attributes(graph, node_attributes)
    graph.graph["crs"] = nodes.crs

    if graph.number_of_nodes() == 0:
        raise ValueError(f"The {city_label} pedestrian graph is empty.")
    largest_component = max(nx.connected_components(graph), key=len)
    walking_graph = graph.subgraph(largest_component)
    network_nodes = (
        nodes[nodes["osmid"].isin(largest_component)]
        .drop_duplicates(subset="osmid")
        .reset_index(drop=True)
    )

    LOGGER.info(
        "%s: %s blocks, %s graph nodes, %s graph edges",
        city_label,
        f"{len(blocks):,}",
        f"{graph.number_of_nodes():,}",
        f"{graph.number_of_edges():,}",
    )
    return CityData(
        network=network,
        label=city_label,
        edges=edges,
        nodes=nodes,
        graph=graph,
        walking_graph=walking_graph,
        network_nodes=network_nodes,
        metro_perimeter=metro_perimeter,
        blocks=blocks,
    )


def select_facilities(
    all_equipment: gpd.GeoDataFrame, city: CityData, equipment: str
) -> gpd.GeoDataFrame:
    """Select one equipment category inside the city's analysis geometry."""

    selected = all_equipment[
        all_equipment["CSIMBOL"].isin(EQUIPMENT_CODES[equipment])
    ].copy()
    selected = (
        gpd.sjoin(
            selected,
            city.metro_perimeter,
            how="inner",
            predicate="intersects",
        )
        .drop(columns="index_right")
        .reset_index(drop=True)
    )
    if selected.empty:
        raise ValueError(
            f"No {EQUIPMENT_LABELS[equipment]['plural']} found in {city.label}."
        )
    return selected


def snap_facilities_to_network(
    selected_facilities: gpd.GeoDataFrame, city: CityData
) -> gpd.GeoDataFrame:
    """Snap facilities to the nearest node in the largest graph component."""

    facilities = selected_facilities[["ID", "CSIMBOL", "geometry"]].copy()
    facilities = facilities.reset_index(drop=True)
    if facilities.crs != city.network_nodes.crs:
        facilities = facilities.to_crs(city.network_nodes.crs)

    node_coordinates = np.column_stack(
        [
            city.network_nodes.geometry.x.to_numpy(),
            city.network_nodes.geometry.y.to_numpy(),
        ]
    )
    facility_coordinates = np.column_stack(
        [facilities.geometry.x.to_numpy(), facilities.geometry.y.to_numpy()]
    )
    node_tree = cKDTree(node_coordinates)
    snap_distance_m, nearest_node_position = node_tree.query(facility_coordinates, k=1)
    network_node_ids = city.network_nodes["osmid"].to_numpy()

    facilities["network_node"] = network_node_ids[nearest_node_position]
    facilities["snap_distance_m"] = snap_distance_m
    facilities["ph_vertex"] = np.arange(len(facilities), dtype=int)
    return facilities


def pairwise_travel_time_matrix(
    facilities: gpd.GeoDataFrame, walking_graph: nx.MultiGraph
) -> np.ndarray:
    """Calculate all facility-to-facility shortest walking times in minutes."""

    facility_nodes = facilities["network_node"].to_numpy()
    n_facilities = len(facilities)
    matrix = np.full((n_facilities, n_facilities), np.inf, dtype=float)
    dijkstra_cache: dict[str, dict[str, float]] = {}

    for source_index, source_node in enumerate(facility_nodes):
        if source_node not in dijkstra_cache:
            dijkstra_cache[source_node] = nx.single_source_dijkstra_path_length(
                walking_graph, source=source_node, weight="impedance_min"
            )
        shortest_times = dijkstra_cache[source_node]
        matrix[source_index] = np.fromiter(
            (shortest_times[target_node] for target_node in facility_nodes),
            dtype=float,
            count=n_facilities,
        )

    np.fill_diagonal(matrix, 0.0)
    if not np.isfinite(matrix).all():
        raise ValueError("The travel-time matrix contains unreachable pairs.")
    if not np.allclose(matrix, matrix.T):
        maximum_asymmetry = float(np.max(np.abs(matrix - matrix.T)))
        raise ValueError(
            "The travel-time matrix is asymmetric "
            f"(maximum difference {maximum_asymmetry:.6g} min)."
        )
    return matrix


def persistence_table(intervals: np.ndarray) -> pd.DataFrame:
    """Build a consistent, sortable table from persistence intervals."""

    values = np.asarray(intervals, dtype=float).reshape(-1, 2)
    table = pd.DataFrame(values, columns=["birth_min", "death_min"])
    if table.empty:
        table["persistence_min"] = pd.Series(dtype=float)
        table["essential"] = pd.Series(dtype=bool)
        return table

    table["persistence_min"] = table["death_min"] - table["birth_min"]
    table["essential"] = np.isinf(table["death_min"])
    return table.sort_values(
        ["persistence_min", "birth_min"], ascending=[False, True]
    ).reset_index(drop=True)


def save_blocks_facilities_figure(
    city: CityData,
    facilities: gpd.GeoDataFrame,
    equipment: str,
    output_path: Path,
) -> None:
    """Save the census-block and facility-location map."""

    label = EQUIPMENT_LABELS[equipment]
    fig, ax = plt.subplots(figsize=(12, 12))
    city.blocks.plot(
        ax=ax, facecolor="#eeeeee", edgecolor="#a8a8a8", linewidth=0.08
    )
    city.metro_perimeter.boundary.plot(ax=ax, color="#333333", linewidth=0.9)
    facilities.plot(
        ax=ax,
        color="#d73027",
        edgecolor="white",
        linewidth=0.3,
        markersize=16,
        alpha=0.9,
        zorder=3,
    )
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor="#d73027",
                markeredgecolor="white",
                markersize=8,
                label=f"DANE {label['singular']}",
            )
        ],
        loc="lower right",
    )
    ax.set_title(
        f"{city.label} network area: census blocks and {label['plural']}",
        fontsize=14,
    )
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def save_persistence_figures(
    persistence: list[tuple[int, tuple[float, float]]],
    city_label: str,
    equipment: str,
    output_dir: Path,
) -> None:
    """Save persistence diagram and barcode as separate PDF files."""

    try:
        import gudhi
    except ImportError as exc:
        raise RuntimeError(
            "The 'gudhi' package is required. Install the project requirements first."
        ) from exc

    plural = EQUIPMENT_LABELS[equipment]["plural"]
    diagram_fig, diagram_ax = plt.subplots(figsize=(7, 6))
    gudhi.plot_persistence_diagram(persistence, axes=diagram_ax)
    diagram_ax.set_title(
        f"{city_label} {plural}: travel-time persistence diagram"
    )
    diagram_fig.tight_layout()
    diagram_fig.savefig(
        output_dir / "persistence_diagram.pdf", format="pdf", bbox_inches="tight"
    )
    plt.close(diagram_fig)

    barcode_fig, barcode_ax = plt.subplots(figsize=(9, 6))
    gudhi.plot_persistence_barcode(persistence, axes=barcode_ax)
    barcode_ax.set_title(f"{city_label} {plural}: travel-time persistence barcode")
    barcode_fig.tight_layout()
    barcode_fig.savefig(
        output_dir / "persistence_barcode.pdf", format="pdf", bbox_inches="tight"
    )
    plt.close(barcode_fig)


def get_rips_snapshot(
    distance_matrix: np.ndarray,
    threshold: float,
    full_vr: bool = False,
) -> tuple[list[tuple[int, int]], list[tuple[int, int, int]]]:
    """Return drawable edges and triangles at a filtration threshold.

    In the default optimized mode, the returned geometry is capped to keep
    the six-panel PDF tractable.  ``full_vr=True`` returns every simplex and
    reproduces the original, potentially very expensive vector figure.
    """

    adjacency = distance_matrix <= threshold
    np.fill_diagonal(adjacency, False)
    row_indices, column_indices = np.where(np.triu(adjacency, k=1))

    if full_vr or len(row_indices) <= MAX_DISPLAY_EDGES:
        displayed_positions = np.arange(len(row_indices))
    else:
        # Evenly sample the ordered edge array so the displayed subset is not
        # concentrated among the first facility identifiers.
        displayed_positions = np.linspace(
            0,
            len(row_indices) - 1,
            MAX_DISPLAY_EDGES,
            dtype=int,
        )

    edge_pairs = list(
        zip(
            row_indices[displayed_positions].tolist(),
            column_indices[displayed_positions].tolist(),
        )
    )
    triangles: list[tuple[int, int, int]] = []

    if full_vr:
        neighbors = [set(np.flatnonzero(adjacency[v])) for v in range(len(adjacency))]
        for first, second in zip(row_indices, column_indices, strict=True):
            for third in neighbors[first] & neighbors[second]:
                if third > second:
                    triangles.append((int(first), int(second), int(third)))
    else:
        # Only discover as many triangle geometries as will actually be drawn.
        # Exact total counts come from the already-built GUDHI simplex tree.
        for first, second in zip(row_indices, column_indices, strict=True):
            third_indices = np.flatnonzero(
                adjacency[first, second + 1 :] & adjacency[second, second + 1 :]
            )
            if len(third_indices) == 0:
                continue
            remaining = MAX_DISPLAY_TRIANGLES - len(triangles)
            triangles.extend(
                (int(first), int(second), int(third + second + 1))
                for third in third_indices[:remaining]
            )
            if len(triangles) >= MAX_DISPLAY_TRIANGLES:
                break

    return edge_pairs, triangles


def simplex_counts_at_thresholds(
    simplex_tree: object, thresholds: np.ndarray
) -> np.ndarray:
    """Count exact edges and triangles at all thresholds in one tree pass."""

    count_increments = np.zeros((len(thresholds) + 1, 2), dtype=np.int64)
    for simplex, filtration_value in simplex_tree.get_filtration():
        dimension = len(simplex) - 1
        if dimension not in (1, 2):
            continue
        first_panel = int(
            np.searchsorted(thresholds, filtration_value, side="left")
        )
        if first_panel < len(thresholds):
            count_increments[first_panel, dimension - 1] += 1

    return np.cumsum(count_increments[:-1], axis=0)


def betti_at_threshold(intervals: np.ndarray, threshold: float) -> int:
    """Return the number of persistence intervals alive at ``threshold``."""

    if len(intervals) == 0:
        return 0
    values = np.asarray(intervals)
    return int(((values[:, 0] <= threshold) & (threshold < values[:, 1])).sum())


def evenly_spaced_filtration_steps(max_filtration_min: float) -> np.ndarray:
    """Return six non-zero, evenly spaced stages through the filtration."""

    return np.linspace(
        max_filtration_min / N_FILTRATION_STAGES,
        max_filtration_min,
        N_FILTRATION_STAGES,
    )


def save_filtration_evolution_figure(
    city: CityData,
    facilities: gpd.GeoDataFrame,
    equipment: str,
    travel_time_matrix: np.ndarray,
    simplex_tree: object,
    h0_intervals: np.ndarray,
    h1_intervals: np.ndarray,
    max_filtration_min: float,
    output_path: Path,
    full_vr: bool = False,
) -> None:
    """Save six evenly spaced Vietoris--Rips filtration snapshots."""

    filtration_steps = evenly_spaced_filtration_steps(max_filtration_min)
    simplex_counts = simplex_counts_at_thresholds(simplex_tree, filtration_steps)
    coordinates = np.column_stack(
        [facilities.geometry.x.to_numpy(), facilities.geometry.y.to_numpy()]
    )
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True, sharey=True)
    axes = np.asarray(axes).reshape(-1)

    xmin, ymin, xmax, ymax = city.metro_perimeter.total_bounds
    x_margin = 0.02 * (xmax - xmin)
    y_margin = 0.02 * (ymax - ymin)

    for panel_index, (ax, threshold) in enumerate(
        zip(axes, filtration_steps, strict=True)
    ):
        edge_pairs, triangles = get_rips_snapshot(
            travel_time_matrix,
            threshold,
            full_vr=full_vr,
        )
        total_edges, total_triangles = simplex_counts[panel_index]
        city.edges.plot(
            ax=ax,
            color="#c8c8c8",
            linewidth=0.25,
            alpha=0.45,
            zorder=1,
            rasterized=not full_vr,
        )
        city.metro_perimeter.boundary.plot(
            ax=ax, color="#4d4d4d", linewidth=0.8, zorder=2
        )

        if triangles:
            triangle_polygons = coordinates[np.asarray(triangles, dtype=int)]
            ax.add_collection(
                PolyCollection(
                    triangle_polygons,
                    facecolor="#fdae61",
                    edgecolor="none",
                    alpha=0.12,
                    zorder=3,
                    rasterized=not full_vr,
                )
            )
        if edge_pairs:
            edge_segments = coordinates[np.asarray(edge_pairs, dtype=int)]
            ax.add_collection(
                LineCollection(
                    edge_segments,
                    colors="#2c7bb6",
                    linewidths=0.7,
                    alpha=0.6,
                    zorder=4,
                    rasterized=not full_vr,
                )
            )

        ax.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            s=15,
            color="#d7191c",
            edgecolor="white",
            linewidth=0.25,
            zorder=5,
            rasterized=not full_vr,
        )
        beta_0 = betti_at_threshold(h0_intervals, threshold)
        beta_1 = betti_at_threshold(h1_intervals, threshold)
        edge_text = f"{total_edges:,} edges"
        triangle_text = f"{total_triangles:,} triangles"
        if len(edge_pairs) < total_edges:
            edge_text += f" ({len(edge_pairs):,} shown)"
        if len(triangles) < total_triangles:
            triangle_text += f" ({len(triangles):,} shown)"
        ax.set_title(
            f"t = {threshold:.2f} minutes\n"
            f"{edge_text} · {triangle_text}\n"
            rf"$\beta_0={beta_0}$ · $\beta_1={beta_1}$",
            fontsize=10,
        )
        ax.set_xlim(xmin - x_margin, xmax + x_margin)
        ax.set_ylim(ymin - y_margin, ymax + y_margin)
        ax.set_aspect("equal")
        ax.set_axis_off()

    label = EQUIPMENT_LABELS[equipment]
    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor="#d7191c",
                markeredgecolor="white",
                markersize=8,
                label=f"{label['singular'].capitalize()} (0-simplex)",
            ),
            Line2D(
                [0],
                [0],
                color="#2c7bb6",
                linewidth=2,
                label="Rips edge (1-simplex)",
            ),
            Patch(
                facecolor="#fdae61",
                alpha=0.3,
                label="Filled triangle (2-simplex)",
            ),
            Line2D(
                [0],
                [0],
                color="#c8c8c8",
                linewidth=2,
                label="Pedestrian network",
            ),
        ],
        loc="lower center",
        ncol=4,
        frameon=False,
    )
    fig.suptitle(
        f"Evolution of the {label['plural']} Vietoris–Rips complex\n"
        f"{city.label}: pedestrian-network travel time",
        fontsize=16,
        y=0.99,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    fig.savefig(output_path, format="pdf", bbox_inches="tight", dpi=200)
    plt.close(fig)


def run_equipment(
    city: CityData,
    all_equipment: gpd.GeoDataFrame,
    equipment: str,
    output_root: Path,
    max_filtration_min: float,
    full_vr: bool = False,
) -> None:
    """Run and save the complete pipeline for one city/equipment pair."""

    try:
        import gudhi
    except ImportError as exc:
        raise RuntimeError(
            "The 'gudhi' package is required. Install the project requirements first."
        ) from exc

    output_dir = output_root / city.network.slug / equipment
    output_dir.mkdir(parents=True, exist_ok=True)
    label = EQUIPMENT_LABELS[equipment]
    LOGGER.info("%s / %s: selecting facilities", city.label, label["title"])
    selected = select_facilities(all_equipment, city, equipment)
    facilities = snap_facilities_to_network(selected, city)
    LOGGER.info(
        "%s / %s: calculating %s × %s travel-time matrix",
        city.label,
        label["title"],
        f"{len(facilities):,}",
        f"{len(facilities):,}",
    )
    travel_time_matrix = pairwise_travel_time_matrix(facilities, city.walking_graph)

    LOGGER.info("%s / %s: building Rips complex", city.label, label["title"])
    rips_complex = gudhi.RipsComplex(
        distance_matrix=travel_time_matrix,
        max_edge_length=max_filtration_min,
    )
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
    persistence = simplex_tree.persistence(homology_coeff_field=2, min_persistence=0)
    h0_intervals = simplex_tree.persistence_intervals_in_dimension(0)
    h1_intervals = simplex_tree.persistence_intervals_in_dimension(1)

    persistence_table(h0_intervals).to_csv(
        output_dir / "h0_persistence.csv", index=False
    )
    persistence_table(h1_intervals).to_csv(
        output_dir / "h1_persistence.csv", index=False
    )
    save_blocks_facilities_figure(
        city, selected, equipment, output_dir / "blocks_facilities.pdf"
    )
    save_persistence_figures(persistence, city.label, equipment, output_dir)
    LOGGER.info(
        "%s / %s: rendering %s filtration figure",
        city.label,
        label["title"],
        "full-vector" if full_vr else "optimized",
    )
    save_filtration_evolution_figure(
        city=city,
        facilities=facilities,
        equipment=equipment,
        travel_time_matrix=travel_time_matrix,
        simplex_tree=simplex_tree,
        h0_intervals=h0_intervals,
        h1_intervals=h1_intervals,
        max_filtration_min=max_filtration_min,
        output_path=output_dir / "filtration_evolution.pdf",
        full_vr=full_vr,
    )
    LOGGER.info(
        "%s / %s complete: %s facilities, %s simplices, %s H0, %s H1",
        city.label,
        label["title"],
        f"{len(facilities):,}",
        f"{simplex_tree.num_simplices():,}",
        f"{len(h0_intervals):,}",
        f"{len(h1_intervals):,}",
    )


def load_shared_inputs(
    shape_dir: Path, table_dir: Path
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, pd.DataFrame, pd.DataFrame, gpd.GeoDataFrame]:
    """Load country-wide inputs once for the entire batch."""

    LOGGER.info("Loading shared country-wide inputs")
    all_perimeters = gpd.read_file(
        shape_dir
        / "Colombia_urban_perimeters_name_pop"
        / "Colombia_urban_perimeters_name_pop.shp"
    ).to_crs(METRIC_CRS)
    all_blocks = gpd.read_file(shape_dir / "MANZANA" / "MGN_URB_MANZANA.shp").to_crs(
        METRIC_CRS
    )
    block_stats = pd.read_csv(
        table_dir / "DANE_2018_personas_manz.txt",
        dtype={"MANZ_CCNCT": "string"},
    )
    speed_table = pd.read_csv(table_dir / "slope_velocity_data.csv", sep=";")
    all_equipment = gpd.read_file(
        shape_dir / "colombia_equipments_ID" / "colombia_equipments_ID.shp"
    ).to_crs(METRIC_CRS)
    all_equipment["CSIMBOL"] = all_equipment["CSIMBOL"].astype(str)
    return all_perimeters, all_blocks, block_stats, speed_table, all_equipment


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Run the notebook's persistent-homology pipeline for every city "
            "listed in a text file and for all six equipment categories."
        )
    )
    parser.add_argument(
        "cities_file",
        type=Path,
        help="UTF-8 text file containing one city name per line.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output root (default: results/vr_1 beside this script). "
            "City and equipment subdirectories are added automatically."
        ),
    )
    parser.add_argument(
        "--max-filtration-min",
        type=float,
        default=DEFAULT_MAX_FILTRATION_MIN,
        help=(
            "Maximum Rips filtration in minutes "
            f"(default: {DEFAULT_MAX_FILTRATION_MIN:g})."
        ),
    )
    parser.add_argument(
        "--full-vr",
        action="store_true",
        default=False,
        help=(
            "Draw every Rips edge and triangle as vector geometry. By default, "
            "dense layers are capped and rasterized for much faster PDFs."
        ),
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop at the first city/equipment failure instead of continuing.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the batch pipeline and return a process exit status."""

    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if not np.isfinite(args.max_filtration_min) or args.max_filtration_min <= 0:
        raise ValueError("--max-filtration-min must be a positive finite number.")
    try:
        import gudhi  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "The 'gudhi' package is required. Install the project requirements first."
        ) from exc

    script_path = Path(__file__).resolve()
    repo_root = find_repo_root(script_path)
    shape_dir = repo_root / "data/input/shape"
    table_dir = repo_root / "data/input/table"
    output_root = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else script_path.parent / "results" / "vr_1"
    )

    networks = discover_city_networks(shape_dir)
    requested_cities = read_city_list(args.cities_file.expanduser().resolve())
    selected_networks = resolve_cities(requested_cities, networks)
    output_root.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Results will be written to %s", output_root)
    LOGGER.info(
        "Filtration figure mode: %s",
        "full vector VR" if args.full_vr else "optimized (capped and rasterized)",
    )
    LOGGER.info(
        "Filtration stages (minutes): %s",
        np.array2string(
            evenly_spaced_filtration_steps(args.max_filtration_min), precision=2
        ),
    )

    shared = load_shared_inputs(shape_dir, table_dir)
    all_perimeters, all_blocks, block_stats, speed_table, all_equipment = shared
    failures: list[tuple[str, str, str]] = []

    for network in selected_networks:
        try:
            city = prepare_city_data(
                network, all_perimeters, all_blocks, block_stats, speed_table
            )
        except Exception as exc:
            LOGGER.exception("Could not prepare city %s", network.slug)
            failures.append((network.slug, "all", str(exc)))
            if args.fail_fast:
                raise
            continue

        for equipment in EQUIPMENT_CODES:
            try:
                run_equipment(
                    city,
                    all_equipment,
                    equipment,
                    output_root,
                    args.max_filtration_min,
                    full_vr=args.full_vr,
                )
            except Exception as exc:
                LOGGER.exception("Failed: %s / %s", network.slug, equipment)
                failures.append((network.slug, equipment, str(exc)))
                if args.fail_fast:
                    raise

    if failures:
        LOGGER.error("Pipeline finished with %s failure(s):", len(failures))
        for city_slug, equipment, message in failures:
            LOGGER.error("  %s / %s: %s", city_slug, equipment, message)
        return 1

    LOGGER.info("Pipeline completed successfully for all requested cities.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        LOGGER.error("%s", exc)
        sys.exit(2)
