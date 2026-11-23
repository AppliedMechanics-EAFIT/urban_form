"""Helpers for extending a city routing graph to its metro municipalities."""

from __future__ import annotations

from pathlib import Path
import re
import unicodedata

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


# Names that intentionally differ between the DANE perimeter layer and inputs.
NETWORK_SLUG_ALIASES = {
    "bogota_d_c": "bogota",
    "guadalajara_de_buga": "buga",
    "ocaa_a": "ocana",
    "la_playa": "playa_blanca",
}

# Keep the original city network when it serves at least 90% of the
# municipality's census-block centroids within a 500 m walking connector.
BASE_NETWORK_SNAP_DISTANCE_M = 500.0
MIN_BASE_NETWORK_COVERAGE_SHARE = 0.90

class NoUsableRoadsError(RuntimeError):
    """OSM returned successfully but supplied no usable roads."""


def _municipality_slug(name: object) -> str:
    """Convert a DANE municipality name to an input-network slug."""
    normalized = unicodedata.normalize("NFKD", str(name))
    normalized = normalized.encode("ascii", "ignore").decode("ascii")
    normalized = re.sub(r"[^a-z0-9]+", "_", normalized.lower())
    return re.sub(r"_+", "_", normalized).strip("_")


def _flat_speed_kmh(speed_table: pd.DataFrame) -> float:
    """Read the project's exact zero-slope speed from its source table."""
    rows = speed_table.loc[
        (speed_table["s1"] == 0) & (speed_table["s2"] == 0),
        "speed",
    ]
    if rows.empty:
        raise ValueError("The slope-speed table has no grade_abs == 0 speed.")
    return float(rows.iloc[0])


def _assign_impedance(
    edges: gpd.GeoDataFrame,
    speed_table: pd.DataFrame,
    flat_speed_kmh: float,
    *,
    force_flat_speed: bool,
) -> gpd.GeoDataFrame:
    """Attach the project impedance, using flat speed for OSM fallback roads."""
    result = edges.copy()
    result["length"] = pd.to_numeric(result["length"], errors="coerce")
    grade = pd.to_numeric(result["grade_abs"], errors="coerce").abs()
    result["grade_abs"] = 0.0 if force_flat_speed else grade.fillna(0.0)

    if force_flat_speed:
        result["speed"] = flat_speed_kmh
    else:
        result["speed"] = 0.0
        for interval in speed_table.itertuples(index=False):
            in_interval = (
                (result["grade_abs"] > interval.s1)
                & (result["grade_abs"] <= interval.s2)
            )
            result.loc[in_interval, "speed"] = interval.speed

        result.loc[result["grade_abs"] == 0, "speed"] = flat_speed_kmh
        result.loc[
            result["grade_abs"] > speed_table["s2"].max(),
            "speed",
        ] = speed_table.iloc[0]["speed"]

    result = result[
        result["length"].notna()
        & result["length"].gt(0)
        & result["speed"].gt(0)
        & result["geometry"].notna()
    ].copy()
    result["impedance_min"] = (
        result["length"] / (result["speed"] * 1_000 / 60)
    )
    result["weight"] = result["impedance_min"]
    return result


def _standardize_network(
    raw_edges: gpd.GeoDataFrame,
    raw_nodes: gpd.GeoDataFrame,
    *,
    metric_crs,
    speed_table: pd.DataFrame,
    flat_speed_kmh: float,
    source: str,
    force_flat_speed: bool,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Return the columns and CRS required by the notebook routing graph."""
    required_edge_columns = {"from", "to", "length", "geometry"}
    missing_edges = required_edge_columns.difference(raw_edges.columns)
    if missing_edges:
        raise ValueError(
            f"{source} edges are missing: {sorted(missing_edges)}"
        )

    required_node_columns = {"osmid", "geometry"}
    missing_nodes = required_node_columns.difference(raw_nodes.columns)
    if missing_nodes:
        raise ValueError(
            f"{source} nodes are missing: {sorted(missing_nodes)}"
        )

    edges = raw_edges.to_crs(metric_crs).copy()
    nodes = raw_nodes.to_crs(metric_crs).copy()

    if "key" not in edges.columns:
        edges["key"] = 0
    if "grade_abs" not in edges.columns:
        edges["grade_abs"] = 0.0
    if "elevation" not in nodes.columns:
        nodes["elevation"] = 0.0

    edges = edges[
        ["from", "to", "key", "length", "grade_abs", "geometry"]
    ].copy()
    nodes = nodes[["osmid", "elevation", "geometry"]].copy()

    edges["from"] = edges["from"].astype(str)
    edges["to"] = edges["to"].astype(str)
    nodes["osmid"] = nodes["osmid"].astype(str)
    edges["network_source"] = source
    nodes["network_source"] = source

    edges = _assign_impedance(
        edges,
        speed_table,
        flat_speed_kmh,
        force_flat_speed=force_flat_speed,
    )
    nodes = nodes[nodes["geometry"].notna()].copy()

    if edges.empty or nodes.empty:
        raise ValueError(f"{source} has no usable road edges or nodes.")

    return edges, nodes


def _load_input_network(
    edge_file: Path,
    *,
    shape_dir: Path,
    metric_crs,
    speed_table: pd.DataFrame,
    flat_speed_kmh: float,
    source: str,
    municipality_geometry,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Read an existing municipality network with its project slope values."""
    network_stem = edge_file.stem.removesuffix("_edges")
    node_file = shape_dir / "nodes" / f"{network_stem}_nodes.shp"
    if not node_file.exists():
        raise FileNotFoundError(
            f"Missing nodes for {edge_file.name}: {node_file}"
        )

    edges, nodes = _standardize_network(
        gpd.read_file(edge_file),
        gpd.read_file(node_file),
        metric_crs=metric_crs,
        speed_table=speed_table,
        flat_speed_kmh=flat_speed_kmh,
        source=source,
        force_flat_speed=False,
    )
    if not edges.geometry.intersects(municipality_geometry).any():
        raise ValueError(
            f"{source} does not intersect the target municipality perimeter."
        )
    return edges, nodes


def _download_osmnx_network(
    municipality_geometry,
    *,
    metric_crs,
    speed_table: pd.DataFrame,
    flat_speed_kmh: float,
    source: str,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Download roads for a DANE urban perimeter through OSMnx."""
    try:
        import osmnx as ox
    except ImportError as error:
        raise RuntimeError(
            "OSMnx is required to fetch an uncovered municipality."
        ) from error

    geometry_wgs84 = gpd.GeoSeries(
        [municipality_geometry],
        crs=metric_crs,
    ).to_crs("EPSG:4326").iloc[0]

    polygon_parts = (
        list(geometry_wgs84.geoms)
        if geometry_wgs84.geom_type == "MultiPolygon"
        else [geometry_wgs84]
    )

    node_frames = []
    edge_frames = []

    for polygon_part in polygon_parts:
        try:
            graph = ox.graph_from_polygon(
                polygon_part,
                network_type="drive",
                simplify=True,
                retain_all=True,
                truncate_by_edge=True,
            )
            if graph.number_of_edges() == 0:
                continue

            nodes, edges = ox.graph_to_gdfs(
                graph,
                nodes=True,
                edges=True,
                fill_edge_geometry=True,
            )
            nodes = nodes.reset_index()
            edges = edges.reset_index().rename(
                columns={"u": "from", "v": "to"}
            )

            if "osmid" not in nodes.columns:
                raise ValueError("OSMnx nodes did not expose an osmid column.")

            node_frames.append(nodes)
            edge_frames.append(edges)
        except Exception as error:
            no_road_error_names = {
                "InsufficientResponseError",
                "EmptyOverpassResponse",
                "NoDataError",
            }
            no_road_messages = (
                "no graph nodes",
                "no street network",
                "no data",
                "no features",
                "no results",
            )
            message = str(error).lower()
            if (
                type(error).__name__ in no_road_error_names
                or any(text in message for text in no_road_messages)
            ):
                continue
            raise RuntimeError(
                f"OSMnx could not download {source}: {type(error).__name__}: "
                f"{error}"
            ) from error

    if not edge_frames:
        raise NoUsableRoadsError("OSM returned no usable road edges.")

    downloaded_nodes = gpd.GeoDataFrame(
        pd.concat(node_frames, ignore_index=True),
        geometry="geometry",
        crs="EPSG:4326",
    )
    downloaded_edges = gpd.GeoDataFrame(
        pd.concat(edge_frames, ignore_index=True),
        geometry="geometry",
        crs="EPSG:4326",
    )

    return _standardize_network(
        downloaded_edges,
        downloaded_nodes,
        metric_crs=metric_crs,
        speed_table=speed_table,
        flat_speed_kmh=flat_speed_kmh,
        source=source,
        force_flat_speed=True,
    )


def _main_municipality_id(
    base_edges: gpd.GeoDataFrame,
    urban_units: gpd.GeoDataFrame,
) -> str:
    """Choose the municipality with the largest road length in the base network."""
    hits = gpd.sjoin(
        base_edges[["length", "geometry"]],
        urban_units[["ID_mpio", "geometry"]],
        how="inner",
        predicate="intersects",
    )
    if hits.empty:
        raise ValueError("Could not identify the selected city's municipality.")
    return str(hits.groupby("ID_mpio")["length"].sum().idxmax())


def _municipality_geometries(
    urban_units: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """Dissolve the original DANE urban units by municipality."""
    records = []
    for (municipality_id, municipality_name), units in urban_units.groupby(
        ["ID_mpio", "Name"],
        dropna=False,
    ):
        records.append(
            {
                "ID_mpio": str(municipality_id),
                "Name": municipality_name,
                "geometry": units.geometry.union_all(),
            }
        )
    return gpd.GeoDataFrame(records, geometry="geometry", crs=urban_units.crs)


def _base_network_coverage(
    municipality_geometry,
    *,
    shape_dir: Path,
    metric_crs,
    node_tree: cKDTree,
) -> dict[str, float | int | bool]:
    """Measure whether the original-city network serves local census blocks."""
    block_path = shape_dir / "MANZANA" / "MGN_URB_MANZANA.shp"
    municipality_wgs84 = gpd.GeoSeries(
        [municipality_geometry],
        crs=metric_crs,
    ).to_crs("EPSG:4326").iloc[0]
    blocks = gpd.read_file(block_path, mask=municipality_wgs84).to_crs(
        metric_crs
    )
    blocks = blocks[blocks.geometry.intersects(municipality_geometry)].copy()

    if blocks.empty:
        # There are no origins to route from in this unit, so no supplemental
        # network is needed for the current block-based analysis.
        return {
            "block_count": 0,
            "covered_share": 1.0,
            "p95_snap_distance_m": 0.0,
            "is_adequately_covered": True,
        }

    centroids = blocks.geometry.centroid
    distances, _ = node_tree.query(
        np.column_stack([centroids.x.to_numpy(), centroids.y.to_numpy()]),
        k=1,
    )
    covered_share = float(
        np.mean(distances <= BASE_NETWORK_SNAP_DISTANCE_M)
    )
    return {
        "block_count": int(len(blocks)),
        "covered_share": covered_share,
        "p95_snap_distance_m": float(np.quantile(distances, 0.95)),
        "is_adequately_covered": bool(
            covered_share >= MIN_BASE_NETWORK_COVERAGE_SHARE
        ),
    }


def extend_metro_network(
    *,
    city: str,
    network_id: str,
    shape_dir: Path,
    metric_crs,
    base_edges: gpd.GeoDataFrame,
    base_nodes: gpd.GeoDataFrame,
    base_city_urban_units: gpd.GeoDataFrame,
    speed_table: pd.DataFrame,
) -> tuple[
    gpd.GeoDataFrame,
    gpd.GeoDataFrame,
    gpd.GeoDataFrame,
    gpd.GeoDataFrame,
    pd.DataFrame,
    float,
]:
    """Append municipal networks without expanding the base-derived metro boundary.

    The starting urban units must be computed from the selected city's network
    before this function is called. Each other municipality receives a local
    input network or OSMnx only when the original network does not adequately
    serve its census blocks. Municipalities with no usable supplemental roads
    are excluded from both the returned polygon and later block/equipment
    selection.
    """
    flat_speed_kmh = _flat_speed_kmh(speed_table)
    base_units = base_city_urban_units.copy()
    primary_municipality_id = _main_municipality_id(base_edges, base_units)
    municipalities = _municipality_geometries(base_units)

    edge_catalog = {
        edge_file.stem.removesuffix("_edges").rsplit("-", 1)[0]: edge_file
        for edge_file in (shape_dir / "edges").glob("*_edges.shp")
    }

    main_edges, main_nodes = _standardize_network(
        base_edges,
        base_nodes,
        metric_crs=metric_crs,
        speed_table=speed_table,
        flat_speed_kmh=flat_speed_kmh,
        source=f"input:{city}-{network_id}",
        force_flat_speed=False,
    )
    edge_frames = [main_edges]
    node_frames = [main_nodes]
    included_ids = {primary_municipality_id}
    report_rows = []
    base_node_tree = cKDTree(
        np.column_stack(
            [
                main_nodes.geometry.x.to_numpy(),
                main_nodes.geometry.y.to_numpy(),
            ]
        )
    )

    print(
        "[metro] Checking original-network coverage for "
        f"{len(municipalities):,} municipalities "
        f"(>= {MIN_BASE_NETWORK_COVERAGE_SHARE:.0%} of block centroids "
        f"within {BASE_NETWORK_SNAP_DISTANCE_M:g} m keeps the base network)."
    )

    for municipality in municipalities.itertuples(index=False):
        municipality_id = str(municipality.ID_mpio)
        municipality_name = municipality.Name
        raw_slug = _municipality_slug(municipality_name)
        source_slug = NETWORK_SLUG_ALIASES.get(raw_slug, raw_slug)

        if municipality_id == primary_municipality_id:
            report_rows.append(
                {
                    "ID_mpio": municipality_id,
                    "Name": municipality_name,
                    "network_source": f"input:{city}-{network_id}",
                    "status": "included_primary_network",
                    "edges": len(main_edges),
                    "nodes": len(main_nodes),
                    "base_block_count": pd.NA,
                    "base_coverage_share": np.nan,
                    "base_p95_snap_m": np.nan,
                    "detail": (
                        "Boundary and primary network selected from this city."
                    ),
                }
            )
            continue

        coverage = _base_network_coverage(
            municipality.geometry,
            shape_dir=shape_dir,
            metric_crs=metric_crs,
            node_tree=base_node_tree,
        )
        coverage_message = (
            f"{coverage['covered_share']:.1%} of "
            f"{coverage['block_count']:,} block centroids within "
            f"{BASE_NETWORK_SNAP_DISTANCE_M:g} m "
            f"(P95 {coverage['p95_snap_distance_m']:.1f} m)"
        )

        if coverage["is_adequately_covered"]:
            print(
                f"[metro] {municipality_name} ({municipality_id}): "
                f"{coverage_message}; using the original city network."
            )
            included_ids.add(municipality_id)
            report_rows.append(
                {
                    "ID_mpio": municipality_id,
                    "Name": municipality_name,
                    "network_source": f"input:{city}-{network_id}",
                    "status": "included_base_coverage",
                    "edges": 0,
                    "nodes": 0,
                    "base_block_count": coverage["block_count"],
                    "base_coverage_share": coverage["covered_share"],
                    "base_p95_snap_m": coverage["p95_snap_distance_m"],
                    "detail": "Original city network is sufficiently close.",
                }
            )
            continue

        print(
            f"[metro] {municipality_name} ({municipality_id}): "
            f"{coverage_message}; supplemental roads required."
        )
        edges = None
        nodes = None
        source = None
        input_error = None

        if source_slug in edge_catalog:
            source = (
                f"input:{edge_catalog[source_slug].stem.removesuffix('_edges')}"
            )
            try:
                print(f"[metro] Loading local network: {source}")
                edges, nodes = _load_input_network(
                    edge_catalog[source_slug],
                    shape_dir=shape_dir,
                    metric_crs=metric_crs,
                    speed_table=speed_table,
                    flat_speed_kmh=flat_speed_kmh,
                    source=source,
                    municipality_geometry=municipality.geometry,
                )
            except Exception as error:
                input_error = f"{type(error).__name__}: {error}"

        if edges is None:
            source = f"osmnx:{source_slug}"
            try:
                print(f"[metro] Downloading roads with OSMnx: {source}")
                edges, nodes = _download_osmnx_network(
                    municipality.geometry,
                    metric_crs=metric_crs,
                    speed_table=speed_table,
                    flat_speed_kmh=flat_speed_kmh,
                    source=source,
                )
            except NoUsableRoadsError as error:
                details = [f"{type(error).__name__}: {error}"]
                if input_error is not None:
                    details.insert(0, f"Local input failed: {input_error}")
                report_rows.append(
                    {
                        "ID_mpio": municipality_id,
                        "Name": municipality_name,
                        "network_source": source,
                        "status": "excluded_no_usable_roads",
                    "edges": 0,
                    "nodes": 0,
                    "base_block_count": coverage["block_count"],
                    "base_coverage_share": coverage["covered_share"],
                    "base_p95_snap_m": coverage["p95_snap_distance_m"],
                    "detail": " | ".join(details),
                }
                )
                print(
                    f"[metro] {municipality_name} ({municipality_id}): "
                    "no usable supplemental roads; excluded."
                )
                continue
            except Exception as error:
                raise RuntimeError(
                    f"Could not establish road coverage for "
                    f"{municipality_name} ({municipality_id}); it was not "
                    f"excluded because the OSMnx request failed."
                ) from error

        edge_frames.append(edges)
        node_frames.append(nodes)
        included_ids.add(municipality_id)
        report_rows.append(
            {
                "ID_mpio": municipality_id,
                "Name": municipality_name,
                "network_source": source,
                "status": (
                    "included_input_network"
                    if source.startswith("input:")
                    else "included_osmnx_network"
                ),
                "edges": len(edges),
                "nodes": len(nodes),
                "base_block_count": coverage["block_count"],
                "base_coverage_share": coverage["covered_share"],
                "base_p95_snap_m": coverage["p95_snap_distance_m"],
                "detail": (
                    "Existing municipal input network."
                    if source.startswith("input:")
                    else (
                        "Downloaded through OSMnx; all edges use "
                        f"{flat_speed_kmh:g} km/h."
                    )
                ),
            }
        )
        print(
            f"[metro] {municipality_name} ({municipality_id}): "
            f"added {len(edges):,} edges and {len(nodes):,} nodes from {source}."
        )

    combined_edges = gpd.GeoDataFrame(
        pd.concat(edge_frames, ignore_index=True),
        geometry="geometry",
        crs=metric_crs,
    )
    combined_nodes = gpd.GeoDataFrame(
        pd.concat(node_frames, ignore_index=True),
        geometry="geometry",
        crs=metric_crs,
    ).drop_duplicates(subset="osmid", keep="first").reset_index(drop=True)

    included_units = base_units[
        base_units["ID_mpio"].astype(str).isin(included_ids)
    ].copy()
    analysis_polygon = gpd.GeoDataFrame(
        {
            "city": [city],
            "network_id": [network_id],
            "geometry": [included_units.geometry.union_all()],
        },
        crs=metric_crs,
    )
    excluded_count = len(municipalities) - len(included_ids)
    print(
        f"[metro] Finished: {len(included_ids):,} municipalities retained, "
        f"{excluded_count:,} excluded; {len(combined_edges):,} edges and "
        f"{len(combined_nodes):,} nodes available for routing."
    )

    return (
        combined_edges,
        combined_nodes,
        included_units,
        analysis_polygon,
        pd.DataFrame(report_rows),
        flat_speed_kmh,
    )
