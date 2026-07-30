"""
medellin_levelset.py - Batch level-set persistent homology for Medellin SIT polygons.

The script reads Medellin polygon geometries from the local shapefile, queries
OSMnx for the drivable street network inside each polygon, then runs the same
street-indicator, level-set, persistent homology, and persistence landscape
pipeline used by the existing batch scripts.

Run from anywhere:
    python src/medellin_levelset.py
    python src/medellin_levelset.py --limit 5
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

# Keep matplotlib/font caches out of unwritable home directories in sandboxed runs.
_TMP_CACHE = Path(tempfile.gettempdir()) / "levelset_urban_cache"
os.environ.setdefault("MPLCONFIGDIR", str(_TMP_CACHE / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_TMP_CACHE / "xdg"))

import matplotlib

matplotlib.use("Agg")
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"

DEFAULT_SHAPEFILE = (
    PROJECT_ROOT
    / "labeled_data"
    / "medellin_polygons"
    / "shape"
    / "AMVA_EOD_2017_poligonos.shp"
)
DEFAULT_OUTPUT_BASE = PROJECT_ROOT / "processed_data" / "medellin_sit"
DEFAULT_CACHE_FOLDER = PROJECT_ROOT / "cache"

EXPECTED_OUTPUTS = (
    "street_binary.pdf",
    "sdf.pdf",
    "levelset_propagation_steps.pdf",
    "arrival_time_map.pdf",
    "persistence_diagrams.pdf",
    "persistence_barcodes.pdf",
    "persistence_landscape.pdf",
    "h1_pairs.csv",
    "persistence_landscape_h1.npz",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the urban level-set PH pipeline for Medellin SIT polygons."
    )
    parser.add_argument(
        "--shapefile",
        type=Path,
        default=DEFAULT_SHAPEFILE,
        help="Input polygon shapefile. Relative paths are resolved from the project root.",
    )
    parser.add_argument(
        "--output-base",
        type=Path,
        default=DEFAULT_OUTPUT_BASE,
        help="Directory where polygon_N output folders will be written.",
    )
    parser.add_argument(
        "--network-type",
        default="drive",
        help="OSMnx network type passed to graph_from_polygon.",
    )
    parser.add_argument("--img-size", type=int, default=1000)
    parser.add_argument("--street-width", type=int, default=2)
    parser.add_argument("--v", type=float, default=1.0)
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--n-steps", type=int, default=2000)
    parser.add_argument("--num-landscapes", type=int, default=50)
    parser.add_argument("--resolution", type=int, default=500)
    parser.add_argument("--min-h1-steps", type=int, default=10)
    parser.add_argument(
        "--start-index",
        type=int,
        default=1,
        help="1-based polygon number to start from.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of polygons to process. Omit to process all selected polygons.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Reprocess polygons even when their status.json is already complete.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Load the shapefile and print selected polygon folders without querying OSMnx.",
    )
    parser.add_argument(
        "--truncate-by-edge",
        action="store_true",
        help="Pass truncate_by_edge=True to OSMnx graph_from_polygon.",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable OSMnx HTTP caching.",
    )
    parser.add_argument(
        "--requests-timeout",
        type=int,
        default=180,
        help="OSMnx requests timeout in seconds.",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop the batch on the first polygon failure.",
    )
    return parser.parse_args()


def resolve_project_path(path: Path) -> Path:
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def validate_args(args: argparse.Namespace) -> None:
    positive_ints = {
        "img_size": args.img_size,
        "street_width": args.street_width,
        "n_steps": args.n_steps,
        "num_landscapes": args.num_landscapes,
        "resolution": args.resolution,
        "start_index": args.start_index,
    }
    for name, value in positive_ints.items():
        if value <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be positive.")

    if args.dt <= 0:
        raise ValueError("--dt must be positive.")
    if args.v == 0:
        raise ValueError("--v must be nonzero.")
    if args.min_h1_steps < 0:
        raise ValueError("--min-h1-steps must be nonnegative.")
    if args.limit is not None and args.limit < 0:
        raise ValueError("--limit must be nonnegative.")


def load_polygons(shapefile: Path):
    try:
        import geopandas as gpd
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Missing dependency 'geopandas'. Install it in the same environment "
            "used to run this script."
        ) from exc

    if not shapefile.exists():
        raise FileNotFoundError(f"Input shapefile not found: {shapefile}")

    gdf = gpd.read_file(shapefile)
    if "geometry" not in gdf:
        raise ValueError(f"No geometry column found in {shapefile}")
    if gdf.crs is None:
        raise ValueError(f"{shapefile} has no CRS. OSMnx requires EPSG:4326 polygons.")
    if gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)

    usable = gdf.geometry.notna() & (~gdf.geometry.is_empty)
    dropped = int((~usable).sum())
    if dropped:
        print(f"Dropping {dropped} rows with empty geometry.")
    gdf = gdf.loc[usable].copy()
    if gdf.empty:
        raise ValueError(f"No usable polygon geometries found in {shapefile}")
    return gdf


def import_pipeline_utils() -> dict[str, Any]:
    if str(SRC_DIR) not in sys.path:
        sys.path.insert(0, str(SRC_DIR))

    try:
        from utils import (
            build_sdf,
            compute_landscape,
            compute_persistence,
            plot_arrival_time_map,
            plot_barcodes,
            plot_landscape,
            plot_persistence_diagrams,
            plot_propagation_steps,
            plot_sdf,
            plot_street_binary,
            propagate,
            rasterize_streets,
        )
    except ModuleNotFoundError as exc:
        missing = exc.name or "unknown"
        raise RuntimeError(
            f"Missing dependency '{missing}'. The full pipeline needs the same "
            "environment as the existing persistent-homology scripts, including gudhi."
        ) from exc

    return {
        "build_sdf": build_sdf,
        "compute_landscape": compute_landscape,
        "compute_persistence": compute_persistence,
        "plot_arrival_time_map": plot_arrival_time_map,
        "plot_barcodes": plot_barcodes,
        "plot_landscape": plot_landscape,
        "plot_persistence_diagrams": plot_persistence_diagrams,
        "plot_propagation_steps": plot_propagation_steps,
        "plot_sdf": plot_sdf,
        "plot_street_binary": plot_street_binary,
        "propagate": propagate,
        "rasterize_streets": rasterize_streets,
    }


def configure_osmnx(args: argparse.Namespace):
    try:
        import osmnx as ox
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Missing dependency 'osmnx'. Install it in the same environment "
            "used to run this script."
        ) from exc

    ox.settings.use_cache = not args.no_cache
    ox.settings.cache_folder = str(DEFAULT_CACHE_FOLDER)
    ox.settings.requests_timeout = args.requests_timeout
    ox.settings.log_console = False
    return ox


def jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(jsonable(payload), f, indent=2, sort_keys=True)
        f.write("\n")


def clean_polygon(geometry: Any) -> Any:
    if geometry is None or geometry.is_empty:
        raise ValueError("Geometry is empty.")
    if not geometry.is_valid:
        geometry = geometry.buffer(0)
    if geometry.is_empty:
        raise ValueError("Geometry became empty after validity repair.")
    if geometry.geom_type == "GeometryCollection":
        from shapely.ops import unary_union

        polygons = [
            part
            for part in geometry.geoms
            if part.geom_type in {"Polygon", "MultiPolygon"} and not part.is_empty
        ]
        if not polygons:
            raise ValueError("GeometryCollection contains no polygonal parts.")
        geometry = unary_union(polygons)
    if geometry.geom_type not in {"Polygon", "MultiPolygon"}:
        raise ValueError(f"Expected Polygon or MultiPolygon, got {geometry.geom_type}.")
    return geometry


def fetch_network_polygon(
    ox: Any,
    polygon: Any,
    network_type: str,
    truncate_by_edge: bool,
):
    graph = ox.graph_from_polygon(
        polygon,
        network_type=network_type,
        retain_all=True,
        truncate_by_edge=truncate_by_edge,
    )
    if len(graph.nodes) == 0 or len(graph.edges) == 0:
        raise ValueError("OSMnx returned an empty graph.")

    graph_proj = ox.project_graph(graph)
    nodes, edges = ox.graph_to_gdfs(graph_proj)
    if nodes.empty or edges.empty:
        raise ValueError("Projected graph has no nodes or edges.")

    x_min = float(nodes.geometry.x.min())
    x_max = float(nodes.geometry.x.max())
    y_min = float(nodes.geometry.y.min())
    y_max = float(nodes.geometry.y.max())
    if x_min == x_max or y_min == y_max:
        raise ValueError("Projected graph has a degenerate bounding box.")

    print(f"  Nodes: {len(nodes)}, Edges: {len(edges)}")
    print(f"  Bounding box (m): dx={x_max - x_min:.0f}, dy={y_max - y_min:.0f}")
    return nodes, edges, x_min, x_max, y_min, y_max


def selected_rows(gdf: Any, start_index: int, limit: int | None):
    start = start_index - 1
    stop = None if limit is None else start + limit
    for polygon_number, (_, row) in enumerate(
        gdf.iloc[start:stop].iterrows(),
        start=start_index,
    ):
        yield polygon_number, row


def row_properties(row: Any) -> dict[str, Any]:
    return {str(key): row[key] for key in row.index if key != "geometry"}


def save_empty_landscape(
    save_dir: Path,
    inf_death_cap: float,
    num_landscapes: int,
    resolution: int,
    plot_landscape: Any,
) -> None:
    x_grid = np.linspace(0.0, inf_death_cap, resolution)
    landscapes = np.zeros((num_landscapes, resolution), dtype=float)
    flat_vec = landscapes.ravel()

    plot_landscape(
        landscapes,
        x_grid,
        save_files=True,
        save_path=str(save_dir / "persistence_landscape.pdf"),
    )
    np.savez(
        save_dir / "persistence_landscape_h1.npz",
        landscape_values=landscapes,
        x_grid=x_grid,
        flat_vector=flat_vec,
        num_landscapes=num_landscapes,
        resolution=resolution,
        empty=True,
    )


def already_complete(save_dir: Path) -> bool:
    status_path = save_dir / "status.json"
    if not status_path.exists():
        return False
    try:
        with status_path.open(encoding="utf-8") as f:
            status = json.load(f)
    except json.JSONDecodeError:
        return False
    if status.get("status") != "complete":
        return False
    return all((save_dir / name).exists() for name in EXPECTED_OUTPUTS)


def run_polygon(
    *,
    ox: Any,
    utils: dict[str, Any],
    args: argparse.Namespace,
    polygon_number: int,
    row: Any,
    shapefile: Path,
    output_base: Path,
) -> dict[str, Any]:
    polygon_id = f"polygon_{polygon_number}"
    save_dir = output_base / polygon_id
    save_dir.mkdir(parents=True, exist_ok=True)

    if already_complete(save_dir) and not args.overwrite:
        print(f"\n[{polygon_id}] Already complete. Skipping.")
        return {
            "polygon_id": polygon_id,
            "status": "skipped",
            "sit_2017": row.get("SIT_2017", ""),
            "message": "already complete",
            "elapsed_seconds": 0.0,
        }

    started = time.time()
    try:
        polygon = clean_polygon(row.geometry)
        metadata = {
            "polygon_id": polygon_id,
            "source_shapefile": shapefile,
            "source_row_index": row.name,
            "properties": row_properties(row),
            "bounds_epsg4326": list(polygon.bounds),
            "network_type": args.network_type,
            "truncate_by_edge": args.truncate_by_edge,
            "img_size": args.img_size,
            "street_width": args.street_width,
            "v": args.v,
            "dt": args.dt,
            "n_steps": args.n_steps,
            "min_h1_steps": args.min_h1_steps,
            "num_landscapes": args.num_landscapes,
            "resolution": args.resolution,
        }
        write_json(save_dir / "metadata.json", metadata)
        write_json(
            save_dir / "status.json",
            {
                "polygon_id": polygon_id,
                "status": "running",
                "started_at_epoch": started,
            },
        )

        print(f"\n[{polygon_id}] Fetching {args.network_type} network from polygon...")
        nodes, edges, x_min, x_max, y_min, y_max = fetch_network_polygon(
            ox,
            polygon,
            args.network_type,
            args.truncate_by_edge,
        )

        street_binary = utils["rasterize_streets"](
            edges,
            nodes,
            args.img_size,
            x_min,
            x_max,
            y_min,
            y_max,
            line_width=args.street_width,
        )
        utils["plot_street_binary"](
            street_binary,
            save_files=True,
            save_path=str(save_dir / "street_binary.pdf"),
        )

        phi0 = utils["build_sdf"](street_binary)
        utils["plot_sdf"](phi0, save_files=True, save_path=str(save_dir / "sdf.pdf"))

        arrival_time_map, snapshots = utils["propagate"](
            phi0,
            args.img_size,
            args.dt,
            args.v,
            args.n_steps,
        )
        utils["plot_propagation_steps"](
            snapshots,
            phi0,
            street_binary,
            args.dt,
            save_files=True,
            save_path=str(save_dir / "levelset_propagation_steps.pdf"),
        )
        utils["plot_arrival_time_map"](
            arrival_time_map,
            save_files=True,
            save_path=str(save_dir / "arrival_time_map.pdf"),
        )

        pairs_h0, pairs_h1 = utils["compute_persistence"](
            arrival_time_map,
            dt=args.dt,
            min_h1_steps=args.min_h1_steps,
        )
        inf_death_cap = (
            float(arrival_time_map.max()) * 1.05
            if float(arrival_time_map.max()) > 0
            else 1.0
        )

        utils["plot_persistence_diagrams"](
            pairs_h0,
            pairs_h1,
            inf_death_cap=inf_death_cap,
            save_files=True,
            save_path=str(save_dir / "persistence_diagrams.pdf"),
        )
        utils["plot_barcodes"](
            pairs_h0,
            pairs_h1,
            inf_death_cap=inf_death_cap,
            save_files=True,
            save_path=str(save_dir / "persistence_barcodes.pdf"),
        )

        h1_array = (
            np.array(pairs_h1, dtype=float)
            if pairs_h1
            else np.empty((0, 2), dtype=float)
        )
        np.savetxt(save_dir / "h1_pairs.csv", h1_array, delimiter=",")

        if pairs_h1:
            landscapes, x_grid, flat_vec, landscape_obj = utils["compute_landscape"](
                pairs_h1,
                args.num_landscapes,
                args.resolution,
            )
            utils["plot_landscape"](
                landscapes,
                x_grid,
                save_files=True,
                save_path=str(save_dir / "persistence_landscape.pdf"),
            )
            np.savez(
                save_dir / "persistence_landscape_h1.npz",
                landscape_values=landscapes,
                x_grid=x_grid,
                flat_vector=flat_vec,
                num_landscapes=landscape_obj.num_landscapes,
                resolution=landscape_obj.resolution,
                empty=False,
            )
        else:
            print(f"  [{polygon_id}] No H1 features found. Saving zero landscape.")
            save_empty_landscape(
                save_dir,
                inf_death_cap,
                args.num_landscapes,
                args.resolution,
                utils["plot_landscape"],
            )

        elapsed = time.time() - started
        status_payload = {
            "polygon_id": polygon_id,
            "status": "complete",
            "elapsed_seconds": elapsed,
            "h0_features": len(pairs_h0),
            "h1_features": len(pairs_h1),
            "output_dir": save_dir,
        }
        write_json(save_dir / "status.json", status_payload)

        print(f"  [{polygon_id}] Done in {elapsed:.1f}s -> {save_dir}")
        return {
            "polygon_id": polygon_id,
            "status": "complete",
            "sit_2017": row.get("SIT_2017", ""),
            "message": "",
            "elapsed_seconds": round(elapsed, 3),
            "h0_features": len(pairs_h0),
            "h1_features": len(pairs_h1),
        }
    except Exception as exc:
        elapsed = time.time() - started
        message = str(exc)
        write_json(
            save_dir / "status.json",
            {
                "polygon_id": polygon_id,
                "status": "failed",
                "elapsed_seconds": elapsed,
                "error": message,
                "output_dir": save_dir,
            },
        )
        print(f"  ERROR [{polygon_id}]: {message}")
        if args.fail_fast:
            raise
        return {
            "polygon_id": polygon_id,
            "status": "failed",
            "sit_2017": row.get("SIT_2017", ""),
            "message": message,
            "elapsed_seconds": round(elapsed, 3),
            "h0_features": "",
            "h1_features": "",
        }
    finally:
        gc.collect()


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fields = [
        "polygon_id",
        "status",
        "sit_2017",
        "elapsed_seconds",
        "h0_features",
        "h1_features",
        "message",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    args.shapefile = resolve_project_path(args.shapefile)
    args.output_base = resolve_project_path(args.output_base)
    validate_args(args)

    if str(SRC_DIR) not in sys.path:
        sys.path.insert(0, str(SRC_DIR))

    gdf = load_polygons(args.shapefile)
    selected = list(selected_rows(gdf, args.start_index, args.limit))
    print(f"Loaded {len(gdf)} polygons from {args.shapefile}")
    print(f"Selected {len(selected)} polygons -> {args.output_base}")

    if args.dry_run:
        for polygon_number, row in selected:
            sit_2017 = row.get("SIT_2017", "")
            print(f"  polygon_{polygon_number}: SIT_2017={sit_2017}")
        return 0
    if not selected:
        print("No polygons selected.")
        return 0

    utils = import_pipeline_utils()
    ox = configure_osmnx(args)
    args.output_base.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, Any]] = []
    try:
        for polygon_number, row in selected:
            summary_rows.append(
                run_polygon(
                    ox=ox,
                    utils=utils,
                    args=args,
                    polygon_number=polygon_number,
                    row=row,
                    shapefile=args.shapefile,
                    output_base=args.output_base,
                )
            )
    finally:
        summary_path = args.output_base / "run_summary.csv"
        write_summary(summary_path, summary_rows)
        if summary_rows:
            print(f"\nSummary saved -> {summary_path}")

    completed = sum(1 for row in summary_rows if row["status"] == "complete")
    skipped = sum(1 for row in summary_rows if row["status"] == "skipped")
    failed = sum(1 for row in summary_rows if row["status"] == "failed")
    print(f"Complete: {completed}, skipped: {skipped}, failed: {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        sys.stdout.flush()
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
