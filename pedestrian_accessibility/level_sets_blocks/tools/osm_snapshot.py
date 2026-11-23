"""Offline OpenStreetMap road snapshot for filling gaps in the base networks.

The base networks (Boeing's OSM street networks, one per urban area) are
bounded by urban-area polygons, not by DANE urban perimeters, so some DANE
blocks are far from any base road. This module turns ONE dated national OSM
extract (e.g. Geofabrik ``colombia-YYMMDD.osm.pbf``) into a GeoPackage of
unsimplified road segments for every DANE perimeter that can enter a metro
analysis. The analysis then reads that file; it never calls a web service.

Segments use OSM node IDs as endpoints, the same IDs as the base networks, so
supplemental roads connect to the base graph where they meet it.

Build once (see build_osm_snapshot.py), commit/keep the output with its JSON
manifest, and cite the OSM timestamp from the manifest.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import LineString

# OSMnx 2.0.6 `_get_network_filter` strings, copied verbatim so the snapshot
# selects the same ways OSMnx would (default access filter included). The base
# networks are drivable networks: on a Medellín test, 16.6% of OSM nodes on
# drivable ways appear in the base graph, versus 0.9% of nodes found only on
# footway/steps/path ways. Use "drive" so supplemented areas are not denser
# than the base network; "walk" is available if the whole study moves to it.
_DEFAULT_ACCESS = '["access"!~"private"]'
OSM_WAY_FILTERS = {
    "drive": (
        f'["highway"]["area"!~"yes"]{_DEFAULT_ACCESS}'
        '["highway"!~"abandoned|bridleway|bus_guideway|construction|corridor|'
        "cycleway|elevator|escalator|footway|no|path|pedestrian|planned|platform|"
        'proposed|raceway|razed|rest_area|service|services|steps|track"]'
        '["motor_vehicle"!~"no"]["motorcar"!~"no"]'
        '["service"!~"alley|driveway|emergency_access|parking|parking_aisle|private"]'
    ),
    "walk": (
        f'["highway"]["area"!~"yes"]{_DEFAULT_ACCESS}'
        '["highway"!~"abandoned|bus_guideway|construction|cycleway|motor|no|planned|'
        'platform|proposed|raceway|razed|rest_area|services"]'
        '["foot"!~"no"]["service"!~"private"]'
        '["sidewalk"!~"separate"]["sidewalk:both"!~"separate"]'
        '["sidewalk:left"!~"separate"]["sidewalk:right"!~"separate"]'
    ),
}

_CLAUSE = re.compile(r'\["([^"]+)"(?:(!~|~|!=|=)"([^"]*)")?\]')


def parse_way_filter(way_filter: str) -> list[tuple[str, str | None, object]]:
    """Parse an Overpass tag filter into (key, operator, value/regex) clauses."""
    clauses = []
    position = 0
    for match in _CLAUSE.finditer(way_filter):
        if match.start() != position:
            raise ValueError(f"Unparsed filter text near: {way_filter[position:]!r}")
        key, operator, value = match.groups()
        if operator in {"~", "!~"}:
            value = re.compile(value)
        clauses.append((key, operator, value))
        position = match.end()
    if position != len(way_filter):
        raise ValueError(f"Unparsed filter text near: {way_filter[position:]!r}")
    return clauses


def way_passes(tags, clauses) -> bool:
    """Evaluate Overpass semantics: `!~`/`!=` also pass when the key is absent."""
    for key, operator, value in clauses:
        tag = tags.get(key)
        if operator is None:
            if tag is None:
                return False
        elif operator == "~":
            if tag is None or not value.search(tag):
                return False
        elif operator == "!~":
            if tag is not None and value.search(tag):
                return False
        elif operator == "=":
            if tag != value:
                return False
        elif operator == "!=":
            if tag == value:
                return False
    return True


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_segment_snapshot(
    pbf_path: Path,
    regions_wgs84: gpd.GeoSeries,
    out_path: Path,
    *,
    network_type: str = "drive",
    region_buffer_m: float = 500.0,
) -> dict:
    """Write unsimplified OSM road segments intersecting the regions.

    regions_wgs84: DANE perimeters (EPSG:4326) that may need supplementing.
    Each output row is one OSM node-to-node segment (from, to = OSM node IDs).
    A JSON manifest with the OSM timestamp and filter is written next to it.
    """
    import osmium

    clauses = parse_way_filter(OSM_WAY_FILTERS[network_type])

    metric = regions_wgs84.estimate_utm_crs()
    region_union = (
        regions_wgs84.to_crs(metric).buffer(region_buffer_m).to_crs("EPSG:4326").union_all()
    )
    shapely.prepare(region_union)
    min_x, min_y, max_x, max_y = region_union.bounds

    from_ids, to_ids, way_ids, highways, segments = [], [], [], [], []
    processor = (
        osmium.FileProcessor(str(pbf_path))
        .with_locations()
        .with_filter(osmium.filter.KeyFilter("highway"))
    )
    for obj in processor:
        if not obj.is_way() or not way_passes(obj.tags, clauses):
            continue
        try:
            refs = [(node.ref, node.lon, node.lat) for node in obj.nodes]
        except osmium.InvalidLocationError:
            refs = [(n.ref, n.lon, n.lat) for n in obj.nodes if n.location.valid()]
        if len(refs) < 2:
            continue
        xs = [r[1] for r in refs]
        ys = [r[2] for r in refs]
        if max(xs) < min_x or min(xs) > max_x or max(ys) < min_y or min(ys) > max_y:
            continue
        if not region_union.intersects(LineString(list(zip(xs, ys)))):
            continue
        highway = obj.tags.get("highway")
        for (u, ux, uy), (v, vx, vy) in zip(refs[:-1], refs[1:]):
            if u == v:
                continue
            from_ids.append(u)
            to_ids.append(v)
            way_ids.append(obj.id)
            highways.append(highway)
            segments.append(((ux, uy), (vx, vy)))

    geometry = shapely.linestrings(np.asarray(segments, dtype=float))
    keep = shapely.intersects(region_union, geometry)
    snapshot = gpd.GeoDataFrame(
        {
            "from": np.asarray(from_ids, dtype=np.int64)[keep].astype(str),
            "to": np.asarray(to_ids, dtype=np.int64)[keep].astype(str),
            "osm_way_id": np.asarray(way_ids, dtype=np.int64)[keep],
            "highway": np.asarray(highways, dtype=object)[keep],
        },
        geometry=geometry[keep],
        crs="EPSG:4326",
    )
    # One undirected segment per node pair (two-way ways are not duplicated).
    pair = pd.Series(
        [f"{min(a, b)}|{max(a, b)}" for a, b in zip(snapshot["from"], snapshot["to"])],
        index=snapshot.index,
    )
    snapshot = snapshot.loc[~pair.duplicated()].reset_index(drop=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    partial = out_path.with_name(out_path.stem + ".partial.gpkg")
    partial.unlink(missing_ok=True)
    snapshot.to_file(partial, layer="segments", driver="GPKG")
    partial.replace(out_path)

    reader = osmium.io.Reader(str(pbf_path), osmium.osm.osm_entity_bits.NOTHING)
    osm_timestamp = reader.header().get("osmosis_replication_timestamp")
    reader.close()
    manifest = {
        "source_file": pbf_path.name,
        "source_md5": _md5(pbf_path),
        "osm_data_timestamp": osm_timestamp,
        "network_type": network_type,
        "way_filter": OSM_WAY_FILTERS[network_type],
        "region_buffer_m": region_buffer_m,
        "segments": int(len(snapshot)),
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    manifest_path(out_path).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def manifest_path(snapshot_path: Path) -> Path:
    return snapshot_path.with_suffix(".json")


def read_manifest(snapshot_path: Path) -> dict:
    return json.loads(manifest_path(snapshot_path).read_text(encoding="utf-8"))


def read_segments(snapshot_path: Path, mask_wgs84) -> gpd.GeoDataFrame:
    """Read snapshot segments intersecting a WGS84 mask (uses the GPKG index)."""
    return gpd.read_file(snapshot_path, layer="segments", mask=mask_wgs84)
