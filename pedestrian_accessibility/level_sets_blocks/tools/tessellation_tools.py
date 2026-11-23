"""Topology repair and Momepy tessellation helpers for Colombian city blocks."""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely


def _polygonal_only(geometry):
    """Return all polygonal parts and discard zero-area line/point remnants."""
    stack = [geometry]
    polygons = []
    while stack:
        current = stack.pop()
        type_id = shapely.get_type_id(current)
        if type_id == 3:  # Polygon
            polygons.append(current)
        elif type_id in (6, 7):  # MultiPolygon or GeometryCollection
            stack.extend(shapely.get_parts(current))
    if not polygons:
        return shapely.GeometryCollection()
    return shapely.union_all(polygons)


def _inset(geometry, distance):
    return shapely.buffer(
        geometry,
        -distance,
        cap_style="flat",
        join_style="mitre",
    )


def _is_inset_safe(geometry, distance):
    inset = _inset(geometry, distance)
    return (not shapely.is_empty(inset)) and shapely.get_type_id(inset) == 3


def _robust_main_component(geometry, distance):
    """Open a narrow neck and retain the largest component surviving the inset."""
    inset = _inset(geometry, distance)
    parts = list(shapely.get_parts(inset))
    if not parts:
        return shapely.GeometryCollection()
    main_inset = max(parts, key=shapely.area)
    opened = shapely.buffer(
        main_inset,
        distance,
        cap_style="flat",
        join_style="mitre",
    )
    return _polygonal_only(shapely.intersection(geometry, opened))


def _positive_overlap_pairs(gdf, tolerance_m2):
    positions = gdf.sindex.query(gdf.geometry, predicate="intersects")
    left, right = positions
    unique = left < right
    left, right = left[unique], right[unique]
    geometry = np.asarray(gdf.geometry.array, dtype=object)
    intersections = shapely.intersection(geometry[left], geometry[right])
    areas = shapely.area(intersections)
    keep = areas > tolerance_m2
    return left[keep], right[keep], intersections[keep], areas[keep]


def repair_block_generators(
    medellin_blocks: gpd.GeoDataFrame,
    medellin_polygon: gpd.GeoDataFrame,
    *,
    target_shrink_m: float = 0.2,
    overlap_tolerance_m2: float = 0.01,
):
    """Create non-overlapping, inset-safe, single-polygon block generators.

    The smaller block owns a shared area unless that subtraction disconnects the
    larger block. Ownership is then reversed. If neither full assignment can stay
    connected after the target inset, the direction losing the least detached area
    is used and that detached fragment is released to the later tessellation.
    """
    city_limit = shapely.set_precision(
        shapely.make_valid(medellin_polygon.geometry.union_all()),
        0.01,
    )
    generators = (
        medellin_blocks.drop_duplicates(subset="MANZ_CCNCT")
        .reset_index(drop=True)
        .copy()
    )
    generators["MANZ_CCNCT"] = generators["MANZ_CCNCT"].astype("string")

    # Keep whole source polygons. Clipping at a jagged municipal edge is what made
    # 27 otherwise valid blocks multipart in the earlier attempt.
    source_geometry = np.asarray(generators.geometry.array, dtype=object)
    city_area = shapely.area(shapely.intersection(source_geometry, city_limit))
    generators = generators[city_area > overlap_tolerance_m2].reset_index(drop=True)
    current = np.asarray(generators.geometry.array, dtype=object).copy()

    # Some DANE blocks contain a tiny, detached polygon fragment. Momepy needs
    # one generator per block, while exploding the geometry would duplicate the
    # block's demographic attributes. Keep the largest polygonal component and
    # report the discarded area instead.
    source_multipart_repairs = 0
    source_multipart_released_m2 = 0.0
    non_polygon = np.flatnonzero(shapely.get_type_id(current) != 3)
    for position in non_polygon:
        polygonal = _polygonal_only(shapely.make_valid(current[position]))
        parts = [
            part
            for part in shapely.get_parts(polygonal)
            if shapely.get_type_id(part) == 3 and shapely.area(part) > 0
        ]
        if not parts:
            block_id = generators.iloc[position]["MANZ_CCNCT"]
            raise RuntimeError(
                f"Block {block_id} has no positive-area polygon component."
            )

        main_part = max(parts, key=shapely.area)
        source_multipart_released_m2 += max(
            0.0,
            float(shapely.area(polygonal) - shapely.area(main_part)),
        )
        current[position] = main_part
        source_multipart_repairs += 1

    generators["geometry"] = current

    # Repair narrow source necks that the selected Momepy inset would split.
    initial_insets = _inset(current, target_shrink_m)
    unsafe_initial = shapely.is_empty(initial_insets) | (
        shapely.get_type_id(initial_insets) != 3
    )
    initial_released_m2 = 0.0
    for position in np.flatnonzero(unsafe_initial):
        repaired = _robust_main_component(current[position], target_shrink_m)
        if shapely.is_empty(repaired):
            block_id = generators.iloc[position]["MANZ_CCNCT"]
            raise RuntimeError(
                f"Block {block_id} collapses at a {target_shrink_m} m inset."
            )
        initial_released_m2 += float(
            shapely.area(current[position]) - shapely.area(repaired)
        )
        current[position] = repaired
    generators["geometry"] = current

    source_area = shapely.area(current)
    block_ids = generators["MANZ_CCNCT"].astype(str).to_numpy()
    left, right, _, overlap_areas = _positive_overlap_pairs(
        generators, overlap_tolerance_m2
    )

    reversals = 0
    crossing_records = []
    for pair_number in np.argsort(-overlap_areas):
        left_position = int(left[pair_number])
        right_position = int(right[pair_number])
        overlap = _polygonal_only(
            shapely.intersection(
                current[left_position],
                current[right_position],
            )
        )
        overlap_area = float(shapely.area(overlap))
        if overlap_area <= overlap_tolerance_m2:
            continue

        left_priority = (source_area[left_position], block_ids[left_position])
        right_priority = (source_area[right_position], block_ids[right_position])
        if left_priority <= right_priority:
            preferred_winner, preferred_loser = left_position, right_position
        else:
            preferred_winner, preferred_loser = right_position, left_position

        preferred_result = _polygonal_only(
            shapely.difference(current[preferred_loser], overlap)
        )
        if _is_inset_safe(preferred_result, target_shrink_m):
            current[preferred_loser] = preferred_result
            continue

        alternate_result = _polygonal_only(
            shapely.difference(current[preferred_winner], overlap)
        )
        if _is_inset_safe(alternate_result, target_shrink_m):
            current[preferred_winner] = alternate_result
            reversals += 1
            continue

        # Crossing corridors cannot both remain connected in a planar partition.
        preferred_main = _robust_main_component(
            preferred_result, target_shrink_m
        )
        alternate_main = _robust_main_component(
            alternate_result, target_shrink_m
        )
        if shapely.is_empty(preferred_main) or shapely.is_empty(alternate_main):
            raise RuntimeError("An overlap fallback collapsed a complete block.")

        preferred_loss = float(
            shapely.area(preferred_result) - shapely.area(preferred_main)
        )
        alternate_loss = float(
            shapely.area(alternate_result) - shapely.area(alternate_main)
        )
        if preferred_loss <= alternate_loss:
            current[preferred_loser] = preferred_main
            rule = "smaller_wins"
            released_m2 = preferred_loss
        else:
            current[preferred_winner] = alternate_main
            rule = "larger_wins"
            released_m2 = alternate_loss
            reversals += 1

        crossing_records.append(
            {
                "left_id": block_ids[left_position],
                "right_id": block_ids[right_position],
                "overlap_m2": overlap_area,
                "rule": rule,
                "released_m2": released_m2,
            }
        )

    generators["geometry"] = current
    residual = _positive_overlap_pairs(generators, overlap_tolerance_m2)
    final_insets = _inset(current, target_shrink_m)
    inset_empty = int(shapely.is_empty(final_insets).sum())
    inset_multipart = int((shapely.get_type_id(final_insets) == 6).sum())
    coverage_overlap_m2 = max(
        0.0,
        float(generators.geometry.area.sum() - generators.geometry.union_all().area),
    )

    if len(residual[0]) or inset_empty or inset_multipart:
        raise RuntimeError(
            "Generator repair failed its overlap/inset validation: "
            f"pairs={len(residual[0])}, empty={inset_empty}, "
            f"multipart={inset_multipart}."
        )

    diagnostics = {
        "input_blocks": len(generators),
        "source_multipart_repairs": source_multipart_repairs,
        "source_multipart_released_m2": source_multipart_released_m2,
        "original_overlap_pairs": len(left),
        "ownership_reversals": reversals,
        "crossing_fallbacks": len(crossing_records),
        "initial_topology_repairs": int(unsafe_initial.sum()),
        "initial_released_m2": initial_released_m2,
        "crossing_released_m2": sum(
            record["released_m2"] for record in crossing_records
        ),
        "remaining_overlap_pairs": len(residual[0]),
        "coverage_overlap_m2": coverage_overlap_m2,
        "inset_empty": inset_empty,
        "inset_multipart": inset_multipart,
    }
    return city_limit, generators, diagnostics, pd.DataFrame(crossing_records)


def build_momepy_tessellation(
    block_generators: gpd.GeoDataFrame,
    city_limit,
    *,
    shrink_m: float = 0.2,
    segment_m: float = 2.0,
):
    """Run the best completed Momepy benchmark and return explicit diagnostics."""
    import momepy

    indexed = (
        block_generators.set_index("MANZ_CCNCT", drop=False)
        .sort_index()
        .copy()
    )
    cells = momepy.morphological_tessellation(
        indexed.geometry,
        clip=city_limit,
        shrink=shrink_m,
        segment=segment_m,
    )
    cells.index.name = "MANZ_CCNCT"

    attributes = indexed.drop(columns=["geometry", "MANZ_CCNCT"], errors="ignore")
    tessellation = cells.join(attributes, how="left", validate="one_to_one")
    tessellation["MANZ_CCNCT"] = tessellation.index.astype(str)
    tessellation["cell_area_km2"] = tessellation.geometry.area / 1_000_000

    union = tessellation.geometry.union_all()
    missing = indexed.index.difference(tessellation.index)
    clipped_sources = shapely.intersection(indexed.geometry.array, city_limit)
    owned_cells = tessellation.geometry.reindex(indexed.index).array
    source_outside_own = float(
        np.nansum(shapely.area(shapely.difference(clipped_sources, owned_cells)))
    )
    diagnostics = {
        "shrink_m": shrink_m,
        "segment_m": segment_m,
        "input_generators": len(indexed),
        "cells": len(tessellation),
        "missing_cells": len(missing),
        "polygon_cells": int(tessellation.geom_type.eq("Polygon").sum()),
        "multipolygon_cells": int(tessellation.geom_type.eq("MultiPolygon").sum()),
        "uncovered_m2": float(city_limit.difference(union).area),
        "outside_city_m2": float(union.difference(city_limit).area),
        "overlap_m2": max(
            0.0,
            float(tessellation.geometry.area.sum() - union.area),
        ),
        "source_outside_own_cell_m2": source_outside_own,
    }
    return indexed, tessellation, diagnostics
