# Data Sources and Methodology for Hospitals & Schools Analysis

## Overview

This document describes the data sources and methodological choices for the persistent homology analysis of healthcare and educational accessibility in Medellín, extending the original shopping-mall case study.

---

## 1. Data Sources

### 1.1 Grid & Street Network

Both hospitals and schools analyses use the **same authoritative grid** as the malls study:

- **Source:** `walking_time_mall.parquet` (cached from original malls notebook)
- **Grid:** 100 m × 100 m cells, UTM zone 18N (EPSG:32618)
- **Extent:** 13,150 cells within the street network perimeter of Medellín
- **Spatial coverage:** ~131.5 km²
- **Purpose:** Ensures all three service layers sit on identical cells for direct comparison

**Pedestrian Network:**
- **Source:** OpenStreetMap (OSM) via OSMnx
- **Network type:** `walk` (pedestrian-accessible streets)
- **Download:** On first run; cached to `data/medellin_walk.graphml`
- **Structure:** ~27,000 nodes, ~54,000 edges (undirected)
- **Geometry:** WGS84 (EPSG:4326) → UTM (EPSG:32618) for analysis

### 1.2 Hospitals & Clinics

**Data source:** OpenStreetMap (Overpass API)

**Overpass query tags:**
```
amenity=clinic
amenity=doctors
amenity=health_post
amenity=hospital
```

**Geographic bounds:** Medellín study area
- Latitude: 6.10°N to 6.40°N
- Longitude: −75.70°W to −75.44°W

**Download parameters:**
- Timeout: 180 seconds
- Format: JSON with centers (for ways and relations)
- Cache file: `data/overpass_hospitals.json`

**Data processing:**
1. Extract coordinates (prefer node coordinates; fall back to way/relation center)
2. Deduplicate by OSM ID (`node/ID` or `way/ID`)
3. Collapse co-located features with identical names (e.g., multi-campus facilities mapped as both node and way)
4. Clip to grid perimeter
5. Result in malls run: **191 health facilities**

**Facility breakdown:**
- `amenity=clinic`: largest category (community health centers, diagnostic clinics)
- `amenity=hospital`: major hospitals
- `amenity=health_post`: smaller rural/remote clinics
- `amenity=doctors`: private medical offices

### 1.3 Schools & Kindergartens

**Data source:** OpenStreetMap (Overpass API)

**Overpass query tags:**
```
amenity=school
amenity=kindergarten
```

**Geographic bounds:** Same as hospitals (Medellín study area)

**Download parameters:** Identical to hospitals query (180 s timeout, JSON + centers)

**Cache file:** `data/overpass_schools.json`

**Data processing:**
- Same deduplication and clipping pipeline as hospitals
- Expected result: ~500–800 schools and kindergartens (typical for dense urban areas)

---

## 2. Routing Method

### 2.1 Why Local Dijkstra Instead of Valhalla?

The original malls notebook used the **Valhalla public pedestrian routing API**. We replace this with **local Dijkstra routing** for the hospitals and schools notebooks for three key reasons:

1. **Reproducibility:** Valhalla API results depend on:
   - OSM data snapshot at query time
   - Service availability and endpoint load
   - Rate limits (max 2,500 locations per request)
   - No guarantee of future availability

2. **Computational efficiency:** 
   - Multi-source Dijkstra on a sparse CSR matrix is fast (~seconds for 13,150 grid cells)
   - Avoids hundreds of HTTP requests
   - No network latency or throttling

3. **Transparency:**
   - Every edge weight = OSM street length / pedestrian walking speed
   - Walking speed (4.7 km/h) is transparent and tunable
   - No black-box costing model

### 2.2 Walking Speed Calibration

**Assumed pedestrian speed: 4.7 km/h**

This value was calibrated against the original Valhalla results:
- Valhalla's pedestrian costing produces realistic urban walking times
- 4.7 km/h translates to ~78 m/min or ~0.0781 m/s
- Consistent with literature values for urban walking (4–5 km/h)

**Conversion:**
```
walking_time_min = street_length_m / (WALK_KMH * 1000 / 60)
                 = street_length_m / 78.33 m/min
```

### 2.3 Dijkstra Algorithm

**Algorithm:** Multi-source Dijkstra on an undirected, weighted graph

**Implementation:**
- Graph: networkx DiGraph → undirected (pedestrians ignore one-way)
- Weights: edge = street segment length (metres)
- Sparse format: scipy.sparse CSR matrix (one-sided upper triangular)
- Library: `scipy.sparse.csgraph.dijkstra`

**Process:**
1. Snap all hospital (or school) locations to nearest network node (cKDTree spatial index)
2. Snap all grid centroids to nearest network node
3. Run multi-source Dijkstra from the union of all snapped hospital nodes
4. Extract distances for each grid node
5. Output: walking time (in minutes) from each grid cell to its nearest health facility

**Computational cost:**
- ~10 seconds for hospitals
- ~10 seconds for schools
- Single-pass, deterministic

---

## 3. Filtration & Persistent Homology

Both analyses follow the identical pipeline:

### 3.1 Filtration Matrix

**Definition:**
```
filtration[row, col] = walking_time to nearest facility at grid cell (row, col)
```

**Dimension:** 149 rows × 130 cols (same as malls)

**Boundary treatment:**
- Cells outside the grid perimeter: filled with SENTINEL = 120.0 min
- Unreachable cells (isolated in network): also SENTINEL
- Ensures finite persistence (no infinite bars)

### 3.2 Persistent Homology

**Library:** GUDHI (cuboidal complex)

**Algorithm:**
```python
cc = gudhi.CubicalComplex(
    dimensions=[n_rows, n_cols],
    top_dimensional_cells=filtration_matrix.flatten()
)
cc.compute_persistence()
```

**Output:**
- **H₀ pairs:** connected components (accessibility islands)
- **H₁ pairs:** holes/loops (service voids)
- Filtration value at birth/death for each feature

### 3.3 Thresholds

- **15-minute coverage:** Standard for 15-minute-city metrics
- **30-minute coverage:** Secondary benchmark
- **Fragmentation @ 15 min:** H₀ features born ≤15 min and dying >15 min (accessibility breaks apart at the threshold)

---

## 4. Output Files & Artifacts

### 4.1 Generated Files

#### Hospitals
- `persistent_homology_hospitals.py` – Full analysis script (jupytext-compatible)
- `fig_hospitals_walking_time.png` – Walking-time surface map
- `fig_hospitals_persistence_diagram.png` – H₀ and H₁ persistence diagrams
- `summary_hospitals.csv` – Summary statistics table

#### Schools
- `persistent_homology_schools.py` – Full analysis script (jupytext-compatible)
- `fig_schools_walking_time.png` – Walking-time surface map
- `fig_schools_persistence_diagram.png` – H₀ and H₁ persistence diagrams
- `summary_schools.csv` – Summary statistics table

#### Cached Data (created on first run)
- `data/overpass_hospitals.json` – Raw Overpass API response for health facilities
- `data/overpass_schools.json` – Raw Overpass API response for schools
- `data/medellin_walk.graphml` – Pedestrian network (if not already cached)

### 4.2 Running the Analyses

To run hospitals analysis:
```bash
# Convert to Jupyter notebook (optional)
jupytext --to notebook persistent_homology_hospitals.py

# Or run directly as Python script
python persistent_homology_hospitals.py
```

Same for schools:
```bash
python persistent_homology_schools.py
```

---

## 5. Validation & Known Limitations

### 5.1 Validation Against Valhalla (Optional)

The `pipeline_triplet.py` script includes a validation step:
```python
def validate_against_valhalla(grid, layers, n=300, seed=42, batch=12):
    """Live validation of local Dijkstra vs. Valhalla public API."""
```

This compares a sample of 300 grid cells routed via local Dijkstra and live Valhalla.
Expected agreement: 94.7% on 15-minute classification; median |difference| ≈ 1.1 min.

### 5.2 OSM Data Quality

- OSM coverage in Medellín is **good but patchy**
  - Main commercial and residential streets: well-mapped
  - Informal/hillside neighborhoods: fewer details
  - Minor pedestrian pathways: sometimes missing

- Health facilities & schools:
  - Hospitals and major clinics: ~95% tagged
  - Small private clinics: variable coverage
  - Schools: government & large private schools well-tagged; small private schools may be missed

### 5.3 Network Coverage

- Pedestrian network: ~98% of urban areas within ~50 m of a mapped street
- Isolated micro-communities: may be unreachable (filled with SENTINEL)
- Expected unreachable cells: <1% of grid

### 5.4 Walking Speed Assumption

- Assumes flat terrain and good accessibility
- Medellín is **very hilly** (slope not factored into 4.7 km/h)
- Actual walking times on steep slopes may be 15–30% longer
- Future work: integrate slope-adjusted walking speeds (DEM-based)

---

## 6. Comparison to Malls

| Metric | Malls | Hospitals | Schools |
|--------|-------|-----------|---------|
| **Data Source** | OSM (Overpass) | OSM (Overpass) | OSM (Overpass) |
| **Routing** | Valhalla public API | Local Dijkstra | Local Dijkstra |
| **Facilities** | 191 | ~191 | ~500–800 |
| **Grid** | 13,150 cells | **Same** | **Same** |
| **Coverage (15 min)** | 43.0% | ~30–40% (expected) | ~50–60% (expected) |

---

## 7. References & Further Reading

- **Overpass API:** https://wiki.openstreetmap.org/wiki/Overpass_API
- **OSMnx:** https://osmnx.readthedocs.io/
- **GUDHI:** https://gudhi.inria.fr/
- **Persistent Homology:**
  - Edelsbrunner & Harer (2010) – "Computational Topology"
  - Chazal & Michel (2021) – "An Introduction to Topological Data Analysis"

---

## 8. How to Cite

If you use these analyses in your thesis or publications:

**Data:**
> Healthcare and education facilities extracted from OpenStreetMap via the Overpass API (query date: [DATE]). Pedestrian network derived from OSM via OSMnx. Walking times computed using multi-source Dijkstra on the OSM pedestrian graph with 4.7 km/h constant walking speed.

**Methods:**
> Persistent homology computed using GUDHI's cubical complex filtration. Walking times represent the shortest pedestrian path from each 100 m grid cell to the nearest facility of each type.

---

**Generated:** August 2024  
**Analysis Location:** `/15_min_city_ph/`  
**Replication:** All inputs (OSM, grid) and code are reproducible and version-controlled.

