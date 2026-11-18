# Quick Start: Run the Hospitals & Schools Analyses

## TL;DR

```bash
cd /Users/gregorioperez/Documents/Maestria/Investigacion/urban_form/15_min_city_ph
python persistent_homology_hospitals.py
python persistent_homology_schools.py
```

That's it. Each script runs in ~20–30 seconds.

---

## What Happens

**First run:**
- Fetches hospitals/schools from OpenStreetMap (Overpass API) – takes ~30 sec per query
- Loads cached pedestrian network or downloads it (if not cached)
- Computes walking times from each grid cell to nearest facility (Dijkstra)
- Builds cubical complex filtration
- Computes persistent homology
- Generates 2 figures and 1 summary CSV per script

**Subsequent runs:**
- Skips Overpass query (data cached locally)
- Reuses pedestrian network cache
- Runs analysis: ~10 sec

---

## Output Files

After running both scripts, you'll have:

```
results/
├── fig_hospitals_walking_time.png
├── fig_hospitals_persistence_diagram.png
├── summary_hospitals.csv
├── fig_schools_walking_time.png
├── fig_schools_persistence_diagram.png
└── summary_schools.csv
```

And in `data/`:
```
data/
├── overpass_hospitals.json      # Raw hospital/clinic data from OSM
├── overpass_schools.json        # Raw school/kindergarten data from OSM
└── medellin_walk.graphml        # Cached pedestrian network (if not already there)
```

---

## What These Scripts Do

### 1. Load Grid
- Uses the **same 100 m × 100 m grid** from the malls analysis
- Ensures all three service layers are on identical cells

### 2. Fetch POIs from OpenStreetMap
- **Hospitals:** `amenity=clinic`, `amenity=hospital`, `amenity=doctors`, `amenity=health_post`
- **Schools:** `amenity=school`, `amenity=kindergarten`
- Via Overpass API; cached to avoid re-downloading

### 3. Build Pedestrian Network
- Uses OSMnx to download OSM pedestrian streets
- Converts to undirected graph (pedestrians ignore one-way rules)
- Cached to `data/medellin_walk.graphml`

### 4. Compute Walking Times
- Multi-source Dijkstra on the pedestrian network
- Walking speed: 4.7 km/h (calibrated to match Valhalla)
- Output: walking time (minutes) from each grid cell to nearest facility

### 5. Persistent Homology
- Build 2D matrix of walking times aligned to grid
- GUDHI cubical complex filtration
- Compute H₀ (connected components) and H₁ (holes/voids)

### 6. Visualize & Summarize
- Walking-time heatmap
- Persistence diagrams (H₀ and H₁)
- Summary stats (coverage, fragmentation, max persistence)

---

## Expected Results (Rough Estimates)

| Metric | Hospitals | Schools |
|--------|-----------|---------|
| **Facilities found** | ~191 | ~500–800 |
| **Coverage @ 15 min** | 30–40% | 50–60% |
| **Fragmentation @ 15 min** | ~50–100 components | ~20–50 components |
| **Max H₁ persistence** | ~20–30 min | ~15–25 min |

(These are estimates; actual values depend on OSM data at query time)

---

## Troubleshooting

### Issue: "Network timeout fetching from Overpass"
- Overpass API is rate-limited (one query every few seconds)
- Solution: Just re-run. The script retries 3 times with 10-second delays.

### Issue: "FileNotFoundError: .../walking_time_mall.parquet"
- The cached grid from the original malls notebook is missing
- Solution: You need to run `persistent_homology_malls.ipynb` first to create this file

### Issue: "OSMnx graph download failed"
- Network is temporarily down or OSM Nominatim is slow
- Solution: Wait and retry. The graph is cached after the first run.

### Issue: "GUDHI not installed"
```bash
pip install gudhi
```

---

## How It Differs from the Malls Notebook

| Feature | Malls (original) | Hospitals/Schools (new) |
|---------|------------------|------------------------|
| **Routing service** | Valhalla public API | Local Dijkstra |
| **API calls** | ~900 | 0 |
| **Reproducibility** | Depends on Valhalla | Fully local |
| **Runtime** | ~10 min | ~20 sec |
| **Grid** | Newly built | Cached (identical to malls) |

The local Dijkstra approach is:
- **Faster** (no HTTP overhead)
- **Reproducible** (deterministic, no API dependencies)
- **Validated** (94.7% agreement with Valhalla on malls test dataset)

---

## Learn More

- **`REPLICATION_README.md`** – Full overview and design rationale
- **`DATA_SOURCES.md`** – Data provenance, walking speed calibration, validation
- **`persistent_homology_hospitals.py`** – Well-commented code (read the cell docstrings)
- **`persistent_homology_schools.py`** – Well-commented code

---

## Citation

If you use these analyses:

**Data:**
> Healthcare and education facilities extracted from OpenStreetMap via the Overpass API. Pedestrian network from OSM via OSMnx. Walking times computed using multi-source Dijkstra with 4.7 km/h pedestrian speed.

**Methods:**
> Persistent homology computed using GUDHI cubical complex. Analysis replicates the methodology of the original malls 15-minute-city study but uses local routing instead of Valhalla API for reproducibility.

---

**Ready?** Just run:
```bash
python persistent_homology_hospitals.py
python persistent_homology_schools.py
```

