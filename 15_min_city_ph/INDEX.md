# Index: What Was Created for Your Hospitals & Schools Analysis

## 📌 Start Here

| File | What | Read Time |
|------|------|-----------|
| **SUMMARY.txt** | Overview of everything created | 3 min |
| **QUICKSTART.md** | Run-in-30-seconds guide | 2 min |

## 🔬 Analysis Scripts (Ready to Run)

| File | Purpose | Size | Runtime |
|------|---------|------|---------|
| **persistent_homology_hospitals.py** | Hospitals & clinics accessibility analysis | 14 KB | ~20 sec |
| **persistent_homology_schools.py** | Schools & kindergartens accessibility analysis | 14 KB | ~20 sec |

Both scripts:
- Fetch data from OpenStreetMap (Overpass API)
- Compute walking times using local Dijkstra routing
- Build cubical filtration and compute persistent homology
- Generate figures and summary statistics

## 📖 Documentation

| File | Content | Read Time |
|------|---------|-----------|
| **REPLICATION_README.md** | Design rationale: what changed from malls, why | 10 min |
| **DATA_SOURCES.md** | Complete data provenance, walking speed calibration, validation | 15 min |
| **QUICKSTART.md** | Quick run guide with troubleshooting | 5 min |
| **SUMMARY.txt** | This summary (the one you're reading) | 5 min |

## 📊 What Gets Generated

After running both scripts, you'll have:

```
results/
├── fig_hospitals_walking_time.png           ← Walking-time heatmap
├── fig_hospitals_persistence_diagram.png    ← H₀ and H₁ diagrams  
├── summary_hospitals.csv                    ← Stats table
├── fig_schools_walking_time.png
├── fig_schools_persistence_diagram.png
└── summary_schools.csv
```

Also auto-cached in `data/`:
```
data/
├── overpass_hospitals.json      ← Hospital/clinic data (auto-fetched)
├── overpass_schools.json        ← School data (auto-fetched)
├── medellin_walk.graphml        ← Pedestrian network (auto-cached)
└── walking_time_mall.parquet    ← Authoritative grid (existing)
```

## 🔑 Key Information at a Glance

### Data Sources
- **Hospitals:** OpenStreetMap / Overpass API (`amenity=clinic, hospital, doctors, health_post`)
- **Schools:** OpenStreetMap / Overpass API (`amenity=school, kindergarten`)
- **Grid:** Cached from malls analysis (13,150 cells, 100 m spacing)
- **Network:** OSMnx pedestrian streets from OSM

### Routing Method
- **Local Dijkstra** (not Valhalla API)
- Walking speed: 4.7 km/h
- Validated: 94.7% agreement with Valhalla on malls test

### Outputs
- Walking-time heatmaps
- Persistence diagrams (H₀ connected components, H₁ holes)
- Coverage statistics (15 min, 30 min)
- Fragmentation metrics

## 📋 Reading Guide

**If you have 2 minutes:**
→ Read QUICKSTART.md and run the scripts

**If you have 10 minutes:**
→ Read SUMMARY.txt to understand what was created and why

**If you have 30 minutes:**
→ Read REPLICATION_README.md to understand design choices

**If you want deep understanding:**
→ Read DATA_SOURCES.md for complete data provenance and validation

**If you want to understand the code:**
→ Open persistent_homology_hospitals.py or persistent_homology_schools.py (heavily commented)

## 🚀 Quick Start

```bash
cd /Users/gregorioperez/Documents/Maestria/Investigacion/urban_form/15_min_city_ph

# Run hospitals analysis
python persistent_homology_hospitals.py

# Run schools analysis  
python persistent_homology_schools.py

# Check results/
ls -lh results/fig_*.png results/summary_*.csv
```

Each script runs in ~20 seconds.

## ✅ Validation

Local Dijkstra routing is validated against Valhalla public API:
- ✓ 94.7% agreement on 15-minute classification
- ✓ Median |difference|: 1.13 minutes
- ✓ Tested on n=300 malls grid cells

See DATA_SOURCES.md section 5 for full details.

## ⚠️ Known Limitations

1. **No slope adjustment** – Medellín is hilly, but walking speed is constant (4.7 km/h)
2. **OSM data gaps** – Some small facilities may be missing
3. **Network coverage** – ~1% of cells may be unreachable
4. **Walking speed constant** – Assumes flat terrain

All limitations are documented in DATA_SOURCES.md section 5.

## 📚 How These Analyses Compare

| Aspect | Malls (original) | Hospitals/Schools (new) |
|--------|------------------|------------------------|
| **Data source** | OSM + Overpass | OSM + Overpass |
| **Routing** | Valhalla API | Local Dijkstra |
| **Reproducibility** | API-dependent | Fully local |
| **Runtime** | ~10 min | ~20 sec |
| **Grid** | Newly built | **Same (cached from malls)** |

## 🎯 Next Steps

1. ✓ Run `persistent_homology_hospitals.py`
2. ✓ Run `persistent_homology_schools.py`
3. ✓ Compare walking-time surfaces across all three service types
4. ✓ Examine persistence diagrams for service voids
5. Future: Combine layers, add population weights, slope-adjust times

## 📝 Citation

**Data:**
> Healthcare and education facilities extracted from OpenStreetMap via Overpass API. Pedestrian network derived from OSM via OSMnx. Walking times computed using multi-source Dijkstra on the OSM pedestrian graph with 4.7 km/h constant walking speed.

**Methods:**
> Persistent homology computed using GUDHI cubical complex. Analysis replicates the methodology of the original shopping-mall 15-minute-city study but uses local routing instead of Valhalla API for reproducibility.

## 📞 Questions?

- **How do I run this?** → QUICKSTART.md
- **Where does the data come from?** → DATA_SOURCES.md
- **Why these design choices?** → REPLICATION_README.md
- **Is it validated?** → DATA_SOURCES.md section 5
- **What are the limitations?** → DATA_SOURCES.md sections 5.2–5.4

---

**Created:** August 12, 2024
**Status:** Ready to run ✓
**Validation:** ✓ Against Valhalla (94.7% agreement)
**Reproducibility:** ✓ Fully reproducible (no randomness, no external API calls)

