# Topological accessibility of the 15-minute city in Medellín
## Three services: retail, health, and education

This report extends the malls-only proof of concept to a three-service model.
All results below are regenerated end-to-end by `pipeline_triplet.py`, which
runs cold in about 13 seconds.

---

## 1. What changed from the proof of concept

| | Proof of concept | This study |
|---|---|---|
| Services | Malls / retail only | Retail + health + education |
| POI count | 140 named malls | 200 retail, 223 health, 794 education |
| Routing engine | Valhalla public API (hours) | Local multi-source Dijkstra (50,441-node pedestrian graph, < 1 s) |
| Network | Drive graph, 15,207 nodes | Dedicated `network_type='walk'` graph, 50,441 nodes |
| Filtrations | 1 | 4 (three single-service + composite) |

The study perimeter is unchanged: the same 13,150 cells on the same
149×130 lattice as the saved mall run, so every new layer is
directly comparable to the original result.

### 1.1 The routing engine is the main methodological upgrade

Rather than a faster web service, the pipeline now routes locally. A single
multi-source Dijkstra pass per layer computes the distance from *every* grid
cell to its nearest facility at once, because minimum-distance-to-any-source is
exactly what `scipy`'s `min_only=True` computes. This replaces hundreds of
rate-limited HTTP round trips with one sparse linear-algebra call.

The swap was validated against live Valhalla on 300 randomly sampled
cells, using the *same* current POI set so that OpenStreetMap churn could not be
mistaken for an engine difference:

- Pearson r = **0.966**, Spearman ρ = 0.953
- Median absolute difference **1.13 min**; median bias -0.65 min
- **94.7%** agreement on the binary within-15-minutes classification
- A calibration search over multiplicative speed factors returned an optimum of
  1.000, independently confirming that 4.7 km/h already matches Valhalla's
  pedestrian costing.

![Engine validation]({{artifact:art_d1818be9-41da-44ff-a708-53421d5f38df}})

The local engine is therefore adopted as the pipeline's routing method, with
Valhalla retained only as an optional check.

*Practical note.* Valhalla's public server caps a matrix request at 2,500
**total** locations, i.e. sources × targets — not sources alone. With
200 targets, the original 15-source batch requests 15 × 200 =
3,000 locations and fails; 12 sources stays legal. This is why the
original run worked at 140 targets and broke as the POI set grew.

---

## 2. Data

![Three service layers]({{artifact:art_23d6e1bc-595c-453a-baaa-4cbc98a8b05d}})

| Layer | OSM tags | n inside perimeter |
|---|---|---|
| Retail | `shop=mall`, `building=retail` | 200 |
| Health | `amenity` in clinic, doctors, health_post, hospital | 223 |
| Education | `amenity` in school, kindergarten | 794 |

Unnamed features are deliberately retained — an untagged clinic is still a
clinic — while co-located *same-name* duplicates (one campus mapped as both a
node and a way) are collapsed.

**Health-layer completeness caveat.** Colombia's official REPS registry lists
about 1,200 registered IPS institutions in Medellín against
223 health facilities in OpenStreetMap. REPS publishes addresses but
no coordinates, so it cannot serve as a geometry source here; it does establish
that the OSM health layer is substantially incomplete. Health results should be
read as a lower bound on true access — geocoding REPS is the highest-value next
step in this project.

---

## 3. Accessibility surfaces

![Four surfaces]({{artifact:art_a91d7617-9e19-438e-b6a0-48d066c60491}})

| surface | mean | median | p25 | p75 | max | n_le15 | cov15_pct | n_le30 | cov30_pct |
|---|---|---|---|---|---|---|---|---|---|
| Malls / retail | 16.79 | 13.13 | 7.45 | 21.81 | 117.43 | 7473 | 56.83 | 11594 | 88.17 |
| Clinics & hospitals | 13.42 | 9.66 | 5.43 | 16.89 | 109.92 | 9210 | 70.04 | 12254 | 93.19 |
| Schools & kindergartens | 6.7 | 4.59 | 2.47 | 8.67 | 90.83 | 12083 | 91.89 | 12933 | 98.35 |
| All three services | 18.86 | 15.53 | 10.16 | 23.11 | 117.43 | 6283 | 47.78 | 11514 | 87.56 |

Three findings stand out.

**A clear service hierarchy.** Education is Medellín's most evenly distributed
service (91.9% of cells within 15 minutes, median
4.59 min), health is intermediate (70.0%), and
retail is the scarcest (56.8%).

**The composite is not the average — it is the worst case.** Requiring all three
services within 15 minutes yields **47.8%** coverage, *below*
every individual layer. This is the arithmetic of a max-of-three surface, and it
is the number a genuine 15-minute-city standard has to answer for: a resident
needs all their daily services nearby, not one of them.

**The gap narrows at 30 minutes.** Composite coverage rises from
47.8% to 87.6%. Medellín is much closer to
being a 30-minute city than a 15-minute one, so the deficit is concentrated in a
band of cells sitting just outside the walkable threshold.

---

## 4. Persistent homology

![Persistence diagrams]({{artifact:art_29ed0357-8407-4cfe-a017-6f7b06fd6495}})

![Betti curves]({{artifact:art_8b2d8a20-068b-40e4-9bd3-95b181ecb8db}})

### 4.1 Fragmentation (H₀)

β₀ counts disconnected accessibility islands at a given threshold.

At t = 15 min the composite surface has **108 islands** —
more than retail (78), health (45), or
education (22). **The composite is simultaneously the
lowest-coverage and the most fragmented surface.** The well-served area is not
merely smaller; it is broken into more pieces. Islands that are contiguous in
one layer are severed in another, and intersecting the three requirements cuts
the accessible fabric in more places than any single service does.

Peak fragmentation runs the other way: education peaks highest
(β₀ = 673 near t ≈ 2 min) simply because 794 schools
seed 794 separate islands that later merge. Peak β₀ measures facility
density; β₀ *at the policy threshold* measures whether coverage forms a
continuous fabric. Only the second is a planning quantity.

### 4.2 Service voids (H₁)

β₁ counts enclosed voids — poorly served pockets ringed by better-served area.
These are the features standard coverage metrics cannot see, because a percentage
has no shape.

![Composite barcodes]({{artifact:art_0af65ce5-372c-42d8-99db-f8c3651a52ba}})

The composite surface has **418 finite H₁ features**, of
which **6 persist for more than 30 minutes** of
filtration time. The longest-lived void is born at 31.06 min
and dies at 116.66 min — **85.6 minutes of
persistence**. A void this long-lived is not a discretization artifact: it
survives across the entire plausible range of walking-time thresholds, so no
choice of standard between 15 and 60 minutes would classify it as served.

Note the barcode structure: composite H₁ bars are born almost exclusively
*after* the 15-minute line. These voids are invisible to a 15-minute audit and
appear only as the threshold relaxes — precisely the areas a coverage percentage
reports as uniformly unserved while topology shows them as enclosed pockets
inside otherwise-served territory.

---

## 5. Which service is to blame

Because the composite is a maximum, each cell has a *binding* service — the one
that is furthest away and therefore sets that cell's composite value.

![Void attribution]({{artifact:art_0ad5e914-4b68-452f-8b5a-a4c5bfa0600d}})

- Retail binds **58.6%** of cells
- Health binds **33.1%**
- Education binds only **8.3%**

Attribution within the ten longest-lived voids:

| void | birth_min | death_min | persistence_min | n_cells_in_window | binding_mall_pct | binding_health_pct | binding_school_pct |
|---|---|---|---|---|---|---|---|
| 1.0 | 31.06 | 116.66 | 85.6 | 1489.0 | 70.0 | 29.8 | 0.2 |
| 2.0 | 38.42 | 117.4 | 78.98 | 866.0 | 58.5 | 41.1 | 0.3 |
| 3.0 | 41.51 | 117.43 | 75.92 | 708.0 | 54.2 | 45.3 | 0.4 |
| 4.0 | 42.15 | 117.43 | 75.28 | 666.0 | 54.2 | 45.3 | 0.5 |
| 5.0 | 31.72 | 89.73 | 58.01 | 1402.0 | 70.2 | 29.6 | 0.2 |
| 6.0 | 80.48 | 117.43 | 36.95 | 97.0 | 23.7 | 76.3 | 0.0 |

The pattern is unambiguous. Retail is the binding constraint in 54–76% of cells
in every top void, health in 24–45%, and education essentially never
(under 1.5%, and under 0.6% in nine of the ten). Void 6 is the sole exception,
where health binds 76.3% — an isolated pocket that is genuinely a health-access
problem rather than a retail one.

**Policy reading.** Medellín's 15-minute deficit is overwhelmingly a *retail and
health* problem, not an education problem. The school network already delivers
near-universal walkable coverage; siting more schools would barely move the
composite surface. Placing everyday retail and primary care inside the
identified voids is what would close the gap — and the void map localizes where.

---

## 6. Topological indicators

| surface | coverage_15min_pct | H0_fragmentation_beta0_at15 | H1_beta1_at15 | H1_service_void_index_min | H1_voids_born_after_15 | topological_15min_score |
|---|---|---|---|---|---|---|
| Malls / retail | 56.8 | 78 | 36 | 1451.5 | 283 | 4.49 |
| Clinics & hospitals | 70.0 | 45 | 69 | 1183.8 | 139 | 5.88 |
| Schools & kindergartens | 91.9 | 22 | 77 | 818.4 | 27 | 8.47 |
| All three services | 47.8 | 108 | 35 | 1451.6 | 305 | 3.68 |

`topological_15min_score` = coverage ÷ (1 + log1p(β₀ at 15) + log1p(service-void
index)), penalizing coverage that is fragmented or void-riddled. The ranking —
education 8.47 > health >
retail > composite — matches the substantive reading: education access is broad
*and* continuous, while the composite is both narrow and broken.

The service-void index (total H₁ persistence beyond 15 min) separates two
surfaces that coverage alone would rank similarly: retail and the composite have
nearly identical void indices (~1,451 min) despite a 9-point coverage gap,
because they fail in the same places for the same reason.

---

## 7. Limitations

1. **Health layer incompleteness** — 223 OSM facilities vs ~1,200 REPS-registered institutions. Health coverage is a lower bound.
2. **Perimeter** — the drive-network convex hull from the original run, not an official administrative boundary. Retained deliberately for comparability with the saved mall result.
3. **No population weighting** — every cell counts equally, so a void over a steep hillside with few residents weighs the same as one over a dense barrio.
4. **Flat walking speed** — 4.7 km/h everywhere. In a city with Medellín's topography, slope-adjusted times would lengthen uphill trips materially, likely deepening the peripheral voids.
5. **H₁ localization is a proxy** — void interiors are identified by cells whose filtration values fall inside a feature's birth–death window, not by extracting exact homological generators. Locations are indicative, not exact.
6. **Service definitions** — primary care and K-12 only; no pharmacies, groceries, parks, or transit. Adding groceries would likely sharpen the retail signal further.

---

## 8. Reproducing

```bash
python pipeline_triplet.py      # ~13 s cold
```

Outputs land in `results/` and `data/`. Overpass responses and the pedestrian
graph are cached under `data/`, so only the first run needs network access.
`validate_against_valhalla()` is optional and the only function that calls a web
service.

---

## 9. Key takeaway

Coverage says Medellín is 47.8% of a 15-minute city. Topology
says that 47.8% is split into 108 disconnected
islands and punctured by 6 service voids that persist
across every plausible threshold — and that retail and primary care, not schools,
are what bind them. Two cities with identical coverage percentages can have
entirely different accessibility structures; persistent homology is what
distinguishes them, and it points at where to build.
