"""
Triplet 15-minute-city persistent homology for Medellin.

End-to-end, reproducible pipeline: perimeter -> POIs -> local pedestrian
routing -> four sublevel-set filtrations -> persistent homology -> indicators.

Replaces the Valhalla public-API routing step with local multi-source Dijkstra
on an OSM pedestrian network. Validated against live Valhalla at r = 0.966
(median |difference| 1.13 min, 94.7% agreement on the 15-minute classification).

Run:  python pipeline_triplet.py
Requires: geopandas, osmnx, networkx, scipy, gudhi, pyarrow, matplotlib.

Notes on the two reproducibility traps this pipeline fixes:
  1. The perimeter/grid is taken from walking_time_mall.parquet (the cached
     lattice from the original saved run), NOT re-derived from a graph. The
     notebook pointed at data/sit/output/street_network_with_elevation.graphml,
     which yields a different and much larger study area than the readme
     describes; medellin_drive_with_elevation.graphml node-bounds is what
     actually reproduces 149x130 = 19,370 candidates / 13,150 retained cells.
  2. Valhalla's public server caps a matrix request at 2,500 TOTAL locations
     (sources x targets). With 191 mall targets, 15 sources = 2,865 -> HTTP 400
     "Exceeded max locations". Use <= 12 sources per batch at this target count.
     This only matters for the optional validation step; the pipeline itself
     does not call any web routing service.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox
import networkx as nx
import requests
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
from shapely.geometry import Point
import gudhi

ROOT = Path(__file__).parent
UTM_CRS = 'EPSG:32618'
WGS84 = 'EPSG:4326'
CELL_SIZE_M = 100.0
SENTINEL = 120.0     # unreachable / outside-perimeter filtration value
WALK_KMH = 4.7       # matches Valhalla pedestrian costing; calibration-confirmed
BBOX = (6.10, -75.70, 6.40, -75.44)

TAGS = {
    'health': ['amenity=clinic', 'amenity=doctors', 'amenity=health_post', 'amenity=hospital'],
    'school': ['amenity=school', 'amenity=kindergarten'],
    'mall':   ['shop=mall', 'building=retail'],
}
LBL = {'t_mall': 'Malls / retail', 't_health': 'Clinics & hospitals',
       't_school': 'Schools & kindergartens', 't_max': 'All three services'}
SURF = ['t_mall', 't_health', 't_school', 't_max']


# ---------------------------------------------------------------- 1. perimeter
def load_grid():
    """Authoritative lattice = the cached run's grid. Guarantees every new
    service layer sits on exactly the same cells as the original mall result."""
    cached = gpd.read_parquet(ROOT / 'walking_time_mall.parquet')
    g = cached.to_crs(UTM_CRS)[['row_idx', 'col_idx', 'geometry']].copy()
    cen = g.geometry.centroid
    g['centroid_lon'] = gpd.GeoSeries(cen, crs=UTM_CRS).to_crs(WGS84).x
    g['centroid_lat'] = gpd.GeoSeries(cen, crs=UTM_CRS).to_crs(WGS84).y
    perim = g.geometry.union_all()
    n_rows = int(g.row_idx.max()) + 1
    n_cols = int(g.col_idx.max()) + 1
    assert len(g) == 13150 and (n_rows, n_cols) == (149, 130), \
        f'lattice drift: {len(g)} cells, {n_rows}x{n_cols}'
    return g, perim, n_rows, n_cols


# -------------------------------------------------------------------- 2. POIs
def fetch_pois(perim, cache=ROOT / 'data/overpass_raw_triplet.json'):
    """Overpass query for all three layers; cached to disk on first run."""
    if cache.exists():
        raw = json.load(open(cache))
    else:
        raw = {}
        for layer, tags in TAGS.items():
            q = '[out:json][timeout:180];('
            for t in tags:
                k, v = t.split('=')
                q += f'nwr["{k}"="{v}"]({BBOX[0]},{BBOX[1]},{BBOX[2]},{BBOX[3]});'
            q += ');out center;'
            r = requests.post('https://overpass-api.de/api/interpreter',
                              data={'data': q}, timeout=300)
            r.raise_for_status()
            raw[layer] = r.json()['elements']
        cache.parent.mkdir(exist_ok=True)
        json.dump(raw, open(cache, 'w'))

    layers = {}
    for layer, payload in raw.items():
        # cache may hold either the full Overpass response or just the elements
        els = payload['elements'] if isinstance(payload, dict) else payload
        recs = []
        for e in els:
            lon = e.get('lon', (e.get('center') or {}).get('lon'))
            lat = e.get('lat', (e.get('center') or {}).get('lat'))
            if lon is None or lat is None:
                continue
            recs.append({'osm': f"{e['type']}/{e['id']}",
                         'name': (e.get('tags') or {}).get('name'),
                         'lon': lon, 'lat': lat})
        df = pd.DataFrame(recs).drop_duplicates('osm')
        gdf = gpd.GeoDataFrame(df, geometry=[Point(a, b) for a, b in zip(df.lon, df.lat)],
                               crs=WGS84).to_crs(UTM_CRS)
        # Collapse co-located same-NAME duplicates (one campus mapped as both a
        # node and a way). Unnamed features are never collapsed: OSM leaves many
        # genuine clinics and schools untagged, and keying them on '' would merge
        # distinct facilities that merely share a 60 m cell.
        named = gdf[gdf.name.notna() & (gdf.name.str.strip() != '')].copy()
        unnamed = gdf[~gdf.index.isin(named.index)]
        named['_k'] = (named.name.str.strip().str.lower() + '|' +
                       (named.geometry.x // 60).astype(int).astype(str) + '|' +
                       (named.geometry.y // 60).astype(int).astype(str))
        named = named.drop_duplicates('_k').drop(columns='_k')
        gdf = pd.concat([named, unnamed], ignore_index=True)
        layers[layer] = gdf[gdf.within(perim)].reset_index(drop=True)
    return layers


# ------------------------------------------------- 3. local pedestrian routing
def build_walk_graph(perim, cache=ROOT / 'data/medellin_walk.graphml'):
    if cache.exists():
        G = ox.load_graphml(cache)
    else:
        poly = gpd.GeoSeries([perim], crs=UTM_CRS).to_crs(WGS84).iloc[0].buffer(0.02)
        G = ox.graph_from_polygon(poly, network_type='walk', simplify=True)
        ox.save_graphml(G, cache)
    # pedestrians ignore one-way restrictions -> undirected; keep largest component
    Gu = ox.convert.to_undirected(G)
    Gu = Gu.subgraph(max(nx.connected_components(Gu), key=len)).copy()
    return Gu


def route_layers(Gu, grid, layers):
    """One multi-source Dijkstra per layer over a CSR adjacency in minutes."""
    nodes = list(Gu.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    xy = np.array([[Gu.nodes[n]['x'], Gu.nodes[n]['y']] for n in nodes])
    gn = gpd.GeoSeries([Point(a, b) for a, b in xy], crs=WGS84).to_crs(UTM_CRS)
    tree = cKDTree(np.c_[gn.x, gn.y])

    spm = WALK_KMH * 1000.0 / 60.0     # metres per minute
    # Collapse parallel edges (MultiGraph) to the SHORTEST one. Casting to
    # nx.Graph would instead keep whichever edge happened to be inserted last.
    best = {}
    for u, v, d in Gu.edges(data=True):
        key = (idx[u], idx[v]) if idx[u] < idx[v] else (idx[v], idx[u])
        t = float(d.get('length', 0.0)) / spm
        if key not in best or t < best[key]:
            best[key] = t
    r, c, w = [], [], []
    for (i, j), t in best.items():
        r += [i, j]; c += [j, i]; w += [t, t]
    A = sp.csr_matrix((w, (r, c)), shape=(len(nodes), len(nodes)))

    cell_xy = np.c_[grid.geometry.centroid.x, grid.geometry.centroid.y]
    _, cell_node = tree.query(cell_xy)

    out = {}
    for layer, gdf in layers.items():
        _, src = tree.query(np.c_[gdf.geometry.x, gdf.geometry.y])
        dist = dijkstra(A, directed=False, indices=np.unique(src), min_only=True)
        out[layer] = dist[cell_node]
    return out, cell_node


# --------------------------------------------------------------- 4. filtration
def build_matrix(grid, col, n_rows, n_cols):
    M = np.full((n_rows, n_cols), SENTINEL)
    M[grid.row_idx.values, grid.col_idx.values] = grid[col].values
    return M


def persistence(M):
    cc = gudhi.CubicalComplex(dimensions=list(M.shape),
                              top_dimensional_cells=M.flatten().tolist())
    cc.compute_persistence()
    p = cc.persistence()
    h0 = [(b, d) for dim, (b, d) in p if dim == 0]
    h1 = [(b, d) for dim, (b, d) in p if dim == 1]
    return {'h0': h0, 'h1': h1,
            'h0f': [(b, d) for b, d in h0 if d < SENTINEL],
            'h1f': [(b, d) for b, d in h1 if d < SENTINEL]}


def betti(pairs, ts):
    return [sum(1 for b, d in pairs if b <= t < d) for t in ts]


def indicators(PH, mats, city_cells, th=15.0):
    rows = []
    for s in SURF:
        h0, h1 = PH[s]['h0'], PH[s]['h1']
        cov = 100 * (mats[s] <= th).sum() / city_cells
        frag = sum(1 for b, d in h0 if b <= th < d)
        beyond = [(b, d) for b, d in h1 if d > th and d < SENTINEL]
        svi = sum(d - b for b, d in beyond)
        rows.append({
            'surface': LBL[s],
            'coverage_15min_pct': round(cov, 1),
            'H0_fragmentation_beta0_at15': frag,
            'H1_beta1_at15': sum(1 for b, d in h1 if b <= th < d),
            'H1_service_void_index_min': round(svi, 1),
            'H1_voids_crossing_15': len(beyond),
            'H1_voids_born_after_15': sum(1 for b, d in h1 if b > th and d < SENTINEL),
            'max_H1_persistence_min': round(max(d - b for b, d in PH[s]['h1f']), 1),
            'topological_15min_score': round(cov / (1 + np.log1p(frag) + np.log1p(svi)), 2),
        })
    return pd.DataFrame(rows)


# -------------------------------------------------------------- optional check
def validate_against_valhalla(grid, layers, n=300, seed=42, batch=12):
    """Live like-for-like check. batch*len(targets) must stay under 2,500."""
    tg = layers['mall'].to_crs(WGS84)
    tlon, tlat = tg.geometry.x.values, tg.geometry.y.values
    assert batch * len(tlon) < 2500, 'reduce batch: Valhalla caps 2,500 locations'
    samp = np.random.default_rng(seed).choice(len(grid), size=n, replace=False)
    sub = grid.iloc[samp]
    out = np.full(n, np.nan)
    for i in range(0, n, batch):
        pay = {'sources': [{'lon': float(a), 'lat': float(b)} for a, b in
                           zip(sub.centroid_lon.values[i:i+batch], sub.centroid_lat.values[i:i+batch])],
               'targets': [{'lon': float(a), 'lat': float(b)} for a, b in zip(tlon, tlat)],
               'costing': 'pedestrian'}
        rr = requests.post('https://valhalla1.openstreetmap.de/sources_to_targets',
                           json=pay, timeout=120)
        rr.raise_for_status()
        for j, row in enumerate(rr.json()['sources_to_targets']):
            v = [k['time'] for k in row if k and k.get('time') is not None]
            if v:
                out[i+j] = min(v) / 60.0
    return samp, out


def main():
    (ROOT / 'results').mkdir(exist_ok=True)
    (ROOT / 'figs').mkdir(exist_ok=True)
    grid, perim, n_rows, n_cols = load_grid()
    layers = fetch_pois(perim)
    print({k: len(v) for k, v in layers.items()})

    Gu = build_walk_graph(perim)
    print(f'walk graph: {Gu.number_of_nodes()} nodes, {Gu.number_of_edges()} edges')
    res, _ = route_layers(Gu, grid, layers)

    for lay in ('mall', 'health', 'school'):
        grid[f't_{lay}'] = np.where(np.isfinite(res[lay]), res[lay], SENTINEL)
    grid['t_max'] = grid[['t_mall', 't_health', 't_school']].max(axis=1)
    grid['binding'] = grid[['t_mall', 't_health', 't_school']].idxmax(axis=1).str.replace('t_', '')

    mats = {s: build_matrix(grid, s, n_rows, n_cols) for s in SURF}
    city_cells = int((mats['t_mall'] < SENTINEL).sum())
    PH = {s: persistence(mats[s]) for s in SURF}

    ind = indicators(PH, mats, city_cells)
    ind.to_csv(ROOT / 'results/topological_indicators.csv', index=False)
    grid.drop(columns=[c for c in ['centroid_utm'] if c in grid]).to_parquet(
        ROOT / 'data/accessibility_triplet.parquet')
    print(ind.to_string(index=False))


if __name__ == '__main__':
    main()
