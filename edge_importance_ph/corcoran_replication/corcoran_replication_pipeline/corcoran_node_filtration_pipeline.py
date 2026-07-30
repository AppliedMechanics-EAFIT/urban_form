"""
Corcoran-style node-first road hierarchy filtration pipeline.

Usage:
    python corcoran_node_filtration_pipeline.py <graphml_path_or_folder> --out <output_dir>
    python corcoran_node_filtration_pipeline.py <graphml_folder> --cities cities.txt --out <output_dir>
"""
import argparse
import ast
import re
from pathlib import Path

import gudhi
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd


HIGHWAY_ORDER = [
    "motorway",
    "trunk",
    "primary",
    "secondary",
    "tertiary",
    "unclassified",
    "residential",
]
N_STEPS = len(HIGHWAY_ORDER)
FILTRATION_TIME = {hw: i for i, hw in enumerate(HIGHWAY_ORDER)}


def slug(value):
    """Normalize city names and file stems for txt-file matching."""
    value = Path(str(value).strip()).stem.lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    value = value.strip("_")
    value = re.sub(r"_(?:osmnx|nx)_dir(?:_\d+)?$", "", value)
    return value


def read_city_filter(cities_path):
    if cities_path is None:
        return None

    cities = set()
    with open(cities_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                cities.add(slug(line))

    if not cities:
        raise ValueError(f"No city names found in {cities_path}")
    return cities


def matches_city_filter(graphml_path, city_filter):
    if city_filter is None:
        return True

    return slug(Path(graphml_path).stem) in city_filter


def normalize_highway(value):
    if isinstance(value, (list, tuple, set)):
        values = list(value)
    elif isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
            values = parsed if isinstance(parsed, list) else [value]
        except (ValueError, SyntaxError):
            values = [value]
    else:
        values = [value]

    values = [str(v) for v in values]
    valid = [v for v in values if v in FILTRATION_TIME]
    if not valid:
        return None
    return min(valid, key=lambda hw: FILTRATION_TIME[hw])


def build_filtration_graph(G):
    edges = ox.graph_to_gdfs(G, nodes=False).copy()
    edges["highway"] = edges["highway"].apply(normalize_highway)
    edges = edges[edges["highway"].isin(HIGHWAY_ORDER)].copy()

    H = nx.Graph()
    H.add_nodes_from(G.nodes())

    for (u, v, _), row in edges.iterrows():
        if u == v:
            continue

        t = FILTRATION_TIME[row["highway"]]
        if H.has_edge(u, v):
            H[u][v]["filtration"] = min(H[u][v]["filtration"], t)
        else:
            H.add_edge(u, v, filtration=t)

    return H


def compute_persistence(H):
    st = gudhi.SimplexTree()
    node_to_int = {node: i for i, node in enumerate(H.nodes())}

    for node in H.nodes():
        st.insert([node_to_int[node]], filtration=0.0)

    for u, v, data in H.edges(data=True):
        st.insert([node_to_int[u], node_to_int[v]], filtration=float(data["filtration"]))

    st.persistence(homology_coeff_field=2, persistence_dim_max=True)
    return st.persistence_intervals_in_dimension(0), st.persistence_intervals_in_dimension(1)


def replace_infinite_deaths(intervals):
    intervals = intervals.copy()
    if len(intervals) > 0:
        intervals[np.isinf(intervals)] = N_STEPS
    return intervals


def persistence_density(intervals, max_persistence=N_STEPS):
    finite = intervals[np.isfinite(intervals).all(axis=1)] if len(intervals) > 0 else np.empty((0, 2))
    i_values = np.arange(1, max_persistence + 1)

    if len(finite) == 0:
        return i_values, np.zeros_like(i_values, dtype=float)

    lengths = (finite[:, 1] - finite[:, 0]).astype(int)
    total = len(finite)
    density = np.array([np.sum(lengths == i) / total for i in i_values])
    return i_values, density


def density_dataframe(intervals, dim):
    intervals = replace_infinite_deaths(intervals)
    steps, density = persistence_density(intervals, max_persistence=N_STEPS)

    if dim == "H1":
        density = density[::-1]

    return pd.DataFrame({"step": steps, "density": density})


def save_densities(h0, h1, out_dir, name):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    h0_df = density_dataframe(h0, "H0")
    h1_df = density_dataframe(h1, "H1")

    h0_df.to_csv(out_dir / f"{name}_H0_density.csv", index=False)
    h1_df.to_csv(out_dir / f"{name}_H1_density.csv", index=False)

    print(f"  Saved {out_dir / f'{name}_H0_density.csv'}")
    print(f"  Saved {out_dir / f'{name}_H1_density.csv'}")

    return h0_df["density"].to_numpy(), h1_df["density"].to_numpy()


def save_average_densities(all_h0, all_h1, out_dir):
    if not all_h0 or not all_h1:
        return

    out_dir = Path(out_dir)
    steps = np.arange(1, N_STEPS + 1)

    for dim, densities in [("H0", all_h0), ("H1", all_h1)]:
        avg = np.mean(np.stack(densities), axis=0)
        df = pd.DataFrame({"step": steps, "density": avg})
        path = out_dir / f"average_{dim}_density.csv"
        df.to_csv(path, index=False)
        print(f"  Saved {path}")


def process_graph(graphml_path, out_dir):
    graphml_path = Path(graphml_path)
    name = slug(graphml_path.stem)
    print(f"\nProcessing: {graphml_path}")

    G = ox.load_graphml(graphml_path)
    H = build_filtration_graph(G)
    h0, h1 = compute_persistence(H)

    print(f"  Nodes: {H.number_of_nodes()}")
    print(f"  Filtered edges: {H.number_of_edges()}")

    return save_densities(h0, h1, out_dir, name)


def collect_graphml_files(input_path, city_filter):
    input_path = Path(input_path)

    if input_path.is_file():
        files = [input_path]
    elif input_path.is_dir():
        files = sorted(input_path.glob("*.graphml"))
    else:
        raise FileNotFoundError(f"Input path not found: {input_path}")

    files = [f for f in files if matches_city_filter(f, city_filter)]
    if not files:
        raise FileNotFoundError("No matching .graphml files found.")

    return files


def main():
    parser = argparse.ArgumentParser(description="Corcoran-style node filtration density pipeline")
    parser.add_argument("input", help="Path to a .graphml file or a folder containing .graphml files")
    parser.add_argument("--out", required=True, help="Folder where density CSVs will be saved")
    parser.add_argument("--cities", default=None, help="Optional txt file with city names or graph stems to process")
    args = parser.parse_args()

    city_filter = read_city_filter(args.cities)
    graphml_files = collect_graphml_files(args.input, city_filter)

    print(f"Found {len(graphml_files)} graph(s).")
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_h0, all_h1 = [], []
    for graphml_path in graphml_files:
        h0_density, h1_density = process_graph(graphml_path, out_dir)
        all_h0.append(h0_density)
        all_h1.append(h1_density)

    save_average_densities(all_h0, all_h1, out_dir)

    print(f"\nDone. Density files saved in: {out_dir}")


if __name__ == "__main__":
    main()
