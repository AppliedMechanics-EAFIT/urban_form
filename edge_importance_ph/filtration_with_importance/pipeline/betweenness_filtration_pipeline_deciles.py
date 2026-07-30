"""
Betweenness Filtration Pipeline
Usage:
    python betweenness_filtration_pipeline_deciles.py <graphml_path_or_folder> [--plot] [--out <output_dir>]
"""
import argparse
import os
from pathlib import Path

import osmnx as ox
import networkx as nx
import numpy as np
import pandas as pd
import gudhi
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

N_STEPS = 10
FINE_STEPS = 5
COARSE_STEPS = N_STEPS - FINE_STEPS
TOP_FRACTION = 0.15

def build_filtration_graph(G):
    G_undir = ox.convert.to_undirected(G)

    raw_values = []
    for u, v, data in G_undir.edges(data=True):
        val = data.get("betweenness", 0)
        try:
            raw_values.append(float(val))
        except (ValueError, TypeError):
            raw_values.append(0.0)
    raw_values = np.array(raw_values)

    log_values = np.log1p(raw_values)

    fine_fraction = TOP_FRACTION / FINE_STEPS
    coarse_fraction = (1.0 - TOP_FRACTION) / COARSE_STEPS
    step_fractions = np.array(
        [fine_fraction] * FINE_STEPS + [coarse_fraction] * COARSE_STEPS,
        dtype=float,
    )

    ascending_quantiles = np.concatenate([
        np.linspace(0.0, 1.0 - TOP_FRACTION, COARSE_STEPS + 1),
        np.linspace(1.0 - TOP_FRACTION + fine_fraction, 1.0, FINE_STEPS),
    ])
    all_edges = np.quantile(log_values, ascending_quantiles)
    top20_lo = all_edges[COARSE_STEPS]

    edge_steps = np.empty(len(log_values), dtype=int)
    order = np.argsort(-log_values, kind="mergesort")
    bounds = np.rint(np.concatenate([[0.0], np.cumsum(step_fractions)]) * len(order)).astype(int)
    bounds[0] = 0
    bounds[-1] = len(order)

    for step in range(1, N_STEPS + 1):
        start = bounds[step - 1]
        stop = bounds[step]
        edge_steps[order[start:stop]] = step

    H = nx.Graph()
    H.add_nodes_from(G_undir.nodes())
    for (u, v, data), step in zip(G_undir.edges(data=True), edge_steps):
        step = int(step)
        if H.has_edge(u, v):
            H[u][v]["filtration"] = min(H[u][v]["filtration"], step)
            H[u][v]["betweenness"] = max(H[u][v]["betweenness"], float(data.get("betweenness", 0)))
        else:
            H.add_edge(u, v, filtration=step, betweenness=float(data.get("betweenness", 0)))

    return H, G_undir, raw_values, log_values, all_edges, top20_lo


def compute_persistence(H):
    st = gudhi.SimplexTree()
    node_to_int = {node: i for i, node in enumerate(H.nodes())}
    for node in H.nodes():
        st.insert([node_to_int[node]], filtration=0.0)
    for u, v, data in H.edges(data=True):
        st.insert([node_to_int[u], node_to_int[v]], filtration=float(data["filtration"]))
    st.persistence(homology_coeff_field=2, persistence_dim_max=True)
    h0 = st.persistence_intervals_in_dimension(0)
    h1 = st.persistence_intervals_in_dimension(1)
    return h0, h1


def persistence_density(intervals, max_persistence=None):
    finite = intervals[np.isfinite(intervals).all(axis=1)] if len(intervals) > 0 else np.empty((0, 2))
    if len(finite) == 0:
        return np.array([]), np.array([])
    lengths = (finite[:, 1] - finite[:, 0]).astype(int)
    total = len(finite)
    if max_persistence is None:
        max_persistence = int(lengths.max())
    i_values = np.arange(1, max_persistence + 1)
    density = np.array([np.sum(lengths == i) / total for i in i_values])
    return i_values, density


def compute_alpha_indices(H):
    rows = []
    alpha_values = []
    for step in range(1, N_STEPS + 1):
        active_edges = [(u, v) for u, v, d in H.edges(data=True) if d["filtration"] <= step]
        G_step = nx.Graph()
        G_step.add_edges_from(active_edges)

        e = len(active_edges)
        v = G_step.number_of_nodes()
        p = nx.number_connected_components(G_step) if v > 0 else 0
        denominator = 2 * v - 5
        alpha = np.nan if denominator == 0 else (e - v + p) / denominator

        rows.append({"step": step, "edges": e, "vertices": v, "components": p, "alpha": alpha})
        alpha_values.append(alpha)

    return pd.DataFrame(rows), np.array(alpha_values, dtype=float)


def save_densities(h0, h1, out_dir, name):
    replacement = N_STEPS + 1
    h0_plot = h0.copy(); h0_plot[np.isinf(h0_plot)] = replacement
    h1_plot = h1.copy(); h1_plot[np.isinf(h1_plot)] = replacement

    densities = {}
    for intervals, dim in [(h0_plot, "H0"), (h1_plot, "H1")]:
        i_vals, dens = persistence_density(intervals, max_persistence=N_STEPS)
        if dim == "H1" and len(dens) > 0:
            dens = dens[::-1]
        if len(i_vals) == 0:
            df = pd.DataFrame(columns=["step", "density"])
            dens = np.zeros(N_STEPS)
        else:
            df = pd.DataFrame({"step": i_vals, "density": dens})
        csv_path = os.path.join(out_dir, f"{name}_{dim}_density.csv")
        df.to_csv(csv_path, index=False)
        print(f"  Saved {csv_path}")
        densities[dim] = dens
    return densities["H0"], densities["H1"]


def save_alpha_indices(alpha_df, out_dir, name):
    csv_path = os.path.join(out_dir, f"{name}_alpha_index.csv")
    alpha_df.to_csv(csv_path, index=False)
    print(f"  Saved {csv_path}")


def save_average_densities(all_h0, all_h1, all_alpha, out_dir):
    i_vals = np.arange(1, N_STEPS + 1)
    for dim, all_dens in [("H0", all_h0), ("H1", all_h1)]:
        avg = np.mean(np.stack(all_dens), axis=0)
        df = pd.DataFrame({"step": i_vals, "density": avg})
        csv_path = os.path.join(out_dir, f"average_{dim}_density.csv")
        df.to_csv(csv_path, index=False)
        print(f"  Saved {csv_path}")

    avg_alpha = np.nanmean(np.stack(all_alpha), axis=0)
    df = pd.DataFrame({"step": i_vals, "alpha": avg_alpha})
    csv_path = os.path.join(out_dir, "average_alpha_index.csv")
    df.to_csv(csv_path, index=False)
    print(f"  Saved {csv_path}")


def save_plots(G, H, G_undir, raw_values, log_values, all_edges, top20_lo, h0, h1, alpha_df, out_dir, name):
    replacement = N_STEPS + 1
    h0_plot = h0.copy(); h0_plot[np.isinf(h0_plot)] = replacement
    h1_plot = h1.copy(); h1_plot[np.isinf(h1_plot)] = replacement

    # 1. Full graph
    fig, ax = ox.plot_graph(G, node_size=0, edge_linewidth=0.4, bgcolor="white",
                            edge_color="#333333", show=False, close=False)
    ax.set_title("Street network (full graph)", fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{name}_full_graph.pdf"))
    plt.close()

    # 2. Betweenness distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    axes[0].hist(raw_values, bins=60, color="#555555", edgecolor="none")
    axes[0].set_xlabel("Raw betweenness"); axes[0].set_ylabel("Edge count")
    axes[0].set_title("Raw betweenness (right-skewed)")
    axes[1].hist(log_values, bins=60, color="#555555", edgecolor="none")
    for t in all_edges[1:-1]:
        color = "red" if t >= top20_lo else "orange"
        axes[1].axvline(t, color=color, linewidth=0.8, linestyle="--", alpha=0.8)
    axes[1].axvline(top20_lo, color="blue", linewidth=1.2, linestyle="-", alpha=0.9,
                    label="top 20% boundary")
    axes[1].set_xlabel("log1p(betweenness)"); axes[1].set_ylabel("Edge count")
    axes[1].set_title("Rank quantiles: red=fine (steps 1-4), orange=coarse (steps 5-10)")
    axes[1].legend(fontsize=8)
    plt.suptitle("Betweenness distribution: raw vs log-scale with refined rank-quantile boundaries", fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{name}_betweenness_distribution.pdf"))
    plt.close()

    # 3. Filtration snapshots (GeoDataFrame)
    nodes_gdf, edges_gdf = ox.graph_to_gdfs(G)
    xmin, ymin, xmax, ymax = edges_gdf.total_bounds
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    axes = axes.flatten()
    for step in range(1, N_STEPS + 1):
        ax = axes[step - 1]
        active_edges = {(u, v) for u, v, d in H.edges(data=True) if d["filtration"] <= step}
        mask = edges_gdf.index.map(lambda idx: (idx[0], idx[1]) in active_edges)
        subset = edges_gdf[mask]
        if not subset.empty:
            subset.plot(ax=ax, linewidth=0.4, color="#333333")
        ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
        ax.set_title(f"Step {step}\n({len(active_edges)} edges)", fontsize=9)
        ax.set_axis_off()
    plt.suptitle("Betweenness filtration snapshots (steps 1–10)", fontsize=12, y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{name}_filtration_snapshots.pdf"))
    plt.close()

    # 4. Persistence diagrams
    for intervals, dim, title in [
        (h0_plot, "H0", "Persistence diagram — H0"),
        (h1_plot, "H1", "Persistence diagram — H1"),
    ]:
        if len(intervals) > 0:
            diag = [(int(dim[1]), (b, d)) for b, d in intervals]
            plt.figure(figsize=(7, 6))
            gudhi.plot_persistence_diagram(diag)
            plt.title(f"{title} (betweenness filtration)")
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"{name}_{dim}_diagram.pdf"))
            plt.close()

    # 5. Persistence barcodes
    for intervals, dim, title in [
        (h0_plot, "H0", "Persistence barcode — H0"),
        (h1_plot, "H1", "Persistence barcode — H1"),
    ]:
        if len(intervals) > 0:
            diag = [(int(dim[1]), (b, d)) for b, d in intervals]
            plt.figure(figsize=(10, 4))
            gudhi.plot_persistence_barcode(diag)
            plt.title(f"{title} (betweenness filtration)")
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"{name}_{dim}_barcode.pdf"))
            plt.close()

    # 6. Persistence densities bar chart
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (intervals, label, dim) in zip(axes, [
        (h0_plot, "Connected component (H0)", "H0"),
        (h1_plot, "Cycle (H1)", "H1"),
    ]):
        i_vals, dens = persistence_density(intervals, max_persistence=N_STEPS)
        if dim == "H1" and len(dens) > 0:
            dens = dens[::-1]
        if len(i_vals) == 0:
            ax.set_title(f"No finite {dim} intervals"); continue
        ax.bar(i_vals, dens, color="#555555", width=0.6)
        ax.set_xlim(0.5, N_STEPS + 0.5); ax.set_ylim(0, 1.0)
        ax.set_xticks(i_vals)
        ax.set_xlabel("Death step" if dim == "H0" else "Birth step of cycle")
        ax.set_ylabel("Density τ(i)"); ax.set_title(label)
    plt.suptitle("Persistence densities (betweenness filtration)", fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{name}_persistence_densities.pdf"))
    plt.close()

    # 7. Alpha index by filtration step
    plt.figure(figsize=(7, 4))
    plt.plot(alpha_df["step"], alpha_df["alpha"], marker="o", color="#555555", linewidth=1.5)
    plt.xticks(np.arange(1, N_STEPS + 1))
    plt.xlabel("Filtration step")
    plt.ylabel("Alpha index")
    plt.title("Alpha index by betweenness filtration step")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{name}_alpha_index.pdf"))
    plt.close()

    # 8. Filtration with nodes
    pos = {node: (data["x"], data["y"]) for node, data in G.nodes(data=True)
           if "x" in data and "y" in data}
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    axes = axes.flatten()
    for step in range(1, N_STEPS + 1):
        ax = axes[step - 1]
        edges_t = [(u, v) for u, v, d in H.edges(data=True) if d["filtration"] <= step]
        G_t = nx.Graph(); G_t.add_nodes_from(H.nodes()); G_t.add_edges_from(edges_t)
        degrees = dict(G_t.degree())
        isolated = [n for n, d in degrees.items() if d == 0 and n in pos]
        connected = [n for n, d in degrees.items() if d > 0 and n in pos]
        valid_edges = [(u, v) for u, v in G_t.edges() if u in pos and v in pos]
        nx.draw_networkx_edges(G_t, pos, edgelist=valid_edges, ax=ax,
                               edge_color="#333333", width=0.25, alpha=0.7)
        if isolated:
            nx.draw_networkx_nodes(G_t, pos, nodelist=isolated, ax=ax,
                                   node_size=0.1, node_color="#cccccc", alpha=0.6)
        if connected:
            deg_vals = np.array([degrees[n] for n in connected], dtype=float)
            nx.draw_networkx_nodes(G_t, pos, nodelist=connected, ax=ax,
                                   node_size=0.1, node_color=deg_vals, cmap="plasma", alpha=0.9)
        ax.set_title(
            f"Step {step}\nisolated={len(isolated)}  connected={len(connected)}  edges={len(valid_edges)}",
            fontsize=9)
        ax.set_axis_off()
    legend_handles = [
        mpatches.Patch(facecolor="#cccccc", edgecolor="grey", label="Isolated node"),
        mpatches.Patch(facecolor="#f0c040", edgecolor="grey", label="Connected node (colour = degree)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2, fontsize=10, bbox_to_anchor=(0.5, -0.02))
    plt.suptitle("Betweenness filtration with nodes", fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{name}_filtration_with_nodes.pdf"))
    plt.close()

    print(f"  Saved plots to {out_dir}/")


def process_graph(graphml_path, out_dir, plot):
    name = Path(graphml_path).stem
    print(f"\nProcessing: {graphml_path}")
    G = ox.load_graphml(graphml_path)
    H, G_undir, raw_values, log_values, all_edges, top20_lo = build_filtration_graph(G)
    h0, h1 = compute_persistence(H)
    alpha_df, alpha_values = compute_alpha_indices(H)
    os.makedirs(out_dir, exist_ok=True)
    dens_h0, dens_h1 = save_densities(h0, h1, out_dir, name)
    save_alpha_indices(alpha_df, out_dir, name)
    if plot:
        save_plots(G, H, G_undir, raw_values, log_values, all_edges, top20_lo, h0, h1, alpha_df, out_dir, name)
    return dens_h0, dens_h1, alpha_values


def main():
    parser = argparse.ArgumentParser(description="Betweenness filtration pipeline")
    parser.add_argument("input", help="Path to a .graphml file or a folder containing .graphml files")
    parser.add_argument("--plot", action="store_true", help="Save all plots as PDFs")
    parser.add_argument("--out", default=None, help="Output directory (default: next to input)")
    args = parser.parse_args()

    input_path = Path(args.input)

    if input_path.is_dir():
        graphml_files = list(input_path.glob("*.graphml"))
        if not graphml_files:
            print("No .graphml files found in the specified folder.")
            return
        out_dir = args.out or str(input_path / "filtration_results")
        all_h0, all_h1, all_alpha = [], [], []
        for f in graphml_files:
            dens_h0, dens_h1, alpha_values = process_graph(str(f), out_dir, args.plot)
            all_h0.append(dens_h0)
            all_h1.append(dens_h1)
            all_alpha.append(alpha_values)
        save_average_densities(all_h0, all_h1, all_alpha, out_dir)
        print(f"\nSaved average densities and alpha index for {len(graphml_files)} graphs.")
    elif input_path.is_file():
        out_dir = args.out or str(input_path.parent / "filtration_results")
        process_graph(str(input_path), out_dir, args.plot)
    else:
        print(f"Input path not found: {input_path}")


if __name__ == "__main__":
    matplotlib.use("Agg")  # non-interactive backend for saving files
    main()
