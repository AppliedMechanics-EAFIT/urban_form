"""Subgraph sampling utilities for cycle-containing urban graph patches."""

from __future__ import annotations

from typing import Optional

import networkx as nx
import numpy as np
from tqdm import tqdm

try:
    import osmnx as ox
except ModuleNotFoundError:
    ox = None


VALID_SAMPLING_MODES = {"radius", "n_nearest"}


def has_cycle(G: nx.Graph) -> bool:
    """Return True when a connected graph has at least one cycle (hole)."""
    if G.number_of_nodes() < 3:
        return False
    if not nx.is_connected(G):
        return False
    return not nx.is_forest(G)


def _clean_and_validate_cycle_subgraph(G_sub: nx.Graph) -> Optional[nx.Graph]:
    """Normalize a sampled graph and reject it unless it contains a cycle."""
    if G_sub.number_of_edges() == 0:
        return None

    if ox is not None:
        G_sub = ox.convert.to_undirected(G_sub)
    else:
        G_sub = G_sub.to_undirected(as_view=False)
    G_sub.remove_nodes_from(list(nx.isolates(G_sub)))

    if G_sub.number_of_edges() == 0:
        return None

    if not nx.is_connected(G_sub):
        largest_cc = max(nx.connected_components(G_sub), key=len)
        G_sub = G_sub.subgraph(largest_cc).copy()

    if not has_cycle(G_sub):
        return None

    return G_sub


def extract_local_graph_by_radius(
    G_large: nx.Graph,
    center_node,
    radius: float,
    weight: str = "length",
) -> Optional[nx.Graph]:
    """Sample a local graph using shortest-path distance cutoff from a center node."""
    nodes_within_radius = nx.single_source_dijkstra_path_length(
        G_large,
        center_node,
        cutoff=radius,
        weight=weight,
    ).keys()

    G_sub = G_large.subgraph(nodes_within_radius).copy()
    return _clean_and_validate_cycle_subgraph(G_sub)


def extract_local_graph_by_n_nearest_nodes(
    G_large: nx.Graph,
    center_node,
    n_nearest_nodes: int,
    weight: str = "length",
) -> Optional[nx.Graph]:
    """Sample a local graph with the n closest reachable nodes from a center node."""
    shortest_paths = nx.single_source_dijkstra_path_length(
        G_large,
        center_node,
        weight=weight,
    )

    if not shortest_paths:
        return None

    ordered_nodes = [
        node
        for node, _ in sorted(shortest_paths.items(), key=lambda item: (item[1], str(item[0])))
    ]

    selected_nodes = ordered_nodes[:n_nearest_nodes]
    G_sub = G_large.subgraph(selected_nodes).copy()

    return _clean_and_validate_cycle_subgraph(G_sub)


def describe_sampling_config(
    sampling_mode: str,
    radius: Optional[float] = None,
    n_nearest_nodes: Optional[int] = None,
) -> str:
    """Return a compact human-readable description for notebook logs."""
    if sampling_mode == "radius":
        return f"radius={radius} m"
    if sampling_mode == "n_nearest":
        return f"n_nearest_nodes={n_nearest_nodes}"
    return f"unknown_mode={sampling_mode}"


def sample_cycle_subgraphs_for_class(
    G_large: nx.Graph,
    class_name: str,
    n_samples: int,
    sampling_mode: str = "radius",
    radius: Optional[float] = None,
    n_nearest_nodes: Optional[int] = None,
    max_attempts_factor: int = 80,
    weight: str = "length",
    rng: Optional[np.random.Generator] = None,
):
    """Sample cycle-containing subgraphs using either radius or n-nearest strategy."""
    if sampling_mode not in VALID_SAMPLING_MODES:
        raise ValueError(
            f"Unsupported sampling_mode='{sampling_mode}'. "
            f"Choose one of {sorted(VALID_SAMPLING_MODES)}."
        )

    if sampling_mode == "radius":
        if radius is None or radius <= 0:
            raise ValueError("radius must be a positive number when sampling_mode='radius'.")
    else:
        if n_nearest_nodes is None or n_nearest_nodes < 3:
            raise ValueError(
                "n_nearest_nodes must be >= 3 when sampling_mode='n_nearest'."
            )

    all_nodes = list(G_large.nodes)
    if len(all_nodes) == 0:
        raise ValueError(f"Class '{class_name}' has no nodes in tract graph.")

    if rng is None:
        rng = np.random.default_rng()

    sampled_graphs = []
    accepted = 0
    attempts = 0
    max_attempts = n_samples * max_attempts_factor

    pbar = tqdm(total=n_samples, desc=class_name, leave=False)

    while accepted < n_samples and attempts < max_attempts:
        center_node = rng.choice(all_nodes)

        if sampling_mode == "radius":
            G_local = extract_local_graph_by_radius(
                G_large,
                center_node,
                radius=radius,
                weight=weight,
            )
        else:
            G_local = extract_local_graph_by_n_nearest_nodes(
                G_large,
                center_node,
                n_nearest_nodes=n_nearest_nodes,
                weight=weight,
            )

        attempts += 1

        if G_local is not None:
            sampled_graphs.append(G_local)
            accepted += 1
            pbar.update(1)

    pbar.close()
    if accepted < n_samples:
        print(
            f"Warning: sampled {accepted}/{n_samples} graphs for '{class_name}' "
            f"after {attempts} attempts (max_attempts={max_attempts})."
        )
    else:
        print(f"Accepted {accepted} graphs for '{class_name}' after {attempts} attempts.")

    return sampled_graphs
