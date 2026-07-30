from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import random
from typing import List, Optional, Sequence, Tuple

import networkx as nx
import numpy as np


@dataclass(frozen=True)
class TUDatasetFiles:
    adjacency: Path
    graph_indicator: Path
    graph_labels: Path
    node_labels: Path
    edge_labels: Optional[Path]


def _read_int_lines(path: Path) -> List[int]:
    """Read a file containing one integer per line."""
    values: List[int] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            values.append(int(line))
    return values

def _read_edge_list(path: Path) -> List[Tuple[int, int]]:
    """Read an edge list file with either 'u, v' or 'u v' per line."""
    edges: List[Tuple[int, int]] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "," in line:
                parts = [p.strip() for p in line.split(",")]
            else:
                parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Malformed edge line in {path}: {line!r}")
            edges.append((int(parts[0]), int(parts[1])))
    return edges


def _dataset_files(dataset_dir: Path, name: str) -> TUDatasetFiles:
    adjacency = dataset_dir / f"{name}_A.txt"
    graph_indicator = dataset_dir / f"{name}_graph_indicator.txt"
    graph_labels = dataset_dir / f"{name}_graph_labels.txt"
    node_labels = dataset_dir / f"{name}_node_labels.txt"
    edge_labels = dataset_dir / f"{name}_edge_labels.txt"

    if not adjacency.exists():
        raise FileNotFoundError(f"Missing adjacency file: {adjacency}")
    if not graph_indicator.exists():
        raise FileNotFoundError(f"Missing graph_indicator file: {graph_indicator}")
    if not graph_labels.exists():
        raise FileNotFoundError(f"Missing graph_labels file: {graph_labels}")
    if not node_labels.exists():
        raise FileNotFoundError(f"Missing node_labels file: {node_labels}")

    return TUDatasetFiles(
        adjacency=adjacency,
        graph_indicator=graph_indicator,
        graph_labels=graph_labels,
        node_labels=node_labels,
        edge_labels=edge_labels if edge_labels.exists() else None,
    )


def _resolve_dataset_dir(dataset_dir: Path) -> Path:
    """Resolve a dataset directory relative to the project root if needed."""
    if dataset_dir.is_absolute():
        return dataset_dir
    project_root = Path(__file__).resolve().parent.parent
    candidate = project_root / dataset_dir
    if candidate.exists():
        return candidate
    if dataset_dir.exists():
        return dataset_dir
    return dataset_dir


def load_tu_dataset(
    dataset_dir: str | Path,
    name: str,
    *,
    relabel_nodes: bool = True,
) -> List[nx.Graph]:
    """Load a TU-format dataset into a list of NetworkX graphs.

    Parameters
    ----------
    dataset_dir : str or Path
        Directory containing the dataset files.
    name : str
        Dataset prefix (e.g., "MUTAG").
    relabel_nodes : bool, optional
        If True, relabel nodes to 0..n_i-1 per graph while preserving
        the original id in the node attribute "global_id".

    Returns
    -------
    list of networkx.Graph
        Graphs ordered by graph_id (1..N). Each graph has:
        - graph.graph["graph_id"]
        - graph.graph["label"] (graph class label)
        Nodes have attribute "label" (node label), and "global_id" when
        relabel_nodes=True.
        Edges have attribute "label" if edge labels are available.
    """
    dataset_dir = _resolve_dataset_dir(Path(dataset_dir))
    files = _dataset_files(dataset_dir, name)

    edges = _read_edge_list(files.adjacency)
    graph_indicator = _read_int_lines(files.graph_indicator)
    node_labels = _read_int_lines(files.node_labels)
    graph_labels = _read_int_lines(files.graph_labels)
    edge_labels = _read_int_lines(files.edge_labels) if files.edge_labels else None

    if len(graph_indicator) != len(node_labels):
        raise ValueError(
            "Mismatched node counts: graph_indicator has "
            f"{len(graph_indicator)} lines, node_labels has {len(node_labels)} lines."
        )
    if edge_labels is not None and len(edge_labels) != len(edges):
        raise ValueError(
            "Mismatched edge counts: adjacency has "
            f"{len(edges)} lines, edge_labels has {len(edge_labels)} lines."
        )

    num_graphs = len(graph_labels)

    # Build a mapping from global node_id -> (graph_id, local_id)
    global_to_local: dict[int, Tuple[int, int]] = {}
    local_counters = [0] * (num_graphs + 1)  # graph_id is 1-based

    for global_id, graph_id in enumerate(graph_indicator, start=1):
        if graph_id < 1 or graph_id > num_graphs:
            raise ValueError(
                f"Invalid graph_id {graph_id} for node_id {global_id}."
            )
        if relabel_nodes:
            local_id = local_counters[graph_id]
            local_counters[graph_id] += 1
        else:
            local_id = global_id
        global_to_local[global_id] = (graph_id, local_id)

    graphs = [nx.Graph() for _ in range(num_graphs)]

    # Add nodes
    for global_id, graph_id in enumerate(graph_indicator, start=1):
        graph = graphs[graph_id - 1]
        _, local_id = global_to_local[global_id]
        graph.add_node(
            local_id,
            label=node_labels[global_id - 1],
            global_id=global_id if relabel_nodes else global_id,
        )

    # Add edges
    for idx, (u, v) in enumerate(edges):
        gu = graph_indicator[u - 1]
        gv = graph_indicator[v - 1]
        if gu != gv:
            raise ValueError(
                f"Edge ({u}, {v}) connects different graphs: {gu} vs {gv}."
            )
        graph = graphs[gu - 1]
        _, lu = global_to_local[u]
        _, lv = global_to_local[v]
        if edge_labels is None:
            graph.add_edge(lu, lv)
        else:
            graph.add_edge(lu, lv, label=edge_labels[idx])

    # Add graph-level labels
    for graph_id, graph in enumerate(graphs, start=1):
        graph.graph["graph_id"] = graph_id
        graph.graph["label"] = graph_labels[graph_id - 1]

    return graphs


def load_mutag(
    dataset_dir: str | Path = "datasets/mutag",
    *,
    relabel_nodes: bool = True,
) -> List[nx.Graph]:
    """Convenience loader for MUTAG in TU format."""
    return load_tu_dataset(dataset_dir, "MUTAG", relabel_nodes=relabel_nodes)


def add_uniform_weights(
    graphs: Sequence[nx.Graph],
    *,
    seed: Optional[int] = None,
) -> List[nx.Graph]:
    """If a graph is unweighted, assign edge weights ~ Uniform(0, 1).

    Weighted graphs are returned unchanged.
    """
    rng = random.Random(seed)
    out: List[nx.Graph] = []
    for graph in graphs:
        g = graph.copy()
        if not nx.is_weighted(g, weight="weight"):
            for u, v in g.edges():
                g[u][v]["weight"] = rng.random()
        out.append(g)
    return out


def ensure_unique_edge_weights(
    graphs: Sequence[nx.Graph],
    *,
    epsilon: float = 1e-8,
) -> List[nx.Graph]:
    """Ensure all edge weights within each graph are unique.

    This adds epsilon to duplicate weights until all weights are distinct.
    """
    if epsilon <= 0:
        raise ValueError("epsilon must be positive.")

    out: List[nx.Graph] = []
    for graph in graphs:
        g = graph.copy()
        for u, v, data in g.edges(data=True):
            if "weight" not in data:
                raise ValueError(
                    "All edges must have a 'weight' attribute to enforce uniqueness."
                )

        used: set[float] = set()
        for u, v, data in g.edges(data=True):
            weight = float(data["weight"])
            if weight in used:
                candidate = weight
                k = 1
                while candidate in used:
                    candidate = weight + k * epsilon
                    k += 1
                data["weight"] = candidate
                used.add(candidate)
            else:
                used.add(weight)

        out.append(g)
    return out


_SUPPORTED_OUTLIER_METRICS = {
    "num_nodes",
    "num_edges",
    "density",
    "avg_degree",
    "num_components",
    "largest_component_ratio",
}


def _graph_metric_value(graph: nx.Graph, metric: str) -> float:
    """Compute a numeric metric used for graph outlier detection."""
    num_nodes = graph.number_of_nodes()
    num_edges = graph.number_of_edges()

    if metric == "num_nodes":
        return float(num_nodes)
    if metric == "num_edges":
        return float(num_edges)
    if metric == "density":
        if num_nodes < 2:
            return 0.0
        return float(nx.density(graph))
    if metric == "avg_degree":
        if num_nodes == 0:
            return 0.0
        return float((2.0 * num_edges) / num_nodes)

    if num_nodes == 0:
        return 0.0

    undirected = graph if not graph.is_directed() else graph.to_undirected(as_view=True)
    if metric == "num_components":
        return float(nx.number_connected_components(undirected))
    if metric == "largest_component_ratio":
        largest = max((len(component) for component in nx.connected_components(undirected)))
        return float(largest / num_nodes)

    raise ValueError(
        f"Unsupported metric {metric!r}. "
        f"Supported metrics: {sorted(_SUPPORTED_OUTLIER_METRICS)}"
    )


def _robust_z_scores(values: np.ndarray) -> np.ndarray:
    """Compute robust z-scores with MAD and stable fallbacks."""
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))

    if mad > 0.0:
        scale = 1.4826 * mad
        return (values - median) / scale

    # Fallback to IQR-based scale when MAD collapses (many repeated values).
    q1, q3 = np.percentile(values, [25.0, 75.0])
    iqr = float(q3 - q1)
    if iqr > 0.0:
        scale = iqr / 1.349
        return (values - median) / scale

    std = float(np.std(values))
    if std > 0.0:
        return (values - float(np.mean(values))) / std

    return np.zeros_like(values, dtype=float)


def remove_graph_outliers(
    graphs: Sequence[nx.Graph],
    *,
    metrics: Sequence[str] = ("num_nodes", "num_edges"),
    z_threshold: float = 3.5,
    log_scale_metrics: Sequence[str] = ("num_nodes", "num_edges"),
    return_indices: bool = False,
) -> List[nx.Graph] | Tuple[List[nx.Graph], List[int], List[int]]:
    """Remove graph outliers based on robust z-scores over graph metrics.

    Parameters
    ----------
    graphs : Sequence[nx.Graph]
        Input graphs.
    metrics : sequence of str, optional
        Metrics to score outliers. Supported values:
        "num_nodes", "num_edges", "density", "avg_degree",
        "num_components", "largest_component_ratio".
    z_threshold : float, optional
        Absolute robust z-score threshold above which a graph is marked
        as an outlier for any selected metric.
    log_scale_metrics : sequence of str, optional
        Metrics to transform with log1p before scoring. This is useful for
        heavy-tailed size metrics such as node/edge counts.
    return_indices : bool, optional
        If True, return `(filtered_graphs, kept_indices, removed_indices)`.
        Indices are 0-based and refer to the original input order.

    Returns
    -------
    list[nx.Graph] or tuple
        Filtered list of graphs, and optional index lists when
        `return_indices=True`.
    """
    if z_threshold <= 0:
        raise ValueError("z_threshold must be positive.")
    if not metrics:
        raise ValueError("metrics must contain at least one metric name.")

    unknown_metrics = sorted(set(metrics) - _SUPPORTED_OUTLIER_METRICS)
    if unknown_metrics:
        raise ValueError(
            f"Unsupported metrics: {unknown_metrics}. "
            f"Supported metrics: {sorted(_SUPPORTED_OUTLIER_METRICS)}"
        )

    unknown_log_metrics = sorted(set(log_scale_metrics) - _SUPPORTED_OUTLIER_METRICS)
    if unknown_log_metrics:
        raise ValueError(
            f"Unsupported log_scale_metrics: {unknown_log_metrics}. "
            f"Supported metrics: {sorted(_SUPPORTED_OUTLIER_METRICS)}"
        )

    if not graphs:
        if return_indices:
            return [], [], []
        return []

    log_scale_set = set(log_scale_metrics)
    outlier_mask = np.zeros(len(graphs), dtype=bool)

    for metric in metrics:
        values = np.array(
            [_graph_metric_value(graph, metric) for graph in graphs],
            dtype=float,
        )
        if metric in log_scale_set:
            values = np.log1p(values)

        z_scores = np.abs(_robust_z_scores(values))
        outlier_mask |= z_scores > z_threshold

    filtered_graphs = [
        graph for graph, is_outlier in zip(graphs, outlier_mask) if not is_outlier
    ]

    if not return_indices:
        return filtered_graphs

    kept_indices = [idx for idx, is_outlier in enumerate(outlier_mask) if not is_outlier]
    removed_indices = [idx for idx, is_outlier in enumerate(outlier_mask) if is_outlier]
    return filtered_graphs, kept_indices, removed_indices
