#!/usr/bin/env python3
"""
Compute the ODD-STh (OS) kernel between two NetworkX graphs using GraKeL.

Notes:
- ODD-STh is a labeled-graph kernel. If your NetworkX graphs don't have node labels,
  this script will assign a default constant label to every node.
- By default, it looks for node labels under node attribute key: 'label'.
"""

from __future__ import annotations

import argparse
from typing import Any, Hashable, Optional, Tuple

import networkx as nx
import numpy as np

from grakel.kernels import OddSth
from grakel.utils import graph_from_networkx


def _ensure_node_labels(
    G: nx.Graph, node_label_attr: str = "label", default_label: Any = "1"
) -> nx.Graph:
    """
    Ensure every node has a label attribute. Returns a (shallow) copy of G.
    """
    H = G.copy()
    for n, data in H.nodes(data=True):
        if node_label_attr not in data:
            data[node_label_attr] = default_label
    return H


def os_kernel_between_two_nx_graphs(
    G1: nx.Graph,
    G2: nx.Graph,
    *,
    h: Optional[int] = None,
    normalize: bool = False,
    node_label_attr: str = "label",
) -> float:
    """
    Compute the scalar ODD-STh kernel k(G1, G2).

    Parameters
    ----------
    h : int or None
        Maximum (single) DAG height for ODD-STh. If None, unrestricted. (GraKeL parameter)
    normalize : bool
        Whether to normalize the kernel.
    node_label_attr : str
        NetworkX node attribute key holding labels.
    """
    G1 = _ensure_node_labels(G1, node_label_attr=node_label_attr)
    G2 = _ensure_node_labels(G2, node_label_attr=node_label_attr)

    # Convert NetworkX -> GraKeL graph objects/iterables
    gk_graphs = list(
        graph_from_networkx([G1, G2], node_labels_tag=node_label_attr, as_Graph=False)
    )

    # Instantiate the kernel (built-in in GraKeL)
    k = OddSth(h=h, normalize=normalize)

    # Fit on first graph, then transform second graph against it => (1x1) matrix
    K_train = k.fit_transform([gk_graphs[0]])
    K_test = k.transform([gk_graphs[1]])

    # Kernel value between G1 and G2 is the single entry
    return float(K_test[0, 0])


def _demo_graphs() -> Tuple[nx.Graph, nx.Graph]:
    """
    Small demo graphs with node labels.
    """
    G1 = nx.Graph()
    G1.add_edges_from([(0, 1), (1, 2), (2, 3)])
    nx.set_node_attributes(G1, {0: "a", 1: "b", 2: "a", 3: "c"}, "label")

    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (0, 2), (2, 3)])
    nx.set_node_attributes(G2, {0: "a", 1: "b", 2: "c", 3: "c"}, "label")

    return G1, G2


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute ODD-STh (OS) kernel using GraKeL.")
    parser.add_argument("--h", type=int, default=None, help="Maximum DAG height (OddSth parameter).")
    parser.add_argument(
        "--normalize", action="store_true", help="Normalize the kernel value (OddSth parameter)."
    )
    parser.add_argument(
        "--node-label-attr",
        type=str,
        default="label",
        help="NetworkX node attribute key to use as node label.",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run a built-in demo instead of loading your own graphs.",
    )

    args = parser.parse_args()

    if args.demo:
        G1, G2 = _demo_graphs()
    else:
        # Replace this with your own loading logic (pickle/graphml/etc.)
        # For now, we just demonstrate with unlabeled graphs (labels will be auto-added).
        G1 = nx.path_graph(5)
        G2 = nx.cycle_graph(5)

    val = os_kernel_between_two_nx_graphs(
        G1,
        G2,
        h=args.h,
        normalize=args.normalize,
        node_label_attr=args.node_label_attr,
    )

    print(val)


if __name__ == "__main__":
    main()
