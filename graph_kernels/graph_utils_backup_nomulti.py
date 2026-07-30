
from typing import Hashable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.linalg import sqrtm


def plot_weighted_graph(G: nx.Graph, seed: int = 42) -> None:
    """Plot an undirected weighted NetworkX graph with edge weights.

    Parameters
    ----------
    G : networkx.Graph
        Undirected weighted graph with numeric "weight" edge attributes.
    seed : int, optional
        Seed for layout reproducibility.
    """
    pos = nx.spring_layout(G, seed=seed)

    nx.draw(
        G,
        pos,
        with_labels=False,
        node_size=20,
        node_color="lightblue",
        font_size=10,
        edge_color="gray"
    )

    
    edge_labels = nx.get_edge_attributes(G, "weight")
    '''
    nx.draw_networkx_edge_labels(
        G,
        pos,
        edge_labels=edge_labels,
        font_size=10
    )
    '''
    plt.show()

def _deterministic_orient(
    u: Hashable,
    v: Hashable,
) -> Tuple[Hashable, Hashable]:
    """Return a stable orientation using node ordering as fallback."""
    if u == v:
        return (u, v)
    try:
        return (u, v) if u < v else (v, u)
    except TypeError:
        return (u, v) if repr(u) < repr(v) else (v, u)


def _undirected_edge_key(u: Hashable, v: Hashable) -> frozenset:
    """Hashable key for undirected edge membership checks."""
    return frozenset((u, v))


def _other_tree_incident_length(
    T: nx.Graph,
    endpoint: Hashable,
    excluded_neighbor: Hashable,
) -> float:
    """Length of the shortest tree edge at endpoint excluding endpoint-excluded_neighbor.

    In the canonical-orientation statement this value is unique ("the other edge").
    When degree > 2 in the tree, we use the shortest candidate to keep orientation
    deterministic.
    """
    candidates = [
        float(T[endpoint][nbr]["weight"])
        for nbr in T.neighbors(endpoint)
        if nbr != excluded_neighbor
    ]
    if not candidates:
        return 0.0
    return min(candidates)


def orient_edge(
    u: Hashable,
    v: Hashable,
    *,
    T: nx.Graph | None = None,
    tree_edge: bool | None = None,
) -> Tuple[Hashable, Hashable]:
    """Return an oriented edge.

    Parameters
    ----------
    u, v : Hashable
        Endpoints of the edge.
    T : networkx.Graph, optional
        Spanning tree used for canonical orientation.
    tree_edge : bool, optional
        Whether (u, v) is a tree edge in T.

    Returns
    -------
    tuple
        Oriented edge (source, target).

    Notes
    -----
    If T/tree_edge are provided, apply canonical orientation:
    1) Tree edge e=(u,v): compute omega(u), omega(v) as lengths of the
       other incident tree edges; orient low->high.
    2) Non-tree edge e=(u,v): compute omega(u), omega(v) as lengths of the
       first tree edges on the unique u-v path in T; orient low->high.
    Ties (non-generic lengths) fall back to deterministic node ordering.
    """
    if T is None or tree_edge is None:
        return _deterministic_orient(u, v)

    if tree_edge:
        omega_u = _other_tree_incident_length(T, endpoint=u, excluded_neighbor=v)
        omega_v = _other_tree_incident_length(T, endpoint=v, excluded_neighbor=u)
    else:
        if u == v:
            return (u, v)
        path = nx.shortest_path(T, source=u, target=v)
        omega_u = float(T[path[0]][path[1]]["weight"])
        omega_v = float(T[path[-1]][path[-2]]["weight"])

    if omega_u < omega_v:
        return (u, v)
    if omega_u > omega_v:
        return (v, u)
    return _deterministic_orient(u, v)

def compute_mst_and_edge_order(
    G: nx.Graph,
) -> Tuple[
    nx.Graph,
    List[Tuple[Hashable, Hashable]],
    List[Tuple[Hashable, Hashable]],
    List[Tuple[Hashable, Hashable, float]],
]:
    """Compute a minimum spanning tree and a deterministic edge ordering.

    Parameters
    ----------
    G : networkx.Graph
        Undirected weighted graph with numeric "weight" edge attributes.

    Returns
    -------
    T : networkx.Graph
        Minimum spanning tree of G.
    ordered_edges : list of tuple
        Oriented edges e1,...,em with tree edges first, then non-tree edges.
    tree_edges : list of tuple
        First n-1 canonically oriented edges (sorted by ascending weight).
    non_tree_edges : list of tuple
        Remaining g edges as (u, v, weight), canonically oriented and sorted
        by ascending weight.
    """
    T = nx.minimum_spanning_tree(G, weight="weight")

    tree_edges = [orient_edge(u, v, T=T, tree_edge=True) for u, v in T.edges()]
    tree_edges.sort(
        key=lambda e: (float(T[e[0]][e[1]]["weight"]), repr(e[0]), repr(e[1]))
    )

    tree_edge_keys = {_undirected_edge_key(u, v) for u, v in T.edges()}

    non_tree_edges = []
    for u, v, d in G.edges(data=True):
        if _undirected_edge_key(u, v) in tree_edge_keys:
            continue
        ou, ov = orient_edge(u, v, T=T, tree_edge=False)
        non_tree_edges.append((ou, ov, float(d["weight"])))

    non_tree_edges.sort(key=lambda x: (x[2], repr(x[0]), repr(x[1])))

    ordered_edges = tree_edges + [(u, v) for u, v, _ in non_tree_edges]

    return T, ordered_edges, tree_edges, non_tree_edges

def fundamental_cycle(
    T: nx.Graph,
    edge: Tuple[Hashable, Hashable],
) -> List[Tuple[Hashable, Hashable]]:
    """Return the oriented edge list for the fundamental cycle of an edge.

    Parameters
    ----------
    T : networkx.Graph
        A tree (typically an MST of the original graph).
    edge : tuple
        Non-tree edge (u, v) that closes a unique cycle with T.

    Returns
    -------
    list of tuple
        Oriented edges forming the fundamental cycle.
    """
    u, v = edge
    path = nx.shortest_path(T, source=u, target=v)

    cycle = []
    for i in range(len(path) - 1):
        cycle.append(orient_edge(path[i], path[i+1], T=T, tree_edge=True))

    cycle.append(orient_edge(u, v, T=T, tree_edge=False))
    return cycle

def cycle_edge_incidence_matrix(
    G: nx.Graph,
) -> Tuple[np.ndarray, List[Tuple[Hashable, Hashable]]]:
    """Build the cycle-edge incidence matrix for the fundamental cycles.

    Parameters
    ----------
    G : networkx.Graph
        Undirected weighted graph with numeric "weight" edge attributes.

    Returns
    -------
    M : numpy.ndarray
        Binary matrix of shape (g, m) marking edges in each fundamental cycle.
    ordered_edges : list of tuple
        Oriented edges e1,...,em consistent with the columns of M.
    """
    T, ordered_edges, tree_edges, non_tree_edges = compute_mst_and_edge_order(G)

    g = len(non_tree_edges)
    m = len(ordered_edges)

    edge_index = {e: i for i, e in enumerate(ordered_edges)}

    M = np.zeros((g, m))

    for i, (u, v, _) in enumerate(non_tree_edges):
        cycle = fundamental_cycle(T, (u, v))

        for e in cycle:
            j = edge_index[e]
            M[i, j] = 1  # orientation fixed, matches paper example

    return M, ordered_edges

def length_matrices(
    G: nx.Graph,
    ordered_edges: Sequence[Tuple[Hashable, Hashable]],
    tree_size: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute diagonal length matrices for tree and non-tree edges.

    Parameters
    ----------
    G : networkx.Graph
        Undirected weighted graph with numeric "weight" edge attributes.
    ordered_edges : sequence of tuple
        Oriented edges in the order e1,...,em.
    tree_size : int
        Number of tree edges (typically n - 1).

    Returns
    -------
    L_T : numpy.ndarray
        Diagonal matrix of tree edge lengths.
    L_g : numpy.ndarray
        Diagonal matrix of non-tree edge lengths.
    """
    lengths = []
    for u, v in ordered_edges:
        lengths.append(G[u][v]["weight"])

    L = np.diag(lengths)

    L_T = L[:tree_size, :tree_size]
    L_g = L[tree_size:, tree_size:]

    return L_T, L_g

def compute_Q(G: nx.Graph) -> np.ndarray:
    """Compute the Q matrix for the graph's fundamental cycle space.

    Parameters
    ----------
    G : networkx.Graph
        Undirected weighted graph with numeric "weight" edge attributes.

    Returns
    -------
    numpy.ndarray
        The matrix Q = M_T L_T M_T^T + L_g.
    """
    M, ordered_edges = cycle_edge_incidence_matrix(G)
    print(M)
    n = G.number_of_nodes()
    tree_size = n - 1

    L_T, L_g = length_matrices(G, ordered_edges, tree_size)

    M_T = M[:, :tree_size]

    Q = M_T @ L_T @ M_T.T + L_g
    return Q

def pad_Q_matrices(Q1: np.ndarray, Q2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Zero-pad the smaller Q matrix so both have the same shape.

    Parameters
    ----------
    Q1, Q2 : numpy.ndarray
        Square matrices to be padded to a common size.

    Returns
    -------
    Q1_padded, Q2_padded : numpy.ndarray
        Padded copies of Q1 and Q2 with matching dimensions.
    """
    n1 = Q1.shape[0]
    n2 = Q2.shape[0]
    g0 = max(n1, n2)

    Q1_pad = np.zeros((g0, g0))
    Q2_pad = np.zeros((g0, g0))

    Q1_pad[:n1, :n1] = Q1
    Q2_pad[:n2, :n2] = Q2

    return Q1_pad, Q2_pad

def ktte(Q1: np.ndarray, Q2: np.ndarray, gamma: float = 1.0) -> float:
    """Tropical Torelli–Euclidean kernel.

    Parameters
    ----------
    Q1, Q2 : numpy.ndarray
        SPD / PSD matrices (square).
    gamma : float, optional
        Kernel bandwidth.

    Returns
    -------
    float
        Kernel value.
    """
    Q1p, Q2p = pad_Q_matrices(Q1, Q2)
    diff = Q1p - Q2p
    frob_sq = np.linalg.norm(diff, ord="fro") ** 2
    return np.exp(-gamma * frob_sq)

def kttw(Q1: np.ndarray, Q2: np.ndarray, gamma: float = 1.0) -> float:
    """Tropical Torelli–Wasserstein kernel using matrix square roots.

    Parameters
    ----------
    Q1, Q2 : numpy.ndarray
        PSD matrices of the same shape.
    gamma : float, optional
        Kernel bandwidth.

    Returns
    -------
    float
        k_TTW(Q1, Q2).
    """
    Q1, Q2 = pad_Q_matrices(Q1, Q2)

    sqrtQ1 = sqrtm(Q1)
    middle = sqrtQ1 @ Q2 @ sqrtQ1
    sqrt_middle = sqrtm(middle)

    W2 = (
        np.trace(Q1)
        + np.trace(Q2)
        - 2.0 * np.trace(sqrt_middle)
    )

    # numerical safety: imaginary parts should be ~0
    W2 = np.real(W2)

    return np.exp(-gamma * W2)
