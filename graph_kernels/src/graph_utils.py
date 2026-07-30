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
        edge_color="gray",
    )

    edge_labels = nx.get_edge_attributes(G, "weight")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10)
    plt.show()

def plot_weighted_graphs(graphs: List[nx.Graph], cols: int = 4, seed: int = 42) -> None:
    """Plot a list of undirected weighted NetworkX graphs in a grid layout.

    Parameters
    ----------
    graphs : list of networkx.Graph
        Graphs with numeric "weight" edge attributes.
    cols : int, optional
        Number of columns in the subplot grid.
    seed : int, optional
        Seed for layout reproducibility.
    """
    n = len(graphs)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 3))
    axes = np.array(axes).flatten()

    for i, G in enumerate(graphs):
        ax = axes[i]
        pos = nx.spring_layout(G, seed=seed)
        nx.draw(
            G,
            pos,
            ax=ax,
            with_labels=False,
            node_size=20,
            node_color="lightblue",
            edge_color="gray",
        )
        title = G.graph.get("name", f"G{i + 1}")
        ax.set_title(title, fontsize=8)

    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
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
    """Length of the shortest tree edge at endpoint excluding the edge to excluded_neighbor.

    When degree > 2 in the tree we use the shortest candidate to keep orientation
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
    """Return a canonically oriented edge.

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
    Canonical orientation (Appendix A.2 of the paper):
    1) Tree edge e=(u,v): omega(x) = length of the other incident tree edge at x
       (0 if leaf). Orient low-omega -> high-omega.
    2) Non-tree edge e=(u,v): omega(x) = length of the first tree edge on the
       unique path from x toward the other endpoint in T. Orient low -> high.
    Ties fall back to deterministic node ordering.
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
        Oriented edges e1,...,em with tree edges first, then non-tree edges
        sorted ascending by weight (Algorithm 1, step 4).
    tree_edges : list of tuple
        First n-1 canonically oriented edges.
    non_tree_edges : list of tuple
        Remaining g edges as (u, v, weight), sorted ascending by weight.
    """
    T = nx.minimum_spanning_tree(G, weight="weight")

    tree_edges = [orient_edge(u, v, T=T, tree_edge=True) for u, v in T.edges()]
    # Tree-edge ordering within T does not affect Q (proven in Appendix A.2),
    # but we sort for reproducibility.
    tree_edges.sort(
        key=lambda e: (float(T[e[0]][e[1]]["weight"]), repr(e[0]), repr(e[1]))
    )

    tree_edge_keys = {_undirected_edge_key(u, v) for u, v in T.edges()}

    non_tree_edges: List[Tuple[Hashable, Hashable, float]] = []
    for u, v, d in G.edges(data=True):
        if _undirected_edge_key(u, v) in tree_edge_keys:
            continue
        ou, ov = orient_edge(u, v, T=T, tree_edge=False)
        non_tree_edges.append((ou, ov, float(d["weight"])))

    # Non-tree edges MUST be sorted ascending by length (Algorithm 1, step 4).
    non_tree_edges.sort(key=lambda x: (x[2], repr(x[0]), repr(x[1])))

    ordered_edges = tree_edges + [(u, v) for u, v, _ in non_tree_edges]

    return T, ordered_edges, tree_edges, non_tree_edges


def fundamental_cycle(
    T: nx.Graph,
    non_tree_edge: Tuple[Hashable, Hashable],
) -> List[Tuple[Hashable, Hashable, float]]:
    """Return signed oriented edges for the fundamental cycle of a non-tree edge.

    The fundamental cycle σ induced by non-tree edge (u, v) is:
        σ = (u→v)  +  tree-path from v back to u

    Each returned entry is (src, dst, sign) where sign ∈ {+1, -1}:
      +1  if the edge's canonical orientation agrees with the cycle traversal,
      -1  otherwise.

    Parameters
    ----------
    T : networkx.Graph
        Minimum spanning tree.
    non_tree_edge : tuple (u, v)
        Canonically oriented non-tree edge that closes the cycle.

    Returns
    -------
    list of (src, dst, sign)
    """
    u, v = non_tree_edge

    # Tree path goes from v back to u (completing the cycle after u→v).
    path = nx.shortest_path(T, source=v, target=u)

    signed_edges: List[Tuple[Hashable, Hashable, float]] = []

    # --- tree-path edges (v → ... → u) ---
    for i in range(len(path) - 1):
        a, b = path[i], path[i + 1]
        canonical = orient_edge(a, b, T=T, tree_edge=True)
        # Traversal direction is a→b; sign = +1 if canonical agrees, else -1.
        sign = 1.0 if canonical == (a, b) else -1.0
        signed_edges.append((canonical[0], canonical[1], sign))

    # --- the non-tree edge itself (u → v), always +1 in its own cycle ---
    signed_edges.append((u, v, 1.0))

    return signed_edges


def cycle_edge_incidence_matrix(
    G: nx.Graph,
) -> Tuple[np.ndarray, List[Tuple[Hashable, Hashable]]]:
    """Build the signed cycle-edge incidence matrix M.

    M[i, j] = +1  if the j-th edge appears in the i-th fundamental cycle
                   with the same orientation as ordered_edges[j],
             = -1  if it appears with opposite orientation,
             =  0  if it does not appear.

    The block structure is M = [M_T | I_g] after the canonical edge ordering
    (Algorithm 1 / eq. 1 of the paper).

    Parameters
    ----------
    G : networkx.Graph
        Connected undirected weighted graph with "weight" edge attributes.

    Returns
    -------
    M : numpy.ndarray, shape (g, m)
        Signed cycle-edge incidence matrix.
    ordered_edges : list of tuple
        Oriented edges e1,...,em consistent with the columns of M.
    """
    T, ordered_edges, tree_edges, non_tree_edges = compute_mst_and_edge_order(G)

    g = len(non_tree_edges)
    m = len(ordered_edges)

    # Build forward lookup: canonical oriented edge -> column index.
    edge_index: dict[Tuple[Hashable, Hashable], int] = {
        e: i for i, e in enumerate(ordered_edges)
    }

    M = np.zeros((g, m))

    for i, (u, v, _) in enumerate(non_tree_edges):
        signed_cycle = fundamental_cycle(T, (u, v))

        for (src, dst, sign) in signed_cycle:
            canonical = (src, dst)
            if canonical not in edge_index:
                raise KeyError(
                    f"Edge {canonical} not found in edge_index. "
                    "This should not happen; check orient_edge consistency."
                )
            j = edge_index[canonical]
            M[i, j] = sign

    return M, ordered_edges


def length_matrices(
    G: nx.Graph,
    ordered_edges: Sequence[Tuple[Hashable, Hashable]],
    tree_size: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute diagonal length matrices L_T and L_g.

    Parameters
    ----------
    G : networkx.Graph
        Undirected weighted graph with "weight" edge attributes.
    ordered_edges : sequence of tuple
        Oriented edges e1,...,em (tree edges first, then non-tree).
    tree_size : int
        Number of tree edges (n - 1 for a connected graph).

    Returns
    -------
    L_T : numpy.ndarray
        Diagonal matrix of tree-edge lengths, shape (tree_size, tree_size).
    L_g : numpy.ndarray
        Diagonal matrix of non-tree-edge lengths, shape (g, g).
    """
    lengths = [G[u][v]["weight"] for u, v in ordered_edges]
    L = np.diag(lengths)
    L_T = L[:tree_size, :tree_size]
    L_g = L[tree_size:, tree_size:]
    return L_T, L_g


def compute_Q(G: nx.Graph, alpha=1.0) -> np.ndarray:
    """Compute the tropical Torelli matrix Q for a connected weighted graph.

    Implements Algorithm 1 from the paper:
        Q = M_T L_T M_T^T + L_g

    Parameters
    ----------
    G : networkx.Graph
        Connected undirected weighted graph with "weight" edge attributes.

    Returns
    -------
    numpy.ndarray
        SPD matrix Q of shape (g, g), where g = m - n + 1 is the graph genus.
        Returns an empty (0, 0) array for trees (g = 0).

    Raises
    ------
    ValueError
        If G is not connected. Use component-wise computation for disconnected
        graphs (the paper handles this via a block-diagonal construction).
    """
    if not nx.is_connected(G):
        raise ValueError(
            "compute_Q requires a connected graph. "
            "For disconnected graphs, call compute_Q on each connected component "
            "and form a block-diagonal matrix."
        )

    n = G.number_of_nodes()
    m = G.number_of_edges()
    g = m - n + 1

    if g == 0:
        # Tree: tropical Torelli map is trivial; all trees map to the zero form.
        return np.zeros((0, 0))

    M, ordered_edges = cycle_edge_incidence_matrix(G)

    tree_size = n - 1
    L_T, L_g = length_matrices(G, ordered_edges, tree_size)

    M_T = M[:, :tree_size]
    Q = M_T @ L_T @ M_T.T + alpha* L_g
    return Q


def compute_Q_disconnected(G: nx.Graph) -> np.ndarray:
    """Compute the block-diagonal tropical Torelli matrix for a (possibly disconnected) graph.

    Parameters
    ----------
    G : networkx.Graph
        Undirected weighted graph with "weight" edge attributes.

    Returns
    -------
    numpy.ndarray
        Block-diagonal SPD matrix Q.
    """
    blocks = []
    for component_nodes in nx.connected_components(G):
        subgraph = G.subgraph(component_nodes).copy()
        Q_block = compute_Q(subgraph)
        if Q_block.size > 0:
            blocks.append(Q_block)

    if not blocks:
        return np.zeros((0, 0))

    # Build block-diagonal matrix
    total = sum(b.shape[0] for b in blocks)
    Q_full = np.zeros((total, total))
    offset = 0
    for b in blocks:
        s = b.shape[0]
        Q_full[offset:offset + s, offset:offset + s] = b
        offset += s

    return Q_full


# ---------------------------------------------------------------------------
# Kernels
# ---------------------------------------------------------------------------

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
    """Tropical Torelli–Euclidean (TTE) kernel.

        k_TTE(G1, G2) = exp(-gamma * ||Q1 - Q2||_F^2)

    Parameters
    ----------
    Q1, Q2 : numpy.ndarray
        SPD / PSD matrices (square).
    gamma : float, optional
        Kernel bandwidth.

    Returns
    -------
    float
        Kernel value in (0, 1].
    """
    Q1p, Q2p = pad_Q_matrices(Q1, Q2)
    frob_sq = np.linalg.norm(Q1p - Q2p, ord="fro") ** 2
    return float(np.exp(-gamma * frob_sq))


def kttw(Q1: np.ndarray, Q2: np.ndarray, gamma: float = 1.0) -> float:
    """Tropical Torelli–Wasserstein (TTW) kernel.

    Uses the Bures–Wasserstein distance (eq. 3 of the paper):
        W^2(Q1, Q2) = Tr(Q1) + Tr(Q2) - 2 Tr(sqrt(sqrt(Q1) Q2 sqrt(Q1)))

        k_TTW(G1, G2) = exp(-gamma * W^2(Q1, Q2))

    Parameters
    ----------
    Q1, Q2 : numpy.ndarray
        PSD matrices (square).
    gamma : float, optional
        Kernel bandwidth.

    Returns
    -------
    float
        Kernel value in (0, 1].
    """
    Q1p, Q2p = pad_Q_matrices(Q1, Q2)

    sqrtQ1 = sqrtm(Q1p)
    middle = sqrtQ1 @ Q2p @ sqrtQ1
    sqrt_middle = sqrtm(middle)

    W2 = np.trace(Q1p) + np.trace(Q2p) - 2.0 * np.trace(sqrt_middle)
    W2 = float(np.real(W2))  # imaginary parts are numerical noise; discard

    return float(np.exp(-gamma * W2))


# ---------------------------------------------------------------------------
# Quick self-test: reproduce the paper's Appendix A.3 example
# ---------------------------------------------------------------------------

def _build_paper_example() -> nx.Graph:
    """Construct the toy graph from Appendix A.3.

    6 nodes, 8 edges with ℓ(eᵢ) = i.
    Expected Q = [[6, 0, 2], [0, 18, 6], [2, 6, 20]].
    """
    G = nx.Graph()
    # Nodes 1-6 (matching the figure labels)
    edges = [
        (1, 2, 1),  # e1
        (2, 3, 2),  # e2
        (1, 3, 3),  # e3  <- non-tree (closes cycle with e1, e2)
        (3, 4, 4),  # e4
        (4, 5, 5),  # e5
        (4, 6, 6),  # e6
        (5, 6, 7),  # e7  <- non-tree (closes cycle with e5, e6)
        (2, 6, 8),  # e8  <- non-tree (closes cycle with e2, e4, e6)
    ]
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)
    return G


if __name__ == "__main__":
    G = _build_paper_example()
    Q = compute_Q(G)
    print("Computed Q:")
    print(Q)
    print()
    expected = np.array([[6, 0, 2], [0, 18, 6], [2, 6, 20]], dtype=float)
    print("Expected Q (Appendix A.3):")
    print(expected)
    print()
    print("Match:", np.allclose(Q, expected))



