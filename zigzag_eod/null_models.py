"""Null models for the zigzag barcode: is there signal beyond degree and volume?

Runs the notebook's own pipeline up to the snapshot complexes `Ks`, then compares the
observed barcode against two nulls at the default configuration:

  order  -- permute the 48 snapshots. Each snapshot's topology is preserved EXACTLY
            (identical per-window betti numbers); only the temporal ordering changes.
            Tests whether smooth temporal evolution -- the thing zigzag is for -- is
            real, or whether the bars would look the same in any order.
  rewire -- per window, degree-preserving double-edge swaps. Preserves each window's
            edge count and every zone's degree; destroys which specific zone pairs are
            joined. Tests whether the cycle structure is spatial or just degree.

Also reports whether the number of H1 bars alive per window is predictable from that
window's static beta_1 -- if it is, the zigzag adds nothing over per-snapshot homology.

Writes outputs/null_models.csv and outputs/null_summary.json. Reproducible via SEED.
"""
import json
import time

import matplotlib
matplotlib.use("Agg")
import networkx as nx
import nbformat as nbf
import numpy as np
import pandas as pd

SEED = 20170101
N_REP = 20

ns = {"__name__": "__main__"}
nb = nbf.read("eod_zigzag_persistence.ipynb", as_version=4)
for cell in nb.cells:
    if cell.cell_type != "code":
        continue
    exec(cell.source, ns)
    # Stop once the interval encoding exists; the observed barcode is recomputed below.
    if "intervals = zigzag_intervals(Ks, centers)" in cell.source:
        break

Ks              = ns["Ks"]
centers         = ns["centers"]
starts          = ns["starts"]
zigzag_intervals = ns["zigzag_intervals"]
zigzag_bars     = ns["zigzag_bars"]
rips_complex    = ns["rips_complex"]
window_edges    = ns["window_edges"]
RIPS_R, MAX_DIM, WINDOW_MIN = ns["RIPS_R"], ns["MAX_DIM"], ns["WINDOW_MIN"]
DEATH_CAP = starts[-1] + WINDOW_MIN


def barcode_stats(Klist, label):
    """Zigzag barcode summary for a snapshot sequence. H0/H1 only -- dim 2 is truncation."""
    t0 = time.time()
    iv = zigzag_intervals(Klist, centers)
    bars = zigzag_bars(iv)
    rows = []
    for dim, b, d in bars:
        open_ = not np.isfinite(d)
        rows.append((dim, b, DEATH_CAP if open_ else d, open_))
    df = pd.DataFrame(rows, columns=["dim", "birth", "death", "open"])
    df["life_h"] = (df.death - df.birth) / 60.0

    out = {"label": label, "secs": round(time.time() - t0, 1),
           "n_simplices": len(iv)}
    for dim in (0, 1):
        s = df[df.dim == dim]
        out[f"h{dim}_bars"]     = int(len(s))
        out[f"h{dim}_open"]     = int(s.open.sum())
        out[f"h{dim}_longest"]  = round(float(s.life_h.max()), 2) if len(s) else 0.0
        out[f"h{dim}_median"]   = round(float(s.life_h.median()), 2) if len(s) else 0.0
        out[f"h{dim}_over_3h"]  = int((s.life_h > 3).sum())
    return out, df


print("=== observed ===", flush=True)
obs, obs_df = barcode_stats(Ks, "observed")
print(json.dumps(obs), flush=True)

# Does the zigzag add anything over per-snapshot beta_1? Compare the number of H1 bars
# alive at each window centre against that window's static beta_1.
h1 = obs_df[obs_df.dim == 1]
alive = [int(((h1.birth <= c) & (h1.death >= c)).sum()) for c in centers]
b1 = [ns["static_betti"](K)[1] for K in Ks] if "static_betti" in ns else None
if b1 is None:
    import gudhi
    def _b1(K):
        st = gudhi.SimplexTree()
        for s in K:
            st.insert(list(s), 0.0)
        st.compute_persistence(persistence_dim_max=True)
        return st.betti_numbers()[1] if len(st.betti_numbers()) > 1 else 0
    b1 = [_b1(K) for K in Ks]
corr = float(np.corrcoef(alive, b1)[0, 1])
print(f"H1-alive vs static beta_1: pearson r = {corr:.4f}", flush=True)

records = [obs]
rng = np.random.default_rng(SEED)

print(f"\n=== null: snapshot order shuffle ({N_REP} reps) ===", flush=True)
for rep in range(N_REP):
    perm = rng.permutation(len(Ks))
    st, _ = barcode_stats([Ks[i] for i in perm], "order")
    st["rep"] = rep
    records.append(st)
    print(f"  rep {rep:2d}: H0 {st['h0_bars']:4d} (open {st['h0_open']:2d}, "
          f"longest {st['h0_longest']:5.2f}h) | H1 {st['h1_bars']:6d} "
          f"(longest {st['h1_longest']:5.2f}h, >3h {st['h1_over_3h']:4d})", flush=True)

print(f"\n=== null: per-window degree-preserving rewire ({N_REP} reps) ===", flush=True)
for rep in range(N_REP):
    Kr = []
    for s in starts:
        E = window_edges(s, s + WINDOW_MIN)
        G = nx.Graph(); G.add_edges_from(E)
        if G.number_of_edges() >= 2:
            # nx needs >=4 nodes and enough edges for swaps to be possible
            try:
                nx.double_edge_swap(G, nswap=10 * G.number_of_edges(),
                                    max_tries=200 * G.number_of_edges(),
                                    seed=int(rng.integers(1 << 31)))
            except (nx.NetworkXError, nx.NetworkXAlgorithmError):
                pass  # too small/constrained to rewire; keep as observed
        Kr.append(rips_complex(set(map(lambda e: tuple(sorted(e)), G.edges)),
                               RIPS_R, MAX_DIM))
    st, _ = barcode_stats(Kr, "rewire")
    st["rep"] = rep
    records.append(st)
    print(f"  rep {rep:2d}: H0 {st['h0_bars']:4d} (open {st['h0_open']:2d}, "
          f"longest {st['h0_longest']:5.2f}h) | H1 {st['h1_bars']:6d} "
          f"(longest {st['h1_longest']:5.2f}h, >3h {st['h1_over_3h']:4d})", flush=True)

res = pd.DataFrame(records)
res.to_csv("outputs/null_models.csv", index=False)

summary = {"seed": SEED, "n_rep": N_REP, "h1_alive_vs_static_beta1_r": round(corr, 4),
           "observed": obs}
for label in ("order", "rewire"):
    sub = res[res.label == label]
    summary[label] = {}
    for col in ("h0_bars", "h0_open", "h0_longest", "h1_bars", "h1_longest",
                "h1_median", "h1_over_3h", "n_simplices"):
        summary[label][col] = {"mean": round(float(sub[col].mean()), 2),
                               "sd": round(float(sub[col].std()), 2),
                               "min": float(sub[col].min()),
                               "max": float(sub[col].max())}
        # one-sided empirical p: how often does the null match or exceed observed?
        summary[label][col]["p_ge_obs"] = round(
            float((sub[col] >= obs[col]).mean()), 3)
json.dump(summary, open("outputs/null_summary.json", "w"), indent=1)

print("\n=== summary ===")
print(json.dumps(summary, indent=1))
