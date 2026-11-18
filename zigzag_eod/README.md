# Zigzag persistent homology of Medellín's 2017 OD survey

Zigzag persistence on the AMVA 2017 origin–destination survey read as a temporal
graph: SIT zones are nodes, a zone pair is an edge while a trip occupies it.
Follows Myers et al., *Temporal network analysis using zigzag persistence*, using
Dionysus 2 as the persistence engine — the same engine Myers used.

## Setup

```bash
pip install -r requirements.txt
```

That is all. Dionysus 2 exposes `zigzag_homology_persistence` directly, so no
compiled wrapper is needed. Building it from source needs CMake and Boost headers
(`conda install cmake boost-cpp`) plus a C++11 compiler; the PyPI package builds
these automatically on install.

The notebook validates the engine against two published reference results before
touching real data — GUDHI's documented zigzag example and the worked toy barcode
in Myers' Figure 2. If either fails, nothing downstream is trustworthy.

Dionysus stores filtration values as **float32**, so reference comparisons use a
1e-4 tolerance. Our own filtration values are multiples of 7.5 minutes, which
float32 represents exactly, so no precision is lost on real data.

## Running

All paths in the notebook are **relative to this directory**, so run it from here:

```bash
jupyter lab eod_zigzag_persistence.ipynb     # interactive
python run_notebook.py                       # headless, in-process
```

`run_notebook.py` exists because `nbconvert --execute` needs to bind a socket for
its ZMQ kernel, which some sandboxes forbid; it execs the cells in a shared
namespace and writes outputs back into the `.ipynb`. Use `jupyter lab` normally.

`build_notebook.py` regenerates the notebook from source. **Edit the builder, not
the notebook**, if you want changes to survive — re-running it overwrites
`eod_zigzag_persistence.ipynb`.

Tables and figures are written to `outputs/`.

## Parameters

All in section 1 of the notebook; everything downstream is a function of them.

| Parameter | Default | Meaning |
|---|---|---|
| `WINDOW_MIN` | 60 | snapshot width, minutes |
| `OVERLAP` | 0.5 | overlap fraction between consecutive windows |
| `RIPS_R` | 1 | Vietoris–Rips scale on the snapshot shortest-path metric |
| `MAX_DIM` | 2 | complex truncation dimension (H0 and H1 interpretable) |
| `MIN_TRIPS` | 100 | **edge support threshold** — see below |
| `TRIP_COUNT` | `"expanded"` | how `MIN_TRIPS` is counted |

At those defaults the notebook produces **16,881 bars** (H0 160, H1 10,569, plus
6,152 uninterpretable dimension-2 bars) from 48 snapshots averaging 1,402 edges.
Section 6 asserts that count, so an accidental change upstream fails loudly
instead of quietly producing different findings.

### `MIN_TRIPS` — the support threshold

An edge `{A,B}` enters a snapshot only if it carries at least `MIN_TRIPS` trips in
that window. Counted two ways:

* **`"surveyed"`** — raw survey rows. Sparse: the median zone pair in a 60-min
  window carries exactly *one* surveyed trip, so `MIN_TRIPS = 3` already discards
  ~91% of edges. Useful range 1–5.
* **`"expanded"`** — sum of the expansion factor `fe_estratos_sin a pie 1 y 2`,
  i.e. population-scaled trips. A row expands to ~73 trips at the median, so
  thresholds live on a different scale: 100 keeps ~44% of edges, 500 keeps ~1%.

The two countings are not interchangeable; a threshold is comparable only against
another run with the same `TRIP_COUNT`. Bar counts are likewise comparable only
within one `WINDOW_MIN` — see the sensitivity sweep in the notebook's last section.

### The window boundary rule

An edge is in window `[lo, hi)` when `t_o < hi and t_d >= lo` — the closed
activity interval `[t_o, t_d]` meeting the half-open window. **Do not weaken the
`>=` to `>`.** Reported times heap on `:00`/`:30` and windows start on `:00`/`:30`,
so excluding trips that arrive exactly on a boundary drops 41,147 trip/window
incidences, thins the mean snapshot from 1,402 to 1,157 edges, and moves the
barcode from 16,881 bars to 14,947. Both are internally consistent; they are
different definitions of the temporal graph, and this project uses the closed one.
If you find a 14,947 figure in older notes, it came from the half-open variant.

The sweep result matters for interpretation, and it cuts both ways. As the expanded
threshold rises, H1 bars collapse (14,361 → 13,629 → 10,569 → 32 at
`MIN_TRIPS` = 1, 50, 100, 500) while H0 bars *rise* (118 → 124 → 160 → 219):
demanding stronger links fragments the graph and starves it of cycles at the same
time. The day-spanning H0 class is essential at thresholds up to 100 but **dies at
21:30 when `MIN_TRIPS = 500`**, so "connected all day" is a statement about ordinary
link strengths, not an invariant. See `outputs/fig_trip_threshold.png` and
`outputs/trip_threshold_sensitivity.csv`.

## Layout

```
data/                    trip workbook + zone shapefile (not tracked)
outputs/                 generated tables and figures
build_notebook.py        regenerates the notebook
run_notebook.py          headless executor
eod_zigzag_persistence.ipynb
concluding_remarks.md    findings, limitations, next steps
```

## Reading the output

Start with `concluding_remarks.md`. Two caveats that constrain interpretation:

* **Dimension-2 bars are not topology.** At `MAX_DIM = 2` there are no
  3-simplices to bound 2-cycles, so those bars count unfilled tetrahedra of the
  truncation. They carry `interpretable = False` in `zigzag_bars.csv`.
* **H1 localization is indicative.** Dionysus' zigzag returns bars without
  representative cycles, so the map shows chordless cycles from the bar's birth
  snapshot — guaranteed to carry genuine H1 classes at `RIPS_R = 1`, but not
  guaranteed to be the class the algorithm paired. H0 localization is exact.
