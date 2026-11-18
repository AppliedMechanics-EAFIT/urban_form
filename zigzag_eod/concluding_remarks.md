# Zigzag persistent homology of Medellín's 2017 OD survey — concluding remarks

Default configuration: `WINDOW_MIN = 60`, `OVERLAP = 0.5` (48 windows), `RIPS_R = 1`,
`MAX_DIM = 2`, no trip filtering. All numbers below come from
`eod_zigzag_persistence.ipynb` at those settings.

## What was built

The survey's 87,613 usable trips become a temporal graph: SIT zones are nodes, and
each interzonal trip makes its zone pair active on the closed interval `[t_o, t_d]`
— the departure-to-arrival reading you chose, so an edge exists for as long as a
trip occupies it. 15,683 intrazonal trips carry no edge, leaving 71,930 activity
intervals over 544 zones and 17,931 distinct undirected pairs. Sliding a 60-minute
window in 30-minute steps gives 48 snapshots; each is replaced by the Vietoris–Rips
complex of its own shortest-path metric at scale 1, truncated at dimension 2.
Adjacent snapshots are interleaved with their unions to form a valid zigzag, and
the sequence is streamed through GUDHI's `Filtered_zigzag_persistence`.

GUDHI's zigzag module is C++-only — no Python binding ships in the wheel — so
`gudhi_zigzag_ext/` is a thin pybind11 wrapper over that class, compiling against
the installed `gudhi` 3.13 headers. It is regression-tested against GUDHI's own
documented example and against the toy barcode Myers states in his Figure 2; both
reproduce exactly. Your GUDHI requirement is met with no algorithmic
reimplementation on my part.

## Findings

**1. Medellín's mobility network has one topological regime, not two.**
The essential H0 bar runs the whole day. Every other H0 bar — 115 of them — is
born before 05:00 or after 22:00, and the longest non-essential one lasts only
2.5 hours. Between roughly 05:00 and 22:00 the network is a single connected
component in 35 of 48 windows. Fragmentation is a night phenomenon exclusively;
there is no midday disconnection anywhere in the metropolitan area.

**2. Fragmentation is broad but shallow.**
238 of 544 zones detach from the core in at least one window, which sounds like
severe peripheral isolation until you look at the duration: only 66 zones detach
in three or more windows, and the maximum for any single zone is 7 of 48. The
zigzag makes this distinction available in a way per-snapshot component counts do
not — a count of 35 components at 03:30 says nothing about whether those are the
same 35 zones as at 03:00. The barcode says they largely are not.

**3. Cycles are abundant but almost all transient.**
14,417 H1 bars, median lifetime 1.0 hour — i.e. a typical cycle of mutual zone
exchange survives one window step and dies. Only 228 bars (1.6%) exceed 3 hours
and only 3 exceed 6. The single longest, 05:45→16:00 (10.25 h), spans the morning
peak through the afternoon. So the mobility network continuously generates and
destroys circulation loops rather than maintaining a stable cyclic backbone. This
is the finding I would flag as most worth following up, because it is exactly what
the static β₁ curve conceals: β₁ sits between 500 and 1,200 for the entire
daytime and peaks at 1,194 at 20:00, which reads as a persistently loopy network.
The zigzag shows that near-constant β₁ is a turnover equilibrium, not persistence.

**4. Part of the barcode is an artifact of the window width.**
Halving the window from 60 to 30 minutes raises H1 bar count from 14,417 to
20,658; doubling it to 120 drops it to 7,968. H0 moves the same way (116 → 180 →
59). Longest-lived H1 persistence, by contrast, *increases* with window width
(10.25 h → 11.5 h at 120 min), because wider windows keep edges alive across gaps
that would otherwise break a cycle. Bar counts are therefore not a comparable
quantity across configurations — only within one. Any cross-city or cross-survey
comparison would need the window fixed by an external argument, not by convenience.

**5. Most of the H1 barcode is survey noise — and the H0 signal is not.**
The `MIN_TRIPS` sweep (Figure 7, `trip_threshold_sensitivity.csv`) settles the
question finding 3 raised. Requiring two surveyed trips instead of one to call a
zone pair connected removes 78% of edges (3,368 → 726 per window on average) and
takes H1 with it: 14,417 bars → 6,175, longest lifetime 10.25 h → 5.75 h. At
`MIN_TRIPS = 3` only 1,681 bars remain with nothing above 3.25 h; at 5, just 133
bars and 2.5 h. Cycle count and longest lifetime fall *together*, in near-lockstep
— the signature of a barcode built from single-report pairs rather than one with a
robust core. **No H1 feature in this survey survives a support threshold, so none
should be interpreted as a durable circulation loop.**

H0 behaves in the opposite, informative way. Bar count *rises* (116 → 237 → 372)
as thresholding disconnects the graph and creates genuine new components, while
the longest bar stays essential — 24 h at `MIN_TRIPS = 1`, 23.25 h at 2, 22.75 h
at 3. The daytime single-component finding is therefore not an artifact of
single-trip edges: the metropolitan core stays connected even when four fifths of
its edges are discarded. Finding 1 is robust; finding 3 was correctly flagged as
the fragile one.

Note that thresholds on `TRIP_COUNT = "expanded"` are not comparable to these:
one survey row expands to ~73 trips at the median, so an equivalent sweep lives
on the scale of tens to hundreds.

## Limitations you should know about before building on this

- **H2 is not interpretable here.** The complex is truncated at dimension 2, so
  there are no 3-simplices to bound 2-cycles; the 81,032 dimension-2 bars count
  unfilled tetrahedra of the truncation, not voids in the mobility network. They
  are written to `zigzag_bars.csv` with `interpretable = False`. Honest H2 needs
  `MAX_DIM = 3`.
- **Cycle localization is indicative, not exact.** GUDHI's zigzag returns bars
  only — no representative cycles — and Myers flags generator extraction as
  unimplemented future work in his own pipeline too. H0 localization *is* exact
  (components are directly enumerable, hence
  `h0_fragmentation_by_zone.csv`). For H1 the notebook reports chordless cycles of
  length ≥ 4 in the bar's birth snapshot, which are guaranteed to carry genuine H1
  classes at `RIPS_R = 1` but are not guaranteed to be the specific class GUDHI
  paired. The strict version — cycles present in *every* window the bar spans —
  returned nothing for the longest bar, so treat Figure 5c's zones (162, 356, 360,
  362, 375) as *a* carrier in the birth window, not *the* generator.
- **Departure times are heaped.** 74% of reported departures fall on :00 or :30.
  The duration-interval edge model absorbs most of this — an edge lives across
  window boundaries rather than appearing only at its departure instant — which is
  the main reason I recommended it over the instantaneous reading.
- **No expansion factors.** Per your instruction, trips are unweighted: the
  topology reflects surveyed trips, not the ~2.7M population-scaled trips
  `fe_estratos_sin a pie 1 y 2` encodes. Since the complex only asks *whether* a
  zone pair is active, weighting would change nothing unless it were used as a
  threshold — which is the filtering discussion you deferred.

## Where I would go next

The filtering question is now answered for the H1 side, and the answer is
negative: the prediction that a modest threshold would leave "a much smaller
barcode whose bars mean something" was wrong. Nothing survives — count and
lifetime fall together at every threshold tested. The productive move is
therefore *not* a finer τ sweep on this construction but a change to the
construction itself. Two candidates:

* **Aggregate before thresholding.** The sparsity is per-pair-per-window; pooling
  trips over a coarser spatial partition (comunas rather than SIT zones) or a
  coarser time base would raise support without discarding information, and would
  test whether the transient cycles are a zone-resolution artifact.
* **Threshold on `"expanded"` counts instead.** Population-scaled trips are much
  less sparse (98% of pairs clear τ = 1, 44% clear 100), so an expanded-count
  sweep probes a genuinely different question — not "was this pair reported twice"
  but "does this pair carry meaningful volume". `TRIP_COUNT = "expanded"` is
  already wired; only the sweep values need changing.

The H0 results need no such rescue and are worth writing up as they stand.

Mode-stratified zigzags (metro vs. bus vs. car, using `modo principal`) would also
be a direct extension, and one where cross-barcode comparison is legitimate because
the window is shared.
