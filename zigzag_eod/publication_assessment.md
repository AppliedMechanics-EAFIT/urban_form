# Is there a paper here? — assessment against literature and null models

Written after (a) locating the nearest prior work and (b) running two null models at the
default configuration (`WINDOW_MIN=60`, `OVERLAP=0.5`, `RIPS_R=1`, `MAX_DIM=2`,
`MIN_TRIPS=100`, `TRIP_COUNT="expanded"`). Numbers from `outputs/null_summary.json`
and `outputs/null_models.csv`, 20 replicates each, `SEED = 20170101`.

## 1. The method is already published — twice over

Myers, Muñoz, Khasawneh & Munch, *Temporal network analysis using zigzag persistence*,
EPJ Data Science 12:6 (2023), doi:10.1140/epjds/s13688-023-00379-5. This is the paper
the notebook already cites as its engine reference, and it does what this notebook
does: sliding windows over a temporal graph, Vietoris–Rips complex per snapshot, unions
as the zigzag bridge, Dionysus for the barcode, `w = 20 min` with 50% overlap.

Their headline application is **the Great Britain transportation network** — rail,
coach and air. So "zigzag persistence applied to a transportation/commuting network"
is not available as a contribution. There is no methods novelty here, and the framing
"first application of zigzag persistence to urban mobility" would not survive review.

What their paper does *not* do, confirmed by full-text search (0 hits for "null model",
0 for "randomi"): any null-model validation. Their own stated future work names
"an analysis on deciding an optimal window size and overlap" and "a method to
incorporate edge weight and directionality" — both of which this notebook has,
in the sweeps and in `MIN_TRIPS`.

One substantive difference is worth keeping in view. Myers et al. use **supply-side**
data: nodes are scheduled destinations, edges are timetabled connections, and the
window width is chosen to *guarantee* connectivity (20 min against a 7-min mean wait).
This notebook uses **demand-side** data: an OD household survey, where edges are
revealed trips carrying an expansion-weighted volume, and where the threshold
deliberately *removes* edges. Nobody has run zigzag persistence on a household travel
survey. That is a real gap — but a gap in application, not in method, and gaps of that
kind carry a paper only when the application returns something the method's own authors
could not have predicted.

## 2. What the nulls say

Two nulls, both at the default configuration:

- **order** — permute the 48 snapshots. Each snapshot's topology is preserved *exactly*
  (identical per-window Betti numbers, identical simplex count: 23,461 in every replicate).
  Isolates whether the *temporal ordering* carries information.
- **rewire** — per window, degree-preserving double-edge swaps. Preserves every window's
  edge count and every zone's degree; destroys which specific zone pairs are joined.
  Isolates whether the structure is *spatial* or merely *degree*.

| statistic | observed | order null (mean ± sd) | rewire null (mean ± sd) | verdict |
|---|---|---|---|---|
| H1 bars living > 3 h | **209** | 0.3 ± 0.7 (max 3) | 0.0 ± 0.0 | **real** |
| longest H1 bar (h) | **5.75** | 2.55 ± 0.41 | 1.68 ± 0.24 | **real** |
| median H1 lifetime (h) | **1.25** | 0.75 ± 0.0 | 0.75 ± 0.0 | **real** |
| total H1 bars | **10,569** | 26,573 ± 992 | 37,809 ± 141 | **real** (observed far *fewer*) |
| longest H0 bar (h) | 23.75 | 17.71 ± 3.64 | **23.75 ± 0.00** | **artifact** |
| essential (open) H0 bars | 32 | 4.7 ± 7.8 | **32.25 ± 0.55** | **artifact** |
| total H0 bars | 160 | 234 ± 7 | 216 ± 10 | real (observed *fewer*) |

### 2a. The H0 findings are a degree artifact

This is the part that changes the write-up. Degree-preserving rewiring reproduces the
day-spanning connected component **exactly** — 23.75 h in all 20 replicates, standard
deviation 0.00 — and reproduces the 32 essential H0 classes (32.25 ± 0.55). Any graph
sequence with Medellín's degree sequence has these properties. That is the expected
result for connectivity, which is largely pinned by the degree sequence, and it means
the current `concluding_remarks.md` finding 1 (persistent connectivity, near-total
`frac_in_largest`) and finding 3 (the threshold-dependent essentiality story) are
**descriptions of the degree distribution**, not discoveries about Medellín's spatial
structure. They should not be presented as findings.

Note the order-shuffle null gives a *lower* longest-H0 (17.71 h) than observed — so the
day-spanning component does require the real temporal ordering. But since rewiring
reproduces it exactly, the ordering is doing the work, not the geography.

### 2b. The H1 lifetime finding is real and it is the strongest result

209 loop classes survive longer than 3 hours. Under the order null the mean is 0.3
(maximum 3 across 20 replicates); under the rewire null it is **exactly zero, every
time**. The observed value is roughly 70× the largest value any null replicate produced.
Both nulls preserve per-window topology or per-zone degree respectively, so this is not
volume, not degree, and not the survey's severe minute-heaping — the order null holds
each snapshot fixed, heaping included.

The direction of the total-count result sharpens it: real mobility has **far fewer**
loops (10,569) than either null (26.6k, 37.8k) but they last **much longer**. Randomised
mobility produces a churn of ephemeral cycles; actual mobility produces a smaller number
of stable ones. That is a defensible, non-obvious structural claim.

### 2c. But the zigzag is redundant for counting

The number of H1 classes alive at each window centre correlates with that window's
static β₁ at **r = 0.978**. Anyone could have computed the per-window loop count with
ordinary snapshot-wise homology and no zigzag at all. The zigzag earns its place *only*
through lifetimes and class identity — tracking that a given loop is the same loop across
windows. So the paper, if written, is about **persistence of cycles**, and every panel
that reports counts rather than lifetimes is doing work the method wasn't needed for.

## 2d. First look at the transit overlay — cycles cluster near Metro, but proximity alone is a weak test

`transit_overlay.py` pulls the operator-published Metro/Metrocable/Tranvía geometry
from Medellín's municipal GIS server (`servidormapas`), keeps the 7 lines open by the
2017 survey year (Metro Líneas A/B; Metrocable H, J, K, L; Tranvía T — excluding
Línea M, Línea P, and "La Aldea", which opened 2019–2021, and the Línea 1/2 BRT
pretroncales, which are buses, not rail), and compares the notebook's own
`indicative_h1_cycles()` output for the top 20 longest H1 bars against transit.

The naive test — does the loop's connecting polygon cross a transit line — is
**useless**: random same-size zone loops cross the transit network 95–99% of the time
by pure geometry, because the lines run the length of the dense valley center where
loop density is also highest. A sharper test uses each cycle's **zone centroids**
(vertices, not the connecting edges) and compares mean distance-to-nearest-line against
a null of random same-size zone tuples drawn from the same "ever active" pool the
cycles come from (`SEED = 20170101`, 3,000 draws per cycle size):

| | observed median | null median (k=4/5/6) | 
|---|---|---|
| mean zone distance to Metro (m) | **1,240** | 2,088 / 2,276 / 2,620 |
| mean zone distance to any 2017 line (m) | **623** | 1,815 / 2,018 / 2,307 |

90% of the 40 examined cycles sit closer to transit than their size-matched null
median; the mean percentile rank is 0.18 against Metro alone and 0.12 against the
full network (0 = closer than every null draw, 1 = farther than every null draw). That
is a real signal, not an artifact of the naive crossing test — but it is roughly a
2× distance reduction, not an order-of-magnitude one, and 40 cycles drawn from 20 bars
is a small, non-independent sample (bars overlap in time and often share zones). This
is a first look, not a validated result.

![Persistent H1 cycles against Metro Líneas A/B]({{artifact:1ed28871-693e-4288-b9b2-b5f9a34ed33a}})
![Persistent H1 cycles against the full 2017 rail network]({{artifact:0240915d-9cfc-4fc3-a2f0-cee989dbb7c7}})

What the two figures show: the five longest-lived H1 bars' indicative carriers cluster
tightly around the Línea A corridor in the dense northern core, with two carriers
extending south along the line's southern extension. Adding Metrocable and Tranvía
(second figure) makes little visual difference — none of the top-5 cycles pass near
those lines, which sit off to the sides of the valley. That asymmetry is itself
informative: if the finding holds up, it says the persistent structure tracks Medellín's
**trunk** commuting corridor specifically, not "rail transit" generically.

**What this first look does not yet establish:**
1. Whether the association is specific to the top 5 (or top 20) bars, or holds across
   the full H1 population — the comparison above uses only the longest-lived cycles.
2. Whether it survives once the same rewire/order nulls from §2 are applied to
   *cycle identity*, not just cycle count — a cycle could sit near Metro simply because
   Metro runs through the highest-degree zones, which is the same degree confound as §2a.
3. Formal significance — the 0.18/0.12 mean percentile ranks are suggestive, not a
   calibrated p-value, and 40 cycles from 20 overlapping bars are not 40 independent
   trials.

## 2e. A parameterized, editable notebook closes gaps 1 and 3 — gap 2 is still open

`metro_overlay.ipynb` (built by `build_metro_overlay.py`) generalizes the §2d
analysis into a standalone notebook that takes any zigzag persistence-diagram CSV
and a population-size knob (`N_EXTRA_SHORT_LIVED`, the number of additional
shorter-lived H1 bars to fold in beyond the `> 4h` core) rather than a hardcoded
top-20. It adds a Monte-Carlo **permutation test** on the population-mean
zone-to-transit distance (`N_PERMUTE` random replicates, default 2,000), which is
what closes gap 3 above.

At the default configuration (core `>4h` bars + 15 extra shorter-lived bars = 63
bars, 121 indicative cycles):

| | observed mean (m) | null mean (m) | fold closer | p-value |
|---|---|---|---|---|
| distance to Metro | 1,530 | 3,390 | 2.22× | **p = 0.0005** |
| distance to any 2017 line | 1,058 | 3,058 | 2.89× | **p = 0.0005** |

Re-running with `N_EXTRA_SHORT_LIVED = 60` (108 bars, 195 cycles — nearly triple the
population) gives 2.42× / 2.96× fold-closer with p = 0.003 on 300 permutation
replicates — the same effect size, holding as the population widens well past the
original top-20/top-5. That answers gap 1: the association is not an artifact of
cherry-picking only the very longest-lived bars.

**Gap 2 remains open.** Both nulls above are spatially random, not degree-preserving
— they draw random same-size zone tuples from the same "ever active" pool, but do
not control for zone degree the way the §2 rewire null does for barcode structure.
A positive result here is consistent with "transit tracks well-connected zones" as
much as with "transit causes/reflects the loop structure specifically." Extending
the degree-preserving rewire null to cycle *identity* (not just count) is the
next step before this is a defensible paper claim, and is not implemented in
`metro_overlay.ipynb`.

## 3. Verdict

**Not a methods paper.** Myers et al. own that ground, including the transportation
application.

**Possibly a short applied paper**, on a much narrower claim than the current write-up
makes: *in a demand-side OD survey, cyclic mobility structure is sparser but far more
temporally persistent than any degree- or volume-matched null, while connectivity
structure is entirely degree-explained.* That is one finding, it is null-robust, and it
carries a genuine methodological caution for the TDA-on-mobility literature — that H0
results in this setting need a degree null before they mean anything. Given that Myers
et al. report H0 findings on GB transportation with no null model, that caution has a
target and is worth stating.

**What it is not yet.** Three gaps stand between this and a submittable paper:

1. **One city, one day, no comparison.** Everything is a cross-section of the 2017
   Medellín survey. There is no second day, no second city, no before/after. Descriptive
   topology of a single OD survey is a workshop note, not a journal article.
2. **The loops are not tied to anything exogenous.** The 209 persistent cycles are the
   result, but nothing yet says *what they are* — whether they trace the Metro/Metrocable
   corridors, ring roads, socioeconomic strata boundaries, or the survey's own zone
   geometry. Until a persistent cycle is matched to a named piece of urban form, the
   finding has no urban-planning content and no reader outside TDA.
3. **The measurement floor is shaky and should be stated.** 74% of departure times heap
   on `:00`/`:30`, against 60-minute windows stepping 30 minutes. The null models show
   the H1 result survives this, but a reviewer will ask, and the `t_d >= lo` boundary
   convention (which moved the bar count from 14,947 to 16,881) has to be in the methods
   section, not a code comment.

**The single highest-value next step** is (2) — overlay the persistent H1 generators on
Medellín's transit network and strata map. It is cheap, the notebook already has
`indicative_h1_cycles()` and a spatial map figure, and it is the difference between "the
barcode has long bars" and "mobility in Medellín sustains a stable cyclic structure along
*these specific corridors*." If those cycles land on the Metro corridors, there is a
paper. If they land on nothing interpretable, the honest conclusion is that this is a
well-built replication with a negative result — worth keeping as a chapter or an
appendix, not worth forcing into a journal.

## Reproducing

```bash
python null_models.py     # ~45 min; writes outputs/null_models.csv, outputs/null_summary.json
```
