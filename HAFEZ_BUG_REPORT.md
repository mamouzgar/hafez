# HAFEZ bug report

> **STATUS: all fixed (2026-09-21).** Verified by `verify_hafez_bugs.R` — 27 checks,
> which fail 15/27 against the pre-fix build and pass 27/27 after. Three corrections
> to the analysis below, found while fixing:
>
> - **Bug 4 is not low severity.** With `use_start_label = NULL` the `character(0)`
>   root does not give "a silent empty root" — it hard-errors in
>   `dplyr::filter(cell.id == CLOSEST_CELL_ID)` (`..1 must be of size 100 or 1, not
>   size 0`). `hafez_TI()` could not run without `use_start_label` **or**
>   `features_for_start_cell_id`.
> - **Bug 2's scope is wider than stated.** "The dataframe form works" held only for
>   the reporter's configuration. A fifth bug, not in this report, broke both forms:
>   `branch_type = c('curve','tree','circle')` reaching `if (branch_type == 'circle')`
>   with no `match.arg()` — a length-3 condition, fatal in R >= 4.2. Every call using
>   the default `branch_type` failed.
> - **Two further `hafez_DBPN` bugs** were found while checking it against its concept
>   (it is a probability-integral transform, and that part is correct): the bin edges
>   hardcoded `-0.01`/`1.01`, silently assuming pseudotime in [0,1] and dropping cells
>   on any other scale; and within-bin interpolation used `scales::rescale`'s default
>   `from = range(x)` (the observed cells) instead of the bin edges — immaterial at
>   `density_bins = 1024`, growing to 0.05 at 16.
>
> The "Not a bug" ggplot2 call at the end was correct. At default settings the DBPN
> fixes shift existing `pseudotime_dbpn` values by <= 0.0015 (rms 0.0002, r = 1.000000).


Four confirmed bugs in `hafez` 0.1.0, found while running landmark TI + DBPN on a
188k-cell B-ALL scRNA dataset. Each one below was **reproduced against the currently
installed package**, not inferred from reading source. A fifth suspected bug was
investigated and ruled out — see "Not a bug" at the end, and please don't re-chase it.

- Environment: R 4.6.0, hafez 0.1.0 (Built 2026-09-21 04:49:55 UTC), ElPiGraph.R 1.0.0,
  dplyr 1.2.1, macOS (arm64, Accelerate BLAS).
- Source of truth for line numbers: `~/phd/hafez/R/`.
- Runnable reproduction for bugs 1–4: **`verify_hafez_bugs.R`** in this directory
  (`Rscript verify_hafez_bugs.R`). It is self-contained and takes a few seconds.

Suggested priority: **1 > 2 > 4 > 3**.

---

## Bug 1 — `hafez_lineages_from_root()` silently returns nothing unless ROOT is a path endpoint

**Severity: high (silent wrong/empty result).** This is the one that actually cost time.

**Where:** `R/traj_inf_fxns.R:150`, with the trigger at `:320` and the call at `:325`.

**What happens:** `hafez_lineages_from_root()` keeps only paths whose *first or last*
element equals `ROOT`:

```r
SelPaths <- Tree_e2e[sapply(Tree_e2e, function(x){any(x[c(1, length(x))] == ROOT)})]
```

`Tree_e2e` comes from `GetSubGraph(Structure = 'end2end')` (`:128`), so its paths run
between **leaf** nodes. If `ROOT` is an interior node, `SelPaths` is empty, `AllPt` is an
empty list, and `bind_cols()` on it (`:170`) yields a **0-column, 0-row data frame** — with
no error and no warning.

Meanwhile `hafez_TI()` chooses `ROOT` as the node nearest the start-label centroid
(`:320`), over **all** nodes. There is nothing constraining that node to be an endpoint,
and for a real trajectory it frequently isn't (a centroid sits inside its cloud, not at
the tip of the curve).

**Observed:**

```
endpoints: 1,2 of 10 nodes
ROOT = interior node 3 -> ncol = 0 , nrow = 0   <-- silently empty
ROOT = endpoint 1      -> ncol = 1 , nrow = 400
```

**Downstream symptom** is misleading. With `return_pseudotime_only = TRUE` the caller just
gets a 0×0 frame; otherwise you hit the generic `'mismatched data and pseudotime
dataframes'` message and a `NULL` return. Neither mentions the root. My first real run
failed as `Error: nrow(res) == nrow(d) is not TRUE`, which points nowhere near the cause.

**Suggested fix** — snap the start label to the nearest *endpoint* rather than the nearest
node, and fail loudly otherwise:

```r
# in hafez_TI, replace the START_NODE_ID computation (~:320)
endpoints <- unique(unlist(lapply(ELPIGRAPH_RES$Tree_e2e,
                                  function(p) as.numeric(names(p))[c(1, length(p))])))
d2 <- as.matrix(dist(ELPIGRAPH_RES$node.df %>% dplyr::select(any_of(FEATURES))))
START_NODE_ID <- endpoints[which.min(d2[endpoints, nrow(d2)])]

# and guard in hafez_lineages_from_root (~:150)
if (length(SelPaths) == 0) {
  stop("ROOT (", ROOT, ") is not an endpoint of any end2end path. Endpoints are: ",
       paste(endpoints, collapse = ", "))
}
```

**Workaround currently in use:** call the internals directly and root at the end-to-end
endpoint nearest the HSC centroid.

---

## Bug 2 — documented `LM_DATA` vector form always errors

**Severity: high (documented API is unusable).**

**Where:** `R/traj_inf_fxns.R`, the `nrow(LM_DATA) == 0` guard just after the
`is.vector(LM_DATA)` branch (~`:215–222`).

The code explicitly supports index vectors — `## LM data can be either a vector of
indexes or a dataframe` — and branches on `is.vector(LM_DATA)`. But the very next guard
calls `nrow(LM_DATA)`, which is `NULL` for a vector, so `if (NULL == 0)` evaluates
`logical(0)`:

```r
hafez_TI(FULL_DATA = df, LM_DATA = 1:100, FEATURES = c("PC1","PC2","PC3"), ...)
#> Error: argument is of length zero
```

The dataframe form works, so the bug is confined to the documented vector path.

**Fix:** compute the training set first, then check it.

```r
if (is.vector(LM_DATA)) FULL_DATA_TRAIN <- FULL_DATA[LM_DATA, , drop = FALSE]
else                    FULL_DATA_TRAIN <- LM_DATA %>% ungroup()
if (nrow(FULL_DATA_TRAIN) == 0) { message('no landmarks found...'); return(NA) }
```

---

## Bug 3 — `hafez_DBPN(bandwidth=)` is accepted and silently ignored

**Severity: medium (silent no-op on a documented knob).**

**Where:** `R/analysis_fxn.R:318` (signature), with the three `density()` calls at
`:360`, `:414`, `:424`.

`bandwidth = 'nrd0'` is in the signature, but no `density()` call ever passes `bw=` —
only `adjust`. Any value the user supplies does nothing:

```
bw='nrd0' vs bw='SJ' identical?  TRUE   <-- ignored
adjust=1 vs adjust=0.1 identical? FALSE  (adjust does work)
```

This matters in practice: with a sharply peaked terminal population the default bandwidth
under-equalizes badly (total-variation distance from uniform 0.319 at `adjust=1` vs 0.116
at `adjust=0.05`), and the natural thing to reach for — `bandwidth` — appears to work but
changes nothing.

**Fix:** either forward it, `density(..., bw = bandwidth, adjust = adjust.value)`, at all
three sites, or remove the argument so it can't mislead.

---

## Bug 4 — `features_for_start_cell_id = NULL` produces `character(0)`, not an error

**Severity: low (latent; masked whenever `use_start_label` is set).**

**Where:** `R/traj_inf_fxns.R:305–306`.

```r
CLOSEST_CELL_IDX = FULL_DATA_TRAIN %>% ungroup() %>%
  dplyr::select(any_of(features_for_start_cell_id)) %>% apply(., 1, mean) %>% which.min(.)
CLOSEST_CELL_ID = FULL_DATA_TRAIN$cell.id[CLOSEST_CELL_IDX]
```

With `features_for_start_cell_id = NULL` (the default), `any_of(NULL)` selects zero
columns, `apply(..., 1, mean)` gives all-`NaN`, `which.min()` returns `integer(0)`, and
`CLOSEST_CELL_ID` becomes `character(0)`. No error is raised. It happens to be harmless
today only because these two lines run *unconditionally* before the `use_start_label`
branch at `:312` overwrites the choice — so anyone using `use_start_label` never notices,
and anyone not using it gets a silent empty root.

The roxygen docs say NULL "will return object and plot for manual node initialization",
which is not what the code does.

**Fix:** move lines 305–306 inside the `else` branch of the `use_start_label` check, and
`stop()` if `features_for_start_cell_id` is NULL there.

---

## Not a bug — ggplot2 imports (investigated, ruled out)

I originally hit `Error in theme_minimal() : could not find function "theme_minimal"`
from inside `hafez_TI_LINEAR_BRANCH` and assumed the NAMESPACE was missing ggplot2
imports. **That diagnosis was wrong.** The NAMESPACE does import them
(`importFrom(ggplot2, theme_minimal)` and nine siblings), and re-tested against the
current install with ggplot2 deliberately *not* attached, it runs clean:

```
ggplot2 attached? FALSE
--- BUG 1: ggplot2 symbol lookup ---
  OK - no error without ggplot2 attached
```

The installed build carries a timestamp partway through my session, so the failure was
almost certainly a **stale install** predating the ggplot2 imports, fixed by a rebuild.
No source change needed. Please don't spend time here.

---

## Minor note, not filed as a bug

`hafez_DBPN`'s SQL range join uses `BETWEEN` on adjacent, edge-sharing intervals
(`R/analysis_fxn.R:454`), which can emit duplicate rows for values landing exactly on a
boundary. It's currently undone by `distinct(cell.id, .keep_all = TRUE)` at `:463`, `:467`
and `:471`. Row counts were correct in my run (141,453 in → 141,453 out, asserted), so
this is only a robustness note: the dedup silently relies on `cell.id` being unique, and
would quietly drop rows if a genuine duplicate `cell.id` ever arrived. Half-open
intervals (`>= min AND < max`) would remove the need for the dedup.
