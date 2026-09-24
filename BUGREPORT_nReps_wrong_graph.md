# Bug report: `hafez_TI_LINEAR_BRANCH` uses the wrong graph when `nReps > 1`

> **STATUS: fixed (2026-09-23).** Every claim above was reproduced before fixing:
> ElPiGraph returns `nReps + 1` graphs with `ReplicaID: 1 2 3 4 5 0` and
> `ProbPoint: 0.6 x5, 1`, and on the pre-fix build `node.df` matched replicate 1
> exactly (max coordinate difference 0.0000) while differing from the consensus.
>
> Fixed by a shared `epg_consensus()` helper selecting on `ReplicaID == 0` with a
> last-element fallback, applied at all five `TreeEPG[[1]]` sites **and** in
> `hafez_TI_circle`, which previously used `[[length()]]` — so the two functions
> can no longer drift apart. Note three of the five sites are in
> `hafez_lineages_from_root`, a different function, which had to be fixed too or
> the projection would have used different node positions than the path
> structure indexing it.
>
> Verified by `verify_nReps_consensus.R`: 20 checks, all passing; the key ones
> fail against the pre-fix build.
>
> **One caveat on the open empirical question.** The report rightly flags that
> whether the consensus gives *better pseudotime* than a good replicate is
> untested. A crude check across 4 seeds found the consensus node positions
> *more* variable than replicate 1's (0.6064 vs 0.5387) on synthetic data, so do
> not assume the consensus is more stable. The fix is about correctness of
> indexing — returning the graph the API designates — not about that.
>
> The three "Also worth doing" items (export/wrapper, default `nReps`/`ProbPoint`,
> `RETURN_ALL_REPLICATES`) are API changes and were deliberately **not** made.


**File:** `R/traj_inf_fxns.R`
**Function:** `hafez_TI_LINEAR_BRANCH` (internal, `@noRd`, reached via `hafez:::`)
**Affects:** both `BRANCH_TYPE = 'curve'` and `BRANCH_TYPE = 'tree'`
**Fires on the default arguments** (`nReps = 5`, `ProbPoint = 0.6`) — see "Severity".
**Found:** 2026-09-23, during BCP-ALL CyTOF B-cell trajectory work
**Package version:** hafez 0.1.0 · branch `fix/ti-dbpn-mahalanobis-bugs` @ `dbe7af5`
**ElPiGraph.R:** as installed at `/Library/Frameworks/R.framework/Versions/4.6/Resources/library/ElPiGraph.R`

---

## Summary

When `nReps > 1`, ElPiGraph returns a list of `nReps + 1` graphs. The **last** element is
the consensus graph built from all replicates; elements `1..nReps` are the individual
bootstrap replicates, each fitted to a random `ProbPoint` fraction of the data.

`hafez_TI_LINEAR_BRANCH` indexes **`[[1]]`**, so it silently returns *the first bootstrap
replicate* instead of the consensus. With the default `ProbPoint = 0.6` that replicate was
fitted to a random ~60% of the landmarks.

`hafez_TI_circle`, in the same file, indexes `[[length(...)]]` and is **correct**. The two
functions disagree with each other.

---

## Evidence 1 — the ElPiGraph return contract

`computeElasticPrincipalCurve` / `computeElasticPrincipalTree` both delegate to
`ElPiGraph.R:::computeElasticPrincipalGraphWithGrammars`. In that function:

Each replicate is fitted to a random subsample when `ProbPoint < 1` (lines 34-41):

```r
if (ProbPoint < 1 & ProbPoint > 0) {
    SelPoints <- lapply(as.list(1:nReps), function(i) {
      return(runif(nrow(X)) <= ProbPoint)      # each replicate sees ~ProbPoint of X
```

Then, **only when `nReps > 1`**, one extra graph is appended (lines 318-357):

```r
if (nReps > 1) {
    print("Constructing average tree")
    AllPoints <- do.call(rbind, lapply(ReturnList[...], "[[", "NodePositions"))
    ...
    ReturnList[[length(ReturnList) + 1]] <- ElPiGraph.R:::computeElasticPrincipalGraph(
        Data = AllPoints, ...)
    ReturnList[[length(ReturnList)]]$SubSetID   <- j
    ReturnList[[length(ReturnList)]]$ReplicaID  <- 0      # 0 marks the consensus
    ReturnList[[length(ReturnList)]]$ProbPoint  <- 1      # consensus spans all points
}
return(ReturnList)
```

So: `length(ReturnList) == nReps + 1` when `nReps > 1`, and the consensus is **last**,
identifiable by `ReplicaID == 0`. When `nReps == 1` the block is skipped and
`length(ReturnList) == 1`, so `[[1]]` is correct in that case only.

## Evidence 2 — the defect in hafez

`R/traj_inf_fxns.R`, all five sites take `[[1]]`:

| line | code |
|---|---|
| 108 | `node.df = TreeEPG[[1]]$NodePositions %>% ...` |
| 127 | `Tree_Graph <- ElPiGraph.R::ConstructGraph(TreeEPG[[1]])` |
| 129 | `NodeLabs <- 1:nrow(TreeEPG[[1]]$NodePositions)` |
| 163 | `PartStruct <- ElPiGraph.R::PartitionData(X = OOS_DATA_features, NodePositions = TreeEPG[[1]]$NodePositions)` |
| 166-167 | `NodePositions = TreeEPG[[1]]$NodePositions`, `Edges = TreeEPG[[1]]$Edges$Edges` |

## Evidence 3 — the same file already does it correctly

`hafez_TI_circle`, lines 39-46, uses the last element:

```r
PartStruct <- ElPiGraph.R::PartitionData(
    X = OOS_DATA,
    NodePositions = CircleEPG.Boot[[length(CircleEPG.Boot)]]$NodePositions)
```

This is the intended pattern. `hafez_TI_LINEAR_BRANCH` does not follow it.

---

## Severity

The **defaults are `nReps = 5, ProbPoint = 0.6`** (line 73). Any user calling
`hafez_TI_LINEAR_BRANCH` without overriding both gets:

- **Wrong graph.** 6 graphs are fitted; the one returned is replicate #1, fitted to a
  random ~60% subsample. The requested bootstrap consensus is computed, then discarded.
- **Non-reproducible without a seed.** Which 60% lands in replicate #1 is drawn from
  `runif`, so repeated calls give materially different trajectories.
- **Wasted compute.** `nReps + 1` graphs are fitted and 1 is used — `(nReps+1)x` the
  necessary work. In my runs the overhead was roughly 10x at the settings I used.

Two regimes, distinguish them when triaging:

| `ProbPoint` | `nReps` | outcome |
|---|---|---|
| `< 1` | `> 1` | **Incorrect** — returns a subsample fit, not the consensus |
| `= 1` | `> 1` | Statistically OK (all replicates see all data) but `(nReps+1)x` wasted |
| any | `= 1` | Correct — no consensus is appended, `[[1]]` is the only graph |

---

## Reproduction

```r
library(ElPiGraph.R)
set.seed(1)
X  <- cbind(PC1 = sort(rnorm(600)), PC2 = rnorm(600))
LM <- data.frame(X, celltype = rep(letters[1:6], each = 100))

E <- ElPiGraph.R::computeElasticPrincipalCurve(
       X = as.matrix(LM[, c("PC1","PC2")]),
       NumNodes = 20, nReps = 5, ProbPoint = 0.6,
       Lambda = 0.01, Mu = 0.01, Do_PCA = FALSE, verbose = FALSE)

length(E)                                  # 6  == nReps + 1
sapply(E, `[[`, "ReplicaID")               # 1 2 3 4 5 0   <- 0 is the consensus
sapply(E, `[[`, "ProbPoint")               # .6 .6 .6 .6 .6 1

# hafez takes E[[1]] (ReplicaID 1, ProbPoint 0.6); it should take E[[6]] (ReplicaID 0).
identical(E[[1]]$NodePositions, E[[length(E)]]$NodePositions)   # FALSE
```

---

## Proposed fix

Resolve the index once, near the top of the post-fit block in
`hafez_TI_LINEAR_BRANCH` (just before current line 108), and use it at all five sites:

```r
## ElPiGraph appends the consensus graph LAST when nReps > 1 (ReplicaID == 0).
## When nReps == 1 no consensus is appended and length(TreeEPG) == 1.
rid <- vapply(TreeEPG, function(g) if (is.null(g$ReplicaID)) NA_real_
                                   else as.numeric(g$ReplicaID), numeric(1))
BEST <- if (any(rid == 0, na.rm = TRUE)) which(rid == 0)[1] else length(TreeEPG)
EPG  <- TreeEPG[[BEST]]
```

then replace `TreeEPG[[1]]` with `EPG` at lines 108, 127, 129, 163, 166, 167.

Selecting on `ReplicaID == 0` rather than bare `length()` is deliberate: it is explicit
about *which* graph is wanted, and stays correct if the return order ever changes.

### Also worth doing

1. **Export it, or add a wrapper.** Callers currently need `hafez:::` to reach this
   function. Either export it or route users through `hafez_TI`.
2. **Reconsider the defaults.** `ProbPoint = 0.6` with `nReps = 5` costs 6 fits. If the
   consensus is the intended product, that is fine; if not, `nReps = 1, ProbPoint = 1`
   is a cheaper honest default.
3. **Bootstrap averaging is not available through this API.** Even after the fix, the
   returned object is one graph. Callers who want replicate-averaged *pseudotime* (which
   is what actually reduces tied values — see below) must still loop externally with
   `nReps = 1`. A `RETURN_ALL_REPLICATES` argument would make that unnecessary.

### Regression test to add

```r
test_that("TI_LINEAR_BRANCH uses the consensus graph, not replicate 1", {
  set.seed(1)
  LM <- data.frame(PC1 = sort(rnorm(300)), PC2 = rnorm(300),
                   celltype = rep(letters[1:3], each = 100))
  out <- hafez:::hafez_TI_LINEAR_BRANCH(
           LM_DATA = LM, FULL_DATA = LM, FEATURES = c("PC1","PC2"),
           CC_PHASE_COLUMN = "celltype", BRANCH_TYPE = "curve",
           NumNodes = 15, nReps = 5, ProbPoint = 0.6, Do_PCA = FALSE)
  expect_equal(out$TreeEPG[[1]]$ReplicaID, 0)     # consensus, not a subsample replicate
})
```

---

## Scope — what I did and did not verify

**Verified** by reading both sources and the reproduction above: the `[[1]]` vs
`[[length()]]` indexing defect, the ElPiGraph return contract, and the inconsistency with
`hafez_TI_circle`.

**Not a projection bug.** `project_point_onto_graph` / `PartitionData` / `getPseudotime`
are ElPiGraph functions and I have no evidence they misbehave. hafez calls them with the
wrong *graph*; the projection maths itself was not tested and is not implicated.

**Not tested:** `hafez_TI_circle` end to end, `hafez_TI`, the `tree` branch specifically
(the indexing defect applies to it by inspection, but I only exercised `curve`), and
whether the consensus graph is actually better *for pseudotime* than a good replicate —
that is a separate empirical question.

## Workaround currently in use

`~/Downloads/BCP-All_trajectory/optimal_analysis/scripts/08b_mk11_final.R` calls with
`nReps = 1, ProbPoint = 0.9` and loops 16 seeded replicates externally, averaging the
resulting **pseudotime** vectors. `nReps = 1` sidesteps this bug entirely (no consensus is
appended, so `[[1]]` is correct). Averaging pseudotime rather than graphs also cut tied
pseudotime values from 29% to 3.6% on 5.1M cells, which averaging node positions does not
achieve.
