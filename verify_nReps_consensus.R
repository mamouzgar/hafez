# =============================================================================
# Verification for BUGREPORT_nReps_wrong_graph.md
#   Rscript --vanilla verify_nReps_consensus.R
#
# --vanilla: the renv project library's ElPiGraph.R imports `distutils`, which
# is not installed there. Environment issue, not a hafez bug.
# =============================================================================
suppressPackageStartupMessages({library(hafez); library(dplyr); library(ElPiGraph.R)})
cat("hafez built:", as.character(packageDescription("hafez")$Built), "\n\n")

PASS <- 0L; FAIL <- 0L
ok <- function(label, cond, detail = "") {
  good <- isTRUE(cond); if (good) PASS <<- PASS + 1L else FAIL <<- FAIL + 1L
  cat(sprintf("  [%s] %-54s %s\n", if (good) "PASS" else "FAIL", label, detail))
}
runq <- function(expr) { v <- NULL
  suppressMessages(suppressWarnings(invisible(capture.output(v <- try(expr, silent = TRUE))))); v }

mkLM <- function(n = 300) data.frame(PC1 = sort(rnorm(n)), PC2 = rnorm(n),
                                     celltype = rep_len(letters[1:3], n),
                                     cell.id = paste0("c", seq_len(n)))
FIT <- function(LM, ...) hafez:::hafez_TI_LINEAR_BRANCH(
  LM_DATA = LM, FULL_DATA = LM, FEATURES = c("PC1","PC2"),
  CC_PHASE_COLUMN = "celltype", BRANCH_TYPE = "curve", NumNodes = 15,
  Lambda = .01, Mu = .01, MaxNumberOfIterations = 10, Do_PCA = FALSE,
  verbose = FALSE, ...)

# -----------------------------------------------------------------------------
cat("1. epg_consensus() selection rule\n")
# -----------------------------------------------------------------------------
HAS_HELPER <- exists("epg_consensus", envir = asNamespace("hafez"), inherits = FALSE)
if (!HAS_HELPER) {
  cat("  [SKIP] epg_consensus() is not present in this build -- section skipped\n")
  FAIL <- FAIL + 1L
}
fake <- list(list(ReplicaID = 1, NodePositions = 1), list(ReplicaID = 2, NodePositions = 2),
             list(ReplicaID = 0, NodePositions = 99))
consensus_of <- function(x) if (HAS_HELPER) hafez:::epg_consensus(x) else x[[1]]  # pre-fix behaviour
if (HAS_HELPER) {
ok("picks ReplicaID == 0", hafez:::epg_consensus(fake)$NodePositions == 99)
ok("falls back to last when untagged",
   hafez:::epg_consensus(list(list(NodePositions = 1), list(NodePositions = 7)))$NodePositions == 7)
ok("single graph (nReps = 1)", hafez:::epg_consensus(list(list(ReplicaID = 1, NodePositions = 5)))$NodePositions == 5)
ok("consensus not last is still found",
   hafez:::epg_consensus(list(list(ReplicaID = 0, NodePositions = 42),
                              list(ReplicaID = 1, NodePositions = 1)))$NodePositions == 42)
ok("empty input errors", inherits(try(hafez:::epg_consensus(list()), silent = TRUE), "try-error"))
}

# -----------------------------------------------------------------------------
cat("\n2. hafez_TI_LINEAR_BRANCH uses the consensus, not replicate 1\n")
# -----------------------------------------------------------------------------
set.seed(1); LM <- mkLM()
set.seed(1); out <- runq(FIT(LM, nReps = 5, ProbPoint = 0.6))
ok("fit succeeds", !inherits(out, "try-error"))

if (!inherits(out, "try-error")) {
  E <- out$TreeEPG
  ok("ElPiGraph returned nReps + 1 graphs", length(E) == 6, sprintf("length = %d", length(E)))
  ok("last graph is the consensus", isTRUE(E[[length(E)]]$ReplicaID == 0),
     paste("ReplicaID:", paste(sapply(E, function(g) g$ReplicaID), collapse = " ")))

  cons <- (if (HAS_HELPER) hafez:::epg_consensus(E) else E[[length(E)]])$NodePositions
  rep1 <- E[[1]]$NodePositions
  used <- as.matrix(out$node.df[, c("X1","X2")])   # node.df names coords X1/X2
  dimnames(used) <- NULL; dimnames(cons) <- NULL; dimnames(rep1) <- NULL

  ok("node.df matches the CONSENSUS graph", isTRUE(all.equal(used, cons, tolerance = 1e-10)))
  ok("node.df does NOT match replicate 1", !isTRUE(all.equal(used, rep1, tolerance = 1e-10)),
     sprintf("max coord diff vs rep1 = %.4f", max(abs(used - rep1))))
  ok("consensus really differs from replicate 1", !identical(cons, rep1))

  # the graph object handed downstream must match too
  ok("Tree_Graph has the consensus node count",
     igraph::vcount(out$Tree_Graph) == nrow(cons))
  ok("NodeLabs length matches consensus", length(out$NodeLabs) == nrow(cons))
}

# -----------------------------------------------------------------------------
cat("\n3. end-to-end projection uses the same graph\n")
# -----------------------------------------------------------------------------
if (!inherits(out, "try-error")) {
  eps  <- unique(unlist(lapply(out$Tree_e2e, function(p) as.numeric(names(p))[c(1, length(p))])))
  pst  <- runq(hafez:::hafez_lineages_from_root(out, OOS_DATA = LM,
                                                FEATURES = c("PC1","PC2"), ROOT = eps[1]))
  ok("pseudotime computes", !inherits(pst, "try-error") && nrow(data.frame(pst)) == nrow(LM))

  # independent reference: project onto the consensus with ElPiGraph directly
  cg   <- if (HAS_HELPER) hafez:::epg_consensus(out$TreeEPG) else out$TreeEPG[[length(out$TreeEPG)]]
  Xm   <- as.matrix(LM[, c("PC1","PC2")])
  part <- ElPiGraph.R::PartitionData(X = Xm, NodePositions = cg$NodePositions)
  proj <- ElPiGraph.R::project_point_onto_graph(X = Xm, NodePositions = cg$NodePositions,
                                                Edges = cg$Edges$Edges, Partition = part$Partition)
  path <- out$Tree_e2e[[which(sapply(out$Tree_e2e,
            function(x) any(as.numeric(names(x))[c(1, length(x))] == eps[1])))[1]]]
  if (as.numeric(names(path))[1] != eps[1]) path <- rev(path)
  ref <- dynutils::scale_minmax(ElPiGraph.R::getPseudotime(ProjStruct = proj,
                                                           NodeSeq = names(path))$Pt)
  ok("matches an independent projection onto the consensus",
     isTRUE(all.equal(as.numeric(data.frame(pst)[[1]]), as.numeric(ref), tolerance = 1e-8)))

  # ... and NOT a projection onto replicate 1. Note the consensus path may not
  # even exist in replicate 1's topology, in which case getPseudotime() errors
  # -- that is itself proof the two graphs differ, so accept either outcome.
  p1    <- out$TreeEPG[[1]]
  part1 <- ElPiGraph.R::PartitionData(X = Xm, NodePositions = p1$NodePositions)
  proj1 <- ElPiGraph.R::project_point_onto_graph(X = Xm, NodePositions = p1$NodePositions,
                                                 Edges = p1$Edges$Edges, Partition = part1$Partition)
  r1 <- try(dynutils::scale_minmax(ElPiGraph.R::getPseudotime(
              ProjStruct = proj1, NodeSeq = names(path))$Pt), silent = TRUE)
  if (inherits(r1, "try-error")) {
    ok("differs from replicate 1", TRUE, "consensus path absent from replicate 1's graph")
  } else {
    d <- suppressWarnings(max(abs(as.numeric(data.frame(pst)[[1]]) - as.numeric(r1)), na.rm = TRUE))
    ok("differs from replicate 1", is.finite(d) && d > 1e-6, sprintf("max diff = %.4f", d))
  }
}

# -----------------------------------------------------------------------------
cat("\n4. nReps = 1 still correct (no consensus is appended)\n")
# -----------------------------------------------------------------------------
set.seed(2); LM2 <- mkLM()
set.seed(2); o1 <- runq(FIT(LM2, nReps = 1, ProbPoint = 1))
ok("nReps = 1 fit succeeds", !inherits(o1, "try-error"))
if (!inherits(o1, "try-error")) {
  ok("only one graph returned", length(o1$TreeEPG) == 1, sprintf("length = %d", length(o1$TreeEPG)))
  u <- as.matrix(o1$node.df[, c("X1","X2")]); dimnames(u) <- NULL
  g <- o1$TreeEPG[[1]]$NodePositions; dimnames(g) <- NULL
  ok("node.df matches the single graph", isTRUE(all.equal(u, g, tolerance = 1e-10)))
}

# -----------------------------------------------------------------------------
cat("\n5. the consensus is the stabler choice across seeds\n")
# -----------------------------------------------------------------------------
spread <- function(which_graph) {
  cs <- lapply(1:4, function(s) { set.seed(s + 100)
    o <- runq(FIT(LM, nReps = 5, ProbPoint = 0.6))
    if (inherits(o, "try-error")) return(NULL)
    np <- if (which_graph == "consensus") (if (HAS_HELPER) hafez:::epg_consensus(o$TreeEPG) else o$TreeEPG[[length(o$TreeEPG)]])$NodePositions
          else o$TreeEPG[[1]]$NodePositions
    np[order(np[,1]), , drop = FALSE]      # order-invariant comparison
  })
  cs <- Filter(Negate(is.null), cs)
  mean(sapply(2:length(cs), function(i) mean(abs(cs[[i]] - cs[[1]]))))
}
## INFORMATIONAL, not an assertion. The bug is that hafez returned a graph the
## ElPiGraph API does not designate as the result -- an indexing defect. Whether
## the consensus also yields *better pseudotime* than a single good replicate is
## a separate empirical question that the bug report explicitly left open, and
## this crude node-position spread does not settle it (independent fits are not
## node-aligned, so the comparison is only indicative).
sc <- spread("consensus"); sr <- spread("rep1")
cat(sprintf("       node-position spread across 4 seeds: consensus %.4f, replicate1 %.4f\n", sc, sr))
cat(sprintf("       -> %s. Informational only; not a pass/fail criterion.\n",
            if (sc < sr) "consensus more stable here" else
              "consensus NOT more stable on this synthetic data"))

# -----------------------------------------------------------------------------
cat("\n6. hafez_TI end to end still works\n")
# -----------------------------------------------------------------------------
set.seed(3); LM3 <- mkLM(400)
r <- runq(hafez_TI(FULL_DATA = LM3, LM_DATA = LM3, FEATURES = c("PC1","PC2"),
                   NumNodes = 15, branch_type = "curve", use_start_label = TRUE,
                   start_label_column_category = c("celltype", "a")))
ok("hafez_TI runs", !inherits(r, "try-error") && nrow(data.frame(r)) == nrow(LM3),
   if (inherits(r, "try-error")) sub("\n.*", "", conditionMessage(attr(r, "condition"))) else "")

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
