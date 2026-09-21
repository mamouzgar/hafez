# =============================================================================
# Verification for the bugs in HAFEZ_BUG_REPORT.md, plus the three found while
# fixing them. Run:  Rscript --vanilla verify_hafez_bugs.R
#
# --vanilla matters: the renv project library's ElPiGraph.R declares
# `Imports: distutils`, which is not installed there, so loading it under the
# project .Rprofile fails with "there is no package called 'distutils'". That
# is an environment problem, not a hafez bug.
#
# Every check below is INDEPENDENT. An earlier failure can no longer silently
# skip a later one -- the previous version nested the root-endpoint test inside
# `if (!inherits(r1, "try-error"))`, so the highest-priority bug was skipped
# with no output at all.
# =============================================================================

suppressPackageStartupMessages({library(hafez); library(dplyr)})
cat("hafez built:", as.character(packageDescription("hafez")$Built), "\n")
cat("ggplot2 attached?", "package:ggplot2" %in% search(), "\n\n")

PASS <- 0L; FAIL <- 0L
ok <- function(label, cond, detail = "") {
  good <- isTRUE(cond)
  if (good) PASS <<- PASS + 1L else FAIL <<- FAIL + 1L
  cat(sprintf("  [%s] %-52s %s\n", if (good) "PASS" else "FAIL", label, detail))
}
quiet <- function(expr) suppressMessages(suppressWarnings(
  invisible(capture.output(val <- try(expr, silent = TRUE))))) || TRUE
runq <- function(expr) {
  val <- NULL
  suppressMessages(suppressWarnings(invisible(capture.output(
    val <- try(expr, silent = TRUE)))))
  val
}
errmsg <- function(x) if (inherits(x, "try-error"))
  sub("\n.*", "", conditionMessage(attr(x, "condition"))) else NA_character_

set.seed(1)
n  <- 400
df <- data.frame(PC1 = c(rnorm(n/2, -3), rnorm(n/2, 3)),
                 PC2 = rnorm(n), PC3 = rnorm(n))
df$gate <- rep(c("A", "B"), each = n/2)
df$cell.id <- paste0("c", seq_len(n))
F_ <- c("PC1", "PC2", "PC3")

# -----------------------------------------------------------------------------
cat("BUG A (new) -- branch_type default was a length-3 `if` condition\n")
# -----------------------------------------------------------------------------
r <- runq(hafez_TI(FULL_DATA = df, LM_DATA = df[1:100, ], FEATURES = F_, NumNodes = 10,
                   use_start_label = TRUE, start_label_column_category = c("gate", "A")))
ok("default branch_type runs", !inherits(r, "try-error"), errmsg(r))
r <- runq(hafez_TI(FULL_DATA = df, LM_DATA = df[1:100, ], FEATURES = F_, NumNodes = 10,
                   branch_type = "nope", use_start_label = TRUE,
                   start_label_column_category = c("gate", "A")))
ok("invalid branch_type errors clearly", inherits(r, "try-error") &&
     grepl("should be one of", errmsg(r)), errmsg(r))

# -----------------------------------------------------------------------------
cat("\nBUG 1 -- ROOT must be an end2end ENDPOINT, not the nearest node\n")
# -----------------------------------------------------------------------------
EPG <- runq(hafez:::hafez_TI_LINEAR_BRANCH(
  LM_DATA = df, FULL_DATA = df, FEATURES = F_, CC_PHASE_COLUMN = "gate",
  BRANCH_TYPE = "curve", NumNodes = 10, Lambda = .01, Mu = .01, nReps = 3,
  ProbPoint = 1, MaxNumberOfIterations = 10, Do_PCA = FALSE, verbose = FALSE))
ok("graph fits without ggplot2 attached", !inherits(EPG, "try-error"), errmsg(EPG))

if (!inherits(EPG, "try-error")) {
  eps <- unique(unlist(lapply(EPG$Tree_e2e, function(x) as.numeric(names(x))[c(1, length(x))])))
  interior <- setdiff(seq_len(nrow(EPG$TreeEPG[[1]]$NodePositions)), eps)[1]
  cat(sprintf("       endpoints %s ; testing interior node %d\n",
              paste(sort(eps), collapse = ","), interior))

  ri <- runq(hafez:::hafez_lineages_from_root(EPG, OOS_DATA = df, FEATURES = F_, ROOT = interior))
  ok("interior ROOT now errors loudly", inherits(ri, "try-error") &&
       grepl("not an endpoint", errmsg(ri)), errmsg(ri))

  re <- runq(hafez:::hafez_lineages_from_root(EPG, OOS_DATA = df, FEATURES = F_, ROOT = eps[1]))
  ok("endpoint ROOT still returns pseudotime", !inherits(re, "try-error") &&
       nrow(data.frame(re)) == n && ncol(data.frame(re)) >= 1,
     sprintf("%d x %d", nrow(data.frame(re)), ncol(data.frame(re))))
}

# hafez_TI must now CHOOSE an endpoint, not just tolerate one
res <- runq(hafez_TI(FULL_DATA = df, LM_DATA = df[1:100, ], FEATURES = F_, NumNodes = 10,
                     branch_type = "curve", use_start_label = TRUE,
                     start_label_column_category = c("gate", "A")))
ok("hafez_TI picks a valid root end to end", !inherits(res, "try-error") &&
     nrow(data.frame(res)) == n, errmsg(res))

# -----------------------------------------------------------------------------
cat("\nBUG 2 -- documented LM_DATA index-vector form\n")
# -----------------------------------------------------------------------------
rv <- runq(hafez_TI(FULL_DATA = df, LM_DATA = 1:100, FEATURES = F_, NumNodes = 10,
                    branch_type = "curve", use_start_label = TRUE,
                    start_label_column_category = c("gate", "A")))
ok("vector LM_DATA runs", !inherits(rv, "try-error"), errmsg(rv))
rd <- runq(hafez_TI(FULL_DATA = df, LM_DATA = df[1:100, ], FEATURES = F_, NumNodes = 10,
                    branch_type = "curve", use_start_label = TRUE,
                    start_label_column_category = c("gate", "A")))
ok("data.frame LM_DATA runs", !inherits(rd, "try-error"), errmsg(rd))
ok("both forms agree", !inherits(rv, "try-error") && !inherits(rd, "try-error") &&
     isTRUE(all.equal(data.frame(rv), data.frame(rd))))
re0 <- runq(hafez_TI(FULL_DATA = df, LM_DATA = integer(0), FEATURES = F_, NumNodes = 10,
                     branch_type = "curve"))
ok("empty landmark vector returns NA, no error",
   !inherits(re0, "try-error") && length(re0) == 1 && is.na(re0))

# -----------------------------------------------------------------------------
cat("\nBUG 3 -- hafez_DBPN(bandwidth=) is forwarded to density()\n")
# -----------------------------------------------------------------------------
dd <- data.frame(cell.id = paste0("c", 1:2000),
                 pt = c(rnorm(1000, .3, .02), rnorm(1000, .8, .15)))
dd$pt <- (dd$pt - min(dd$pt)) / diff(range(dd$pt))
a  <- runq(hafez_DBPN(dd, column_to_normalize = "pt", bandwidth = "nrd0", new_dbp_name = "P"))
b  <- runq(hafez_DBPN(dd, column_to_normalize = "pt", bandwidth = "SJ",   new_dbp_name = "P"))
c_ <- runq(hafez_DBPN(dd, column_to_normalize = "pt", adjust.value = 0.1, new_dbp_name = "P"))
ok("bw='nrd0' vs bw='SJ' now DIFFER", !isTRUE(all.equal(a$P, b$P)))
ok("adjust still works",              !isTRUE(all.equal(a$P, c_$P)))
ok("numeric bandwidth accepted",
   !inherits(runq(hafez_DBPN(dd, column_to_normalize = "pt", bandwidth = 0.05,
                             new_dbp_name = "P")), "try-error"))

# -----------------------------------------------------------------------------
cat("\nBUG 4 -- features_for_start_cell_id = NULL no longer yields character(0)\n")
# -----------------------------------------------------------------------------
rn <- runq(hafez_TI(FULL_DATA = df, LM_DATA = df[1:100, ], FEATURES = F_, NumNodes = 10,
                    branch_type = "curve"))
ok("both root selectors NULL -> clear error", inherits(rn, "try-error") &&
     grepl("no way to choose a root", errmsg(rn)), errmsg(rn))
rf <- runq(hafez_TI(FULL_DATA = df, LM_DATA = df[1:100, ], FEATURES = F_, NumNodes = 10,
                    branch_type = "curve", features_for_start_cell_id = "PC1"))
ok("features_for_start_cell_id path runs", !inherits(rf, "try-error"), errmsg(rf))
rb <- runq(hafez_TI(FULL_DATA = df, LM_DATA = df[1:100, ], FEATURES = F_, NumNodes = 10,
                    branch_type = "curve", features_for_start_cell_id = "NOT_A_COLUMN"))
ok("unknown start feature errors clearly", inherits(rb, "try-error") &&
     grepl("not found in the landmark data", errmsg(rb)), errmsg(rb))

# -----------------------------------------------------------------------------
cat("\nBUG 5 (new) -- hafez_DBPN assumed pseudotime was on a 0-1 scale\n")
# -----------------------------------------------------------------------------
for (mult in c(1, 10, 100)) {
  d2 <- dd; d2$pt <- d2$pt * mult
  o <- runq(hafez_DBPN(d2, column_to_normalize = "pt", new_dbp_name = "P"))
  ok(sprintf("pseudotime scaled x%-4g", mult),
     !inherits(o, "try-error") && nrow(o) == nrow(d2) && !any(is.na(o$P)),
     if (inherits(o, "try-error")) errmsg(o) else sprintf("rows %d -> %d", nrow(d2), nrow(o)))
}
d3 <- dd; d3$pt <- d3$pt * 4 - 2          # spans negatives too
o3 <- runq(hafez_DBPN(d3, column_to_normalize = "pt", new_dbp_name = "P"))
ok("pseudotime spanning negatives", !inherits(o3, "try-error") &&
     nrow(o3) == nrow(d3) && !any(is.na(o3$P)), errmsg(o3))

# -----------------------------------------------------------------------------
cat("\nBUG 6 (new) -- within-bin interpolation uses the BIN edges\n")
# -----------------------------------------------------------------------------
for (nb in c(1024L, 64L, 16L)) {
  o <- runq(hafez_DBPN(dd, column_to_normalize = "pt", density_bins = nb,
                       new_dbp_name = "P", RETURN_DENSITY_COORDINATES = TRUE))
  if (inherits(o, "try-error")) { ok(sprintf("density_bins=%d", nb), FALSE, errmsg(o)); next }
  x  <- as.data.frame(o$df.dbpn); g <- o$density_df_01
  ref <- stats::approx(x = g$x, y = g$cumulative.sum, xout = x$pt, rule = 2)$y
  ok(sprintf("density_bins=%-5d matches the true CDF", nb), max(abs(x$P - ref)) < 0.005,
     sprintf("max|diff| = %.4f", max(abs(x$P - ref))))
}

# -----------------------------------------------------------------------------
cat("\nRegression -- DBPN still does what it is supposed to do\n")
# -----------------------------------------------------------------------------
tv <- function(v, nbin = 50) {
  h <- hist(v, breaks = seq(0, 1, length.out = nbin + 1), plot = FALSE)$counts
  sum(abs(h/sum(h) - 1/nbin)) / 2
}
o1 <- runq(hafez_DBPN(dd, column_to_normalize = "pt", adjust.value = 0.05, new_dbp_name = "P"))
ok("row count preserved", nrow(o1) == nrow(dd), sprintf("%d -> %d", nrow(dd), nrow(o1)))
ok("output is monotone in the input",
   !is.unsorted(o1$P[order(o1$pt)]))
ok("output is ~uniform at low smoothing", tv(o1$P) < 0.05, sprintf("TV = %.3f", tv(o1$P)))
ok("output lies in [0,1]", min(o1$P) >= 0 && max(o1$P) <= 1,
   sprintf("[%.3f, %.3f]", min(o1$P), max(o1$P)))

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
