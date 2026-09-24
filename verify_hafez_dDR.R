# =============================================================================
# Verification for hafez_dDR(): both back-ends, on synthetic and real data.
#   Rscript --vanilla verify_hafez_dDR.R
# =============================================================================
suppressPackageStartupMessages(library(hafez))
PASS <- 0L; FAIL <- 0L
ok <- function(l, c, d = "") { g <- isTRUE(c); if (g) PASS <<- PASS + 1L else FAIL <<- FAIL + 1L
  cat(sprintf("  [%s] %-50s %s\n", if (g) "PASS" else "FAIL", l, d)) }
rq <- function(e) { v <- NULL
  suppressMessages(suppressWarnings(invisible(capture.output(v <- try(e, silent = TRUE))))); v }
em <- function(x) if (inherits(x, "try-error")) sub("\n.*", "", conditionMessage(attr(x, "condition"))) else ""

set.seed(1)
G <- 40; per <- 400
cells <- do.call(rbind, lapply(seq_len(G), function(i) {
  mu <- 0.2 + 0.6 * (i - 1) / (G - 1); sdv <- 0.05 + 0.10 * ((i %% 5) / 4)
  data.frame(grp = sprintf("g%02d", i), arm = c("ctl", "trt")[1 + (i %% 2)],
             pst = pmin(pmax(rnorm(per, mu, sdv), 0), 1))
}))

cat("1. defaults\n")
r <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst", groups = c("grp", "arm")))
ok("runs with no method given", !inherits(r, "try-error"), em(r))
ok("defaults to w2", identical(r$method, "w2"))
ok("no clustering by default", is.null(r$dDR$dDR_cluster))
ok("returns 2 dims + groups", all(c("dDR1","dDR2","grp","arm","n_cells") %in% names(r$dDR)))
ok("one row per group", nrow(r$dDR) == G, sprintf("%d rows", nrow(r$dDR)))
ok("reports variance explained", !is.null(r$var_explained),
   sprintf("%.1f%% in 2 dims", 100 * sum(r$var_explained)))

cat("\n2. the method switch\n")
rw <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst", groups = c("grp","arm"), method = "w2"))
rd <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst", groups = c("grp","arm"), method = "dtw"))
ok("method = 'w2' runs",  !inherits(rw, "try-error"), em(rw))
ok("method = 'dtw' runs", !inherits(rd, "try-error"), em(rd))
ok("dtw returns a distmat, w2 returns Q",
   !is.null(rd$distmat) && is.null(rd$Q) && !is.null(rw$Q) && is.null(rw$distmat))
ok("dtw reports no variance explained", is.null(rd$var_explained))
ok("bad method errors clearly",
   inherits(bad <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst",
                                groups = "grp", method = "nope")), "try-error") &&
     grepl("should be one of", em(bad)), em(bad))
## These synthetic groups vary in BOTH location and dispersion, with location
## dominant. The back-ends are expected to DISAGREE on which axis leads: W2
## puts location first, DTW is shift-invariant and demotes it. That difference
## is the reason w2 is the default, so assert it rather than assuming agreement.
if (!inherits(rw, "try-error") && !inherits(rd, "try-error")) {
  st <- stats::aggregate(pst ~ grp, cells, function(v) c(m = mean(v), s = stats::sd(v)))
  tru <- data.frame(g = st$grp, mu = st$pst[, "m"], sdv = st$pst[, "s"])
  w <- rw$dDR[match(tru$g, rw$dDR$grp), ]; d <- rd$dDR[match(tru$g, rd$dDR$grp), ]
  ok("w2 dDR1 = location", abs(cor(w$dDR1, tru$mu)) > 0.99,
     sprintf("|r| = %.3f (vs sd %.3f)", abs(cor(w$dDR1, tru$mu)), abs(cor(w$dDR1, tru$sdv))))
  ok("w2 dDR2 = dispersion", abs(cor(w$dDR2, tru$sdv)) > 0.99,
     sprintf("|r| = %.3f", abs(cor(w$dDR2, tru$sdv))))
  ok("dtw demotes location off dDR1", abs(cor(d$dDR1, tru$mu)) < abs(cor(w$dDR1, tru$mu)),
     sprintf("dtw dDR1 vs mean |r| = %.3f, vs sd |r| = %.3f",
             abs(cor(d$dDR1, tru$mu)), abs(cor(d$dDR1, tru$sdv))))
}

cat("\n3. w2 leading axis is cell-cycle position\n")
if (!inherits(rw, "try-error")) {
  mu <- tapply(cells$pst, paste(cells$grp, cells$arm, sep = "\r"), mean)
  d1 <- rw$dDR$dDR1[match(names(mu), rw$dDR$group)]
  ok("dDR1 tracks mean pseudotime", cor(d1, as.numeric(mu)) > 0.99,
     sprintf("r = %.4f", cor(d1, as.numeric(mu))))
}

cat("\n4. opt-in clustering\n")
rk <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst", groups = c("grp","arm"), k = 4))
ok("k = 4 clusters", !inherits(rk, "try-error") && length(unique(rk$dDR$dDR_cluster)) == 4, em(rk))
ok("centroids are returned", !is.null(rk$centroids) && nrow(rk$centroids) == 4)
ok("k >= n_groups errors",
   inherits(bk <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst",
                               groups = "grp", k = 999)), "try-error"), em(bk))

cat("\n5. ts_input path (no raw cells)\n")
tsin <- do.call(rbind, lapply(split(cells, paste(cells$grp, cells$arm)), function(d) {
  dd <- density(d$pst, n = 100, from = 0, to = 1)
  data.frame(x = dd$x, y = dd$y, grp = d$grp[1], arm = d$arm[1])
}))
rt <- rq(hafez_dDR(ts_input = tsin, groups = c("grp","arm"), approach = "density"))
ok("ts_input + w2 runs", !inherits(rt, "try-error"), em(rt))
rt2 <- rq(hafez_dDR(ts_input = tsin, groups = c("grp","arm"), approach = "density", method = "dtw"))
ok("ts_input + dtw runs", !inherits(rt2, "try-error"), em(rt2))
if (!inherits(rt, "try-error") && !inherits(rw, "try-error")) {
  a <- rw$dDR[order(rw$dDR$group), "dDR1"]; b <- rt$dDR[order(rt$dDR$group), "dDR1"]
  ok("cells and ts_input agree", abs(cor(a, b, method = "spearman")) > 0.99,
     sprintf("|rho| = %.4f", abs(cor(a, b, method = "spearman"))))
}
ok("neither input errors", inherits(rq(hafez_dDR()), "try-error"))

cat("\n6. reproducibility and grid size\n")
r1 <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst", groups = "grp"))
r2 <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst", groups = "grp"))
ok("deterministic across runs", isTRUE(all.equal(r1$dDR$dDR1, r2$dDR$dDR1)))
r64  <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst", groups = "grp", m = 64))
r512 <- rq(hafez_dDR(cells = cells, pseudotime_column = "pst", groups = "grp", m = 512))
ok("m = 64 matches m = 512", abs(cor(r64$dDR$dDR1, r512$dDR$dDR1)) > 0.9999,
   sprintf("|r| = %.6f", abs(cor(r64$dDR$dDR1, r512$dDR$dDR1))))

cat("\n7. real data (if present)\n")
TAB <- "/Users/amouzgar/postdoc_projects/Hafez_manuscript/cell_cycle_t_cell_atlas/datasets_investigation/tables"
f <- file.path(TAB, "downsampled_cells_CAR_T_dasatinib.rds")
if (file.exists(f)) {
  d <- readRDS(f); d <- d[!is.na(d$pseudotime_dbpn), ]
  t1 <- system.time(rr <- rq(hafez_dDR(cells = d, pseudotime_column = "pseudotime_dbpn",
                                       groups = c("condition_f","car_type_f","timepoint_f","donor_f"))))
  ok("real data, w2", !inherits(rr, "try-error"),
     sprintf("%d groups, %.1f%% in 2 dims, %.2f s", nrow(rr$dDR),
             100*sum(rr$var_explained), t1[["elapsed"]]))
  t2 <- system.time(rr2 <- rq(hafez_dDR(cells = d, pseudotime_column = "pseudotime_dbpn",
                                        groups = c("condition_f","car_type_f","timepoint_f","donor_f"),
                                        method = "dtw")))
  ok("real data, dtw", !inherits(rr2, "try-error"), sprintf("%.2f s", t2[["elapsed"]]))
  if (!inherits(rr,"try-error") && !inherits(rr2,"try-error"))
    ok("real data: both agree on ordering",
       abs(cor(rr$dDR$dDR1[order(rr$dDR$group)], rr2$dDR$dDR1[order(rr2$dDR$group)],
               method="spearman")) > 0.8,
       sprintf("|rho| = %.3f", abs(cor(rr$dDR$dDR1[order(rr$dDR$group)],
               rr2$dDR$dDR1[order(rr2$dDR$group)], method="spearman"))))
} else cat("  [SKIP] real data not available here\n")

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
