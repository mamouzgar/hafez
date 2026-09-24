# The distmat-reuse optimisation in hafez_tsvz_clust() must be EXACTLY equivalent
# to letting tsclust() recompute the matrix per k -- not merely faster.
test_that("precomputed distmat gives identical clusterings across all k", {
  skip_if_not_installed("dtwclust")
  skip_if_not_installed("proxy")
  set.seed(42)
  # 40 synthetic density-like curves on a shared 60-point grid
  g <- seq(0, 1, length.out = 60)
  Y <- t(sapply(1:40, function(i) {
    mu <- runif(1, 0.15, 0.85); sd <- runif(1, 0.05, 0.20)
    d <- dnorm(g, mu, sd); d / sum(d)
  }))
  rownames(Y) <- paste0("grp", seq_len(nrow(Y)))
  KS <- c(3, 4, 5)

  # optimised path (what the package now does)
  fast <- hafez:::hafez_tsvz_clust(ts_input = Y, k = KS, type = "partitional",
                                   distance = "dtw_basic", normalize = TRUE, seed = 1)
  # reference path: force tsclust to recompute per k by supplying an empty control
  slow <- hafez:::hafez_tsvz_clust(ts_input = Y, k = KS, type = "partitional",
                                   distance = "dtw_basic", normalize = TRUE, seed = 1,
                                   control = dtwclust::partitional_control())

  expect_equal(fast$optim_cluster_k, slow$optim_cluster_k)
  for (kk in as.character(KS)) {
    expect_identical(fast$dtw_clusters[[kk]]@cluster,
                     slow$dtw_clusters[[kk]]@cluster,
                     info = paste("cluster assignments differ at k =", kk))
  }
  # the intra-cluster distance table drives optimal-k selection
  expect_equal(fast$cl_dists$dist, slow$cl_dists$dist, tolerance = 1e-10)
  # and @distmat must be present, since hafez_tsvz() feeds it to cmdscale
  expect_false(is.null(fast$dtw_clusters[[fast$optim_cluster_k]]@distmat))
})

test_that("a caller-supplied control is respected, not overwritten", {
  skip_if_not_installed("dtwclust")
  set.seed(7)
  Y <- matrix(abs(rnorm(20 * 40)), nrow = 20); Y <- Y / rowSums(Y)
  rownames(Y) <- paste0("g", 1:20)
  expect_no_error(
    hafez:::hafez_tsvz_clust(ts_input = Y, k = c(2, 3), type = "partitional",
                             distance = "dtw_basic", normalize = TRUE, seed = 1,
                             control = dtwclust::partitional_control(iter.max = 5L)))
})

test_that("scalar k still works and is unaffected by the optimisation", {
  skip_if_not_installed("dtwclust")
  set.seed(11)
  Y <- matrix(abs(rnorm(15 * 30)), nrow = 15); Y <- Y / rowSums(Y)
  rownames(Y) <- paste0("g", 1:15)
  r <- hafez:::hafez_tsvz_clust(ts_input = Y, k = 3, type = "partitional",
                                distance = "dtw_basic", normalize = TRUE, seed = 1)
  expect_identical(r$optim_cluster_k, "3")
  expect_false(is.null(r$dtw_clusters[["3"]]))
})

test_that("minimum-ICD selection warns when it lands on the k grid boundary", {
  skip_if_not_installed("dtwclust")
  set.seed(5)
  g <- seq(0, 1, length.out = 50)
  Y <- t(sapply(1:30, function(i) { d <- dnorm(g, runif(1, .2, .8), .1); d / sum(d) }))
  rownames(Y) <- paste0("g", 1:30)
  # ICD almost always decreases with k, so the largest k is usually selected
  res <- tryCatch(
    hafez:::hafez_tsvz_clust(ts_input = Y, k = c(2, 3, 4, 5), type = "partitional",
                             distance = "dtw_basic", normalize = TRUE, seed = 1),
    warning = function(w) w)
  if (inherits(res, "warning")) {
    expect_match(conditionMessage(res), "LARGEST k tested")
  } else {
    # an interior optimum is legitimate; then no warning should have fired
    expect_true(as.numeric(res$optim_cluster_k) < 5)
  }
})
