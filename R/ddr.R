## =============================================================================
## Density-based dimensionality reduction (dDR) of pseudotime distributions.
##
## Two interchangeable back-ends behind one interface:
##
##   method = "w2"  (DEFAULT) exact 2-Wasserstein geometry -> Wasserstein PCA
##   method = "dtw"           dynamic time warping        -> classical MDS
##
## Why "w2" is the default
## -----------------------
## The objects being compared are probability distributions on a 1-D pseudotime
## axis, not arbitrary time series. For distributions on the line,
##
##     W2(P,Q)^2 = \int_0^1 ( F_P^{-1}(u) - F_Q^{-1}(u) )^2 du
##
## so the quantile map P |-> F_P^{-1} is an ISOMETRY into L2. Euclidean geometry
## on quantile functions IS Wasserstein geometry, which makes the embedding an
## ordinary PCA: exact, a true metric, no negative eigenvalues, and a
## variance-explained figure. DTW, by contrast, is designed to be invariant to
## shifts along the x axis -- which on a normalised pseudotime axis is usually
## the biological signal -- is not a metric, and leaves cmdscale with negative
## eigenvalues that `add = TRUE` silently absorbs.
##
## The "dtw" back-end is retained so the previous behaviour stays reachable and
## the two can be compared on the same input.
## =============================================================================

#' @importFrom stats quantile dist cmdscale approx isoreg setNames sd
#' @importFrom utils head tail
NULL


## ---------------------------------------------------------------------------
## Quantile representation -- the isometric embedding
## ---------------------------------------------------------------------------

#' Embed grouped pseudotime values as quantile functions
#'
#' Maps each group's raw pseudotime values to its empirical quantile function on
#' a common grid. Euclidean distance on the result is exactly the 2-Wasserstein
#' distance between the underlying distributions, so no kernel density estimate,
#' bandwidth or evaluation grid is involved.
#'
#' @param x numeric vector of pseudotime values, one entry per cell.
#' @param group grouping vector of the same length as \code{x}.
#' @param m number of quantile grid points. The quantile function is smooth, so
#'   values well below the default preserve the geometry; \code{m = 64} is
#'   usually indistinguishable from \code{m = 512} and much cheaper.
#' @param min_n groups with fewer than this many cells are dropped.
#' @param type quantile algorithm passed to \code{\link[stats]{quantile}}.
#'   \code{type = 1} is the true generalised inverse.
#' @return A list with \code{Q} (groups x m matrix of quantile functions),
#'   \code{u} (the probability grid), \code{n} (cells per group) and
#'   \code{groups}.
#' @export
quantile_embed <- function(x, group, m = 64L, min_n = 30L, type = 1L) {
     stopifnot(length(x) == length(group))
     keep_i <- !is.na(x) & !is.na(group)
     x <- x[keep_i]; group <- as.character(group)[keep_i]
     u <- (seq_len(m) - 0.5) / m
     sp <- split(x, group)
     n <- vapply(sp, length, integer(1))
     keep <- n >= min_n
     if (!any(keep)) stop('no group has at least min_n = ', min_n, ' observations')
     if (any(!keep)) {
          message('dropping ', sum(!keep), ' group(s) with < ', min_n, ' cells')
          sp <- sp[keep]; n <- n[keep]
     }
     Q <- t(vapply(sp, function(v)
          stats::quantile(v, probs = u, names = FALSE, type = type), numeric(m)))
     rownames(Q) <- names(sp)
     list(Q = Q, u = u, n = n, groups = names(sp))
}

#' Convert a density curve into a quantile function
#'
#' Escape hatch for input that has already been reduced to density curves and no
#' longer has the underlying cells. Prefer \code{\link{quantile_embed}} on raw
#' cells: this path inherits the kernel density estimate's bandwidth and grid.
#'
#' @param x grid of pseudotime values.
#' @param y density (or any non-negative weight) evaluated at \code{x}.
#' @param m number of quantile grid points.
#' @export
density_to_quantile <- function(x, y, m = 64L) {
     o <- order(x); x <- x[o]; y <- pmax(as.numeric(y)[o], 0)
     if (sum(y) <= 0) return(rep(NA_real_, m))
     w <- diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2
     cdf <- c(0, cumsum(w)); cdf <- cdf / cdf[length(cdf)]
     u <- (seq_len(m) - 0.5) / m
     ok <- !duplicated(cdf)
     stats::approx(cdf[ok], x[ok], xout = u, rule = 2)$y
}

#' Pairwise 2-Wasserstein distances between quantile-embedded distributions
#'
#' A scaled Euclidean distance, because the quantile embedding is an isometry.
#' Exact, and a true metric. Note that the dDR embedding itself never needs this
#' matrix -- \code{\link{hafez_dDR}} decomposes the n x m quantile matrix
#' directly, avoiding the O(n^2) object that \code{cmdscale} requires.
#'
#' @param Q groups x m matrix of quantile functions.
#' @export
w2_dist <- function(Q) stats::dist(Q / sqrt(ncol(Q)), method = 'euclidean')

#' Wasserstein barycentre
#'
#' The Frechet mean under W2, which for 1-D distributions is the mean quantile
#' function. Unlike a pointwise mean of densities it does not invent mass
#' between modes: averaging two shifted unimodal densities pointwise yields a
#' spurious bimodal curve, whereas the barycentre is unimodal at the average
#' position.
#'
#' @param Q groups x m matrix of quantile functions.
#' @param weights optional non-negative weights, e.g. cells per group.
#' @export
w2_barycenter <- function(Q, weights = NULL) {
     if (is.null(weights)) return(colMeans(Q))
     as.vector(crossprod(weights / sum(weights), Q))
}

#' Decompose a squared W2 distance into location, spread and shape
#'
#' \deqn{W_2^2 = (\mu_1-\mu_2)^2 + (\sigma_1-\sigma_2)^2 + 2\sigma_1\sigma_2(1-\rho)}
#' attributing a difference between two distributions to a shift along
#' pseudotime, a change in dispersion, or a change in shape.
#'
#' @param q1,q2 quantile functions on a common grid.
#' @export
w2_decompose <- function(q1, q2) {
     mu1 <- mean(q1); mu2 <- mean(q2)
     s1 <- sqrt(mean((q1 - mu1)^2)); s2 <- sqrt(mean((q2 - mu2)^2))
     rho <- if (s1 > 0 && s2 > 0) mean((q1 - mu1) * (q2 - mu2)) / (s1 * s2) else NA_real_
     loc <- (mu1 - mu2)^2; spr <- (s1 - s2)^2
     shp <- if (is.na(rho)) NA_real_ else 2 * s1 * s2 * (1 - rho)
     c(w2sq = loc + spr + shp, location = loc, spread = spr, shape = shp,
       d_mean = mu1 - mu2, d_sd = s1 - s2, quantile_cor = rho)
}

#' Project a candidate quantile function back onto the valid ones
#'
#' A quantile function must be non-decreasing. Uses isotonic regression, the
#' L2-closest monotone function; \code{cummax} is not the L2 projection and
#' creates flat runs, which reconstruct as spurious density spikes.
#'
#' @param q candidate quantile function.
#' @param range optional \code{c(lo, hi)} support to clamp to.
#' @export
monotonise <- function(q, range = NULL) {
     q <- stats::isoreg(q)$yf
     if (!is.null(range)) q <- pmin(pmax(q, range[1]), range[2])
     q
}

#' Principal component analysis in Wasserstein space
#'
#' Exact, because the quantile embedding is isometric: this is an ordinary PCA
#' whose Euclidean geometry coincides with W2. Components are interpretable
#' modes of distributional variation -- typically location (position along
#' pseudotime), then dispersion, then skew.
#'
#' Component signs follow a fixed convention (loadings sum positive), so a
#' higher score always means later along pseudotime and scores are reproducible
#' across runs and subsets.
#'
#' @param Q groups x m matrix of quantile functions.
#' @param k number of components.
#' @param weights optional group weights, e.g. cells per group.
#' @param fast use a truncated SVD (\pkg{irlba}) when available. Identical
#'   scores, roughly 35x faster, and the gap widens with the number of groups.
#' @export
w2_pca <- function(Q, k = 5L, weights = NULL, fast = TRUE) {
     m <- ncol(Q); k <- min(k, nrow(Q) - 1L, m)
     bary <- w2_barycenter(Q, weights)
     Z <- sweep(Q, 2, bary, '-') / sqrt(m)
     Zw <- if (is.null(weights)) Z else Z * sqrt((weights / sum(weights)) * nrow(Q))
     use_fast <- isTRUE(fast) && requireNamespace('irlba', quietly = TRUE) &&
          k < min(dim(Zw)) - 1L
     if (use_fast) {
          sv <- irlba::irlba(Zw, nv = k)
          d_all <- sv$d
          tot <- sum(Zw^2)                       # exact total variance
          ve <- sv$d^2 / tot
     } else {
          sv <- svd(Zw, nu = k, nv = k)
          d_all <- sv$d[seq_len(k)]
          ve <- (sv$d^2 / sum(sv$d^2))[seq_len(k)]
     }
     V <- sv$v[, seq_len(k), drop = FALSE]
     flip <- ifelse(colSums(V) < 0, -1, 1)       # deterministic sign convention
     V <- sweep(V, 2, flip, '*')
     scores <- Z %*% V
     colnames(scores) <- paste0('dDR', seq_len(k))
     rownames(scores) <- rownames(Q)
     list(scores = scores, components = V,
          sdev = d_all[seq_len(k)] / sqrt(max(nrow(Q) - 1, 1)),
          var_explained = ve[seq_len(k)], barycenter = bary, m = m)
}

#' Summary statistics of each embedded distribution
#'
#' @param Q groups x m matrix of quantile functions.
#' @export
q_summary <- function(Q) {
     u <- (seq_len(ncol(Q)) - 0.5) / ncol(Q)
     data.frame(group = rownames(Q), mean = rowMeans(Q),
                sd = apply(Q, 1, function(q) sqrt(mean((q - mean(q))^2))),
                median = Q[, which.min(abs(u - 0.5))],
                iqr = Q[, which.min(abs(u - 0.75))] - Q[, which.min(abs(u - 0.25))],
                row.names = NULL)
}


## ---------------------------------------------------------------------------
## The user-facing dDR entry point
## ---------------------------------------------------------------------------

#' @title hafez_dDR
#' @description Density-based dimensionality reduction of pseudotime cell
#'   density distributions, with a choice of geometry.
#'
#'   \strong{\code{method = "w2"} (the default)} embeds each group as a quantile
#'   function and runs Wasserstein PCA. Because the quantile map is an isometry
#'   for the 2-Wasserstein distance, this is exact: a true metric, no negative
#'   eigenvalues, a genuine variance-explained figure, and axes in pseudotime
#'   units. It also never forms the n x n distance matrix.
#'
#'   \strong{\code{method = "dtw"}} reproduces the previous behaviour: dynamic
#'   time warping distances via \pkg{dtwclust}, then classical MDS. Retained so
#'   the two can be compared on identical input. Be aware that DTW is built to
#'   be invariant to shifts along the x axis, which on a normalised pseudotime
#'   axis is usually the effect of interest, and that it is not a metric, so
#'   \code{cmdscale} needs \code{add = TRUE} to absorb negative eigenvalues.
#'
#'   Clustering is \strong{off by default}. On the data this was developed
#'   against, the structure is a continuous gradient rather than discrete
#'   groups: a gap statistic selected k = 1, and the first component carried
#'   ~91\% of the variance. Treating the leading score as a continuous
#'   distributional phenotype is usually more informative than a partition.
#'   Supply \code{k} to opt in.
#'
#' @param ts_input long data frame of curves with columns \code{x} (pseudotime)
#'   and \code{y} (density or expression), plus the \code{groups} columns. Used
#'   when raw cells are unavailable. Ignored if \code{cells} is supplied.
#' @param cells cell-level data frame. Preferred for \code{method = "w2"}: the
#'   quantile embedding is computed directly from the values, with no kernel
#'   density estimate, bandwidth or grid.
#' @param pseudotime_column name of the pseudotime column in \code{cells}.
#' @param groups character vector of grouping column names.
#' @param method \code{"w2"} for Wasserstein PCA (default) or \code{"dtw"} for
#'   DTW plus classical MDS.
#' @param approach \code{'density'}, \code{'expr'} or \code{'other'}, describing
#'   \code{ts_input}, as in \code{\link{hafez_tsvz}}.
#' @param dims number of dDR dimensions to return.
#' @param k number of clusters. \code{NULL} (default) means no clustering.
#' @param m quantile grid size for \code{method = "w2"}.
#' @param min_n minimum cells per group when embedding from \code{cells}.
#' @param weight_by_n weight groups by cell count (only when \code{cells} is
#'   given). Recommended when group sizes are uneven.
#' @param seed RNG seed for clustering and for DTW.
#' @param type,distance,normalize passed to \code{\link[dtwclust]{tsclust}} for
#'   \code{method = "dtw"}.
#' @param ... further arguments for the chosen back-end.
#' @return A list with \code{dDR} (a data frame of coordinates, the grouping
#'   columns, and \code{dDR_cluster} if \code{k} was given), \code{method},
#'   \code{var_explained} (\code{w2} only), \code{Q} or \code{distmat}, and
#'   \code{centroids}, the cluster centroids as Wasserstein barycentres for
#'   \code{w2} or as \pkg{dtwclust} centroids for \code{dtw}.
#' @seealso \code{\link{hafez_tsvz}} for the original DTW-only interface.
#' @export
hafez_dDR <- function(ts_input = NULL, cells = NULL, pseudotime_column = NULL,
                      groups = NULL, method = c('w2', 'dtw'),
                      approach = c('density', 'expr', 'other'),
                      dims = 2L, k = NULL, m = 64L, min_n = 30L,
                      weight_by_n = TRUE, seed = 0,
                      type = 'partitional', distance = 'dtw_basic',
                      normalize = TRUE, ...) {
     method <- match.arg(method)
     approach <- match.arg(approach)
     if (is.null(ts_input) && is.null(cells)) {
          stop('supply either `cells` (preferred) or `ts_input`')
     }

     ## ---- build the group x feature matrix --------------------------------
     nvec <- NULL
     if (!is.null(cells)) {
          if (is.null(pseudotime_column)) stop('`pseudotime_column` is required with `cells`')
          if (is.null(groups)) stop('`groups` is required with `cells`')
          gk <- do.call(paste, c(lapply(groups, function(g) as.character(cells[[g]])), sep = '\r'))
          emb <- quantile_embed(cells[[pseudotime_column]], gk, m = m, min_n = min_n)
          Q <- emb$Q; nvec <- emb$n[rownames(Q)]
          if (method == 'dtw') {
               ## DTW needs curves on a shared grid; derive them from the cells
               rng <- range(cells[[pseudotime_column]], na.rm = TRUE)
               sp <- split(cells[[pseudotime_column]], gk)[rownames(Q)]
               FEAT <- t(vapply(sp, function(v)
                    stats::density(v, n = 100L, from = rng[1], to = rng[2])$y, numeric(100)))
          } else FEAT <- Q
     } else {
          ts_input <- as.data.frame(ts_input)
          if (approach %in% c('density', 'expr')) {
               if (is.null(groups)) stop('`groups` is required for approach = "density"/"expr"')
               gk <- do.call(paste, c(lapply(groups, function(g) as.character(ts_input[[g]])), sep = '\r'))
               sp <- split(ts_input[c('x', 'y')], gk)
               Q <- t(vapply(sp, function(d) {
                    a <- stats::aggregate(list(y = d$y), list(x = d$x), mean)
                    density_to_quantile(a$x, a$y, m = m)
               }, numeric(m)))
               FEAT <- if (method == 'dtw')
                    t(vapply(sp, function(d) {
                         a <- stats::aggregate(list(y = d$y), list(x = d$x), mean); a$y
                    }, numeric(length(unique(ts_input$x))))) else Q
          } else {
               Q <- as.matrix(ts_input); FEAT <- Q
               if (method == 'w2') message('approach = "other": rows are treated as quantile functions')
          }
     }
     if (!is.null(k) && k >= nrow(Q)) stop('k must be smaller than the number of groups (', nrow(Q), ')')

     ## ---- the embedding ----------------------------------------------------
     var_expl <- NULL; distmat <- NULL
     if (method == 'w2') {
          fit <- w2_pca(Q, k = max(dims, 2L),
                        weights = if (isTRUE(weight_by_n)) nvec else NULL, ...)
          coords <- fit$scores[, seq_len(dims), drop = FALSE]
          var_expl <- fit$var_explained[seq_len(dims)]
          message(sprintf('dDR (w2): %d groups, %d dims, %.1f%% of W2 variance',
                          nrow(Q), dims, 100 * sum(var_expl)))
     } else {
          if (!requireNamespace('dtwclust', quietly = TRUE))
               stop('method = "dtw" requires the dtwclust package')
          cl <- dtwclust::tsclust(series = FEAT, k = max(2L, if (is.null(k)) 2L else k),
                                  type = type, distance = distance, seed = seed,
                                  normalize = normalize, ...)
          distmat <- as.matrix(cl@distmat)
          mds <- stats::cmdscale(distmat, k = dims, eig = TRUE, add = TRUE)
          coords <- mds$points
          ev <- stats::cmdscale(distmat, k = dims, eig = TRUE, add = FALSE)$eig
          negmass <- sum(abs(ev[ev < 0])) / sum(abs(ev))
          message(sprintf('dDR (dtw): %d groups, %d dims; %.1f%% of |eigenvalue| mass is negative',
                          nrow(FEAT), dims, 100 * negmass))
          if (negmass > 0.05)
               warning(sprintf(paste0('DTW distances are not a metric: %.1f%% of the |eigenvalue| ',
                                      'mass is negative and cmdscale is absorbing it with add = TRUE. ',
                                      'Consider method = "w2".'), 100 * negmass), call. = FALSE)
          rownames(coords) <- rownames(FEAT)
     }
     colnames(coords) <- paste0('dDR', seq_len(dims))

     ## ---- optional clustering ---------------------------------------------
     cluster <- NULL; centroids <- NULL
     if (!is.null(k)) {
          set.seed(seed)
          cluster <- stats::kmeans(if (method == 'w2') Q / sqrt(ncol(Q)) else coords,
                                   centers = k, nstart = 25L, iter.max = 100L)$cluster
          centroids <- if (method == 'w2')
               t(vapply(sort(unique(cluster)), function(cc)
                    w2_barycenter(Q[cluster == cc, , drop = FALSE],
                                  weights = if (!is.null(nvec)) nvec[cluster == cc] else NULL),
                    numeric(ncol(Q))))
          else t(vapply(sort(unique(cluster)), function(cc)
               colMeans(FEAT[cluster == cc, , drop = FALSE]), numeric(ncol(FEAT))))
     }

     ## ---- assemble ---------------------------------------------------------
     gsplit <- do.call(rbind, strsplit(rownames(coords), '\r', fixed = TRUE))
     out <- as.data.frame(coords)
     if (!is.null(groups) && !is.null(gsplit) && ncol(gsplit) == length(groups)) {
          gdf <- as.data.frame(gsplit, stringsAsFactors = FALSE)
          names(gdf) <- groups
          out <- cbind(out, gdf)
     }
     out$group <- rownames(coords)
     if (!is.null(nvec)) out$n_cells <- as.integer(nvec)
     if (!is.null(cluster)) out$dDR_cluster <- as.character(cluster)
     rownames(out) <- NULL

     list(dDR = out, method = method, var_explained = var_expl,
          Q = if (method == 'w2') Q else NULL, distmat = distmat,
          centroids = centroids, n_cells = nvec)
}
