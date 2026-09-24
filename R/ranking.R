#' Rank treatments while retaining uncertainty information
#'
#' @param fit An `nmaflow_fit`.
#' @param small_values Direction setting passed to netmeta.
#' @param probabilities For multinma, return posterior rank probabilities instead of ranks.
#' @param ... Backend arguments.
#' @export
nma_rank <- function(fit, small_values = "undesirable", probabilities = TRUE, ...) {
  if (!inherits(fit, "nmaflow_fit")) stop("Expected `nmaflow_fit`.", call. = FALSE)
  if (fit$engine %in% c("netmeta")) return(netmeta::netrank(fit$fit, small.values = small_values, ...))
  if (fit$engine == "multinma") {
    if (isTRUE(probabilities)) return(multinma::posterior_rank_probs(fit$fit, ...))
    return(multinma::posterior_ranks(fit$fit, ...))
  }
  if (fit$engine == "NMA") return(NMA::nmarank(fit$setup, ...))
  if (fit$engine %in% c("MBNMAdose", "MBNMAtime")) {
    funs <- getNamespaceExports(fit$engine)
    cand <- intersect(c("rank", "rank.mbnma", "rank.mbnma.predict"), funs)
    if (!length(cand)) stop("Ranking entry point is backend/version specific; use the fitted backend object directly.", call. = FALSE)
    return(getExportedValue(fit$engine, cand[[1L]])(fit$fit, ...))
  }
  stop("Ranking not implemented for this backend.", call. = FALSE)
}

#' Summarize uncertainty in treatment ranking
#'
#' @param x Result from [nma_rank()] or [nma_boot_rank()].
#' @return A backend-neutral summary where feasible.
#' @export
nma_rank_uncertainty <- function(x) {
  if (inherits(x, "nmaflow_boot")) {
    reps <- x$replicates[x$success, , drop = FALSE]
    return(data.frame(
      treatment = colnames(reps),
      mean = colMeans(reps, na.rm = TRUE),
      sd = apply(reps, 2L, stats::sd, na.rm = TRUE),
      q025 = apply(reps, 2L, stats::quantile, probs = 0.025, na.rm = TRUE),
      q975 = apply(reps, 2L, stats::quantile, probs = 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }
  if (inherits(x, "nma_rank_probs")) return(as.data.frame(x))
  if (is.matrix(x) || is.data.frame(x)) return(as.data.frame(x))
  if (is.list(x)) return(x)
  stop("Unsupported ranking object.", call. = FALSE)
}

#' Multi-criteria decision utility calculated after NMA
#'
#' The utility layer is deliberately separate from statistical synthesis. It never changes
#' treatment-effect estimates.
#'
#' @param estimates Data frame with one row per treatment.
#' @param treatment Treatment-name column.
#' @param criteria Numeric criterion columns.
#' @param weights Named numeric weights matching `criteria`.
#' @param directions Named vector with `higher` or `lower` for each criterion.
#' @param standardize Standardize each criterion to 0-1 before weighting.
#' @return Data frame sorted by utility.
#' @export
nma_utility <- function(estimates, treatment, criteria, weights,
                        directions = stats::setNames(rep("higher", length(criteria)), criteria),
                        standardize = TRUE) {
  .nma_assert_data(estimates); .nma_assert_cols(estimates, c(treatment, criteria))
  if (is.null(names(weights)) || !all(criteria %in% names(weights))) stop("`weights` must be named for all criteria.", call. = FALSE)
  if (!all(criteria %in% names(directions))) stop("`directions` must be named for all criteria.", call. = FALSE)
  w <- weights[criteria] / sum(abs(weights[criteria]))
  z <- estimates[criteria]
  for (nm in criteria) {
    v <- z[[nm]]
    if (standardize) {
      r <- range(v, na.rm = TRUE)
      v <- if (diff(r) == 0) rep(0.5, length(v)) else (v - r[1L]) / diff(r)
    }
    if (directions[[nm]] == "lower") v <- 1 - v
    z[[nm]] <- v
  }
  utility <- Reduce(`+`, Map(function(v, ww) v * ww, z, w))
  out <- data.frame(treatment = estimates[[treatment]], utility = utility, stringsAsFactors = FALSE)
  out[order(out$utility, decreasing = TRUE), , drop = FALSE]
}

#' Decision-threshold analysis wrapper
#' @param fit Fitted NMA model.
#' @param ... Passed to an available `nmathresh` entry function.
#' @export
nma_threshold <- function(fit, ...) {
  .nma_require("nmathresh", "NMA decision threshold analysis")
  funs <- getNamespaceExports("nmathresh")
  cand <- intersect(c("nma_thresh", "nmathresh", "threshold", "thresholds"), funs)
  if (!length(cand)) {
    stop("The installed nmathresh version does not export a supported threshold entry function. Available exports: ", paste(funs, collapse = ", "), call. = FALSE)
  }
  obj <- if (inherits(fit, "nmaflow_fit")) fit$fit else fit
  args <- list(...)
  if (cand[[1L]] == "nma_thresh") {
    if (!all(c("mean.dk", "lhood", "post") %in% names(args))) {
      stop("`nma_threshold()` with nmathresh requires `mean.dk`, `lhood` and `post` through `...`.", call. = FALSE)
    }
    return(do.call(nmathresh::nma_thresh, args))
  }
  getExportedValue("nmathresh", cand[[1L]])(obj, ...)
}
