.nma_resample_studies <- function(data, study, draw) {
  chunks <- vector("list", length(draw))
  for (i in seq_along(draw)) {
    z <- data[as.character(data[[study]]) == draw[[i]], , drop = FALSE]
    z[[study]] <- paste0(draw[[i]], "__boot", i)
    chunks[[i]] <- z
  }
  do.call(rbind, chunks)
}

.nma_refit_netmeta <- function(fit, data) {
  sp <- fit$specification
  nma_fit_freq(data, sm = fit$measure, common = sp$common, random = sp$random,
               reference = sp$reference, method_tau = sp$method_tau)
}

#' Cluster bootstrap for network meta-analysis
#'
#' Studies, not rows, are sampled with replacement. All contrasts from a sampled multi-arm
#' study therefore remain together. Resampled copies receive unique bootstrap study IDs so
#' that dependence is not collapsed accidentally.
#'
#' @param fit A fitted `nmaflow_fit`, currently with automatic refitting for `netmeta`.
#' @param statistic Function taking a refitted `nmaflow_fit` and returning a named numeric vector.
#' @param B Number of bootstrap replicates selected by the user.
#' @param seed Reproducible random seed.
#' @param study Study column in the contrast data.
#' @param conf Confidence level for percentile intervals.
#' @param keep_fits Keep successful fitted models. Memory intensive.
#' @return An `nmaflow_boot` object.
#' @export
nma_boot <- function(fit, statistic = function(z) stats::coef(z$fit), B = 1999,
                     seed = getOption("nmaFlow.seed", 260819), study = "study",
                     conf = 0.95, keep_fits = FALSE) {
  if (!inherits(fit, "nmaflow_fit")) stop("Expected `nmaflow_fit`.", call. = FALSE)
  if (B < 2L) stop("`B` must be at least 2.", call. = FALSE)
  if (fit$engine != "netmeta") stop("Automatic bootstrap refitting is currently implemented for netmeta-backed fits.", call. = FALSE)
  d <- fit$data; .nma_assert_cols(d, study)
  ids <- unique(as.character(d[[study]]))
  if (length(ids) < 2L) stop("At least two studies are required for study-level bootstrap.", call. = FALSE)
  set.seed(seed)
  theta0 <- statistic(fit)
  theta0 <- as.numeric(theta0); names(theta0) <- names(statistic(fit)) %||% paste0("stat", seq_along(theta0))
  sims <- matrix(NA_real_, nrow = B, ncol = length(theta0), dimnames = list(NULL, names(theta0)))
  failures <- character(B); fits <- if (keep_fits) vector("list", B) else NULL
  for (b in seq_len(B)) {
    draw <- sample(ids, length(ids), replace = TRUE)
    db <- .nma_resample_studies(d, study, draw)
    ans <- tryCatch({
      fb <- .nma_refit_netmeta(fit, db)
      sb <- as.numeric(statistic(fb))
      if (length(sb) != ncol(sims)) stop("Statistic length changed across bootstrap replicates.")
      sims[b, ] <- sb
      if (keep_fits) fits[[b]] <- fb
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(ans)) failures[b] <- ans
  }
  ok <- apply(sims, 1L, function(z) all(is.finite(z)))
  alpha <- (1 - conf) / 2
  ci <- t(apply(sims[ok, , drop = FALSE], 2L, stats::quantile,
                probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE))
  ci <- data.frame(term = rownames(ci), estimate = theta0, lower = ci[, 1L], upper = ci[, 2L],
                   conf = conf, stringsAsFactors = FALSE)
  .nma_new("nmaflow_boot", original = theta0, replicates = sims, success = ok,
           failure_messages = failures[!ok], success_rate = mean(ok), ci = ci,
           B = B, seed = seed, fits = fits, unit = "study")
}

#' @export
print.nmaflow_boot <- function(x, ...) {
  cat("nmaFlow study-level bootstrap\n")
  cat("  B:", x$B, " success rate:", sprintf("%.1f%%", 100 * x$success_rate), "\n")
  print(x$ci, row.names = FALSE)
  if (x$success_rate < 0.95) cat("  Warning: substantial refit failure; inspect nma_boot_diagnose().\n")
  invisible(x)
}

#' Bootstrap coefficient intervals
#' @param fit A netmeta-backed fit.
#' @param B Number of replicates.
#' @param seed Seed.
#' @param ... Passed to [nma_boot()].
#' @export
nma_boot_coef <- function(fit, B = 1999, seed = getOption("nmaFlow.seed", 260819), ...) {
  stat <- function(z) {
    nm <- z$fit
    ref <- nm$reference.group
    if (is.null(ref)) ref <- nm$trts[1]
    te <- nm$TE.random
    trts <- nm$trts[nm$trts != ref]
    vals <- te[ref, trts]
    names(vals) <- paste0(trts, " vs ", ref)
    vals
  }
  nma_boot(fit, statistic = stat, B = B, seed = seed, ...)
}

#' Bootstrap ranking statistics
#' @param fit A netmeta-backed fit.
#' @param B Number of replicates.
#' @param seed Seed.
#' @param small_values Whether lower values indicate better outcomes.
#' @param ... Passed to [nma_boot()].
#' @export
nma_boot_rank <- function(fit, B = 1999, seed = getOption("nmaFlow.seed", 260819),
                          small_values = "undesirable", ...) {
  stat <- function(z) {
    r <- netmeta::netrank(z$fit, small.values = small_values)
    ans <- r$Pscore %||% r$ranking.random %||% r$ranking.common
    if (is.null(ans)) stop("Could not extract a ranking statistic from netmeta::netrank().")
    ans
  }
  nma_boot(fit, statistic = stat, B = B, seed = seed, ...)
}

#' Bootstrap a user-defined prediction statistic
#' @param fit A netmeta-backed fit.
#' @param prediction_function Function returning a named numeric prediction vector.
#' @param B Number of replicates.
#' @param seed Seed.
#' @param ... Passed to [nma_boot()].
#' @export
nma_boot_predict <- function(fit, prediction_function, B = 1999,
                             seed = getOption("nmaFlow.seed", 260819), ...) {
  if (!is.function(prediction_function)) stop("`prediction_function` must be a function.", call. = FALSE)
  nma_boot(fit, statistic = prediction_function, B = B, seed = seed, ...)
}

#' Diagnose bootstrap stability
#' @param x An `nmaflow_boot` object.
#' @return List with success rate, unique failure messages and Monte Carlo summaries.
#' @export
nma_boot_diagnose <- function(x) {
  if (!inherits(x, "nmaflow_boot")) stop("Expected `nmaflow_boot`.", call. = FALSE)
  reps <- x$replicates[x$success, , drop = FALSE]
  list(
    B = x$B,
    successful = sum(x$success),
    failed = sum(!x$success),
    success_rate = x$success_rate,
    failure_messages = sort(table(x$failure_messages), decreasing = TRUE),
    bootstrap_sd = apply(reps, 2L, stats::sd, na.rm = TRUE),
    bootstrap_mcse_mean = apply(reps, 2L, stats::sd, na.rm = TRUE) / sqrt(nrow(reps))
  )
}
