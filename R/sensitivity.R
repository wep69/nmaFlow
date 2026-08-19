#' Leave-one-study-out influence analysis
#'
#' @param fit A netmeta-backed `nmaflow_fit`.
#' @param statistic Function extracting a named numeric statistic.
#' @param study Study column name.
#' @return An `nmaflow_influence` object.
#' @export
nma_influence <- function(fit, statistic = function(z) stats::coef(z$fit), study = "study") {
  if (!inherits(fit, "nmaflow_fit") || fit$engine != "netmeta") stop("Currently implemented for netmeta-backed fits.", call. = FALSE)
  d <- fit$data; .nma_assert_cols(d, study)
  ids <- unique(as.character(d[[study]])); base <- as.numeric(statistic(fit)); nm <- names(statistic(fit)) %||% paste0("stat", seq_along(base))
  out <- matrix(NA_real_, nrow = length(ids), ncol = length(base), dimnames = list(ids, nm))
  for (i in seq_along(ids)) {
    z <- d[as.character(d[[study]]) != ids[i], , drop = FALSE]
    ff <- tryCatch(.nma_refit_netmeta(fit, z), error = function(e) NULL)
    if (!is.null(ff)) out[i, ] <- as.numeric(statistic(ff))
  }
  delta <- sweep(out, 2L, base, "-")
  .nma_new("nmaflow_influence", baseline = stats::setNames(base, nm), leave_one_out = out, delta = delta,
           max_abs_change = apply(abs(delta), 1L, max, na.rm = TRUE))
}

#' Outlier-detection orchestration
#' @param fit A fitted model.
#' @param engine `NMAoutlier` or `influence`.
#' @param ... Backend arguments.
#' @export
nma_outliers <- function(fit, engine = c("influence", "NMAoutlier"), ...) {
  engine <- match.arg(engine)
  if (engine == "influence") return(nma_influence(fit, ...))
  .nma_require("NMAoutlier", "network outlier detection")
  funs <- getNamespaceExports("NMAoutlier")
  candidate <- intersect(c("NMAoutlier", "forward.search", "nmaoutlier"), funs)
  if (!length(candidate)) stop("Could not identify the public NMAoutlier entry function in the installed version.", call. = FALSE)
  getExportedValue("NMAoutlier", candidate[[1L]])(fit$fit, ...)
}

#' Run a structured sensitivity set
#'
#' @param fits Named list of defensible `nmaflow_fit` models.
#' @param statistic Function extracting the quantity to compare.
#' @return Long data frame of model-by-statistic values.
#' @export
nma_sensitivity <- function(fits, statistic = function(z) stats::coef(z$fit)) {
  if (!is.list(fits) || !length(fits)) stop("Provide a named list of fitted models.", call. = FALSE)
  if (is.null(names(fits))) names(fits) <- paste0("model", seq_along(fits))
  out <- lapply(names(fits), function(nm) {
    v <- statistic(fits[[nm]]); data.frame(model = nm, term = names(v) %||% paste0("stat", seq_along(v)), value = as.numeric(v), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Missing-data sensitivity scenario registry
#' @param x An `nmaflow_data` object.
#' @param scenarios Named list of functions that each take and return an `nmaflow_data` object.
#' @return Named list of scenario data objects with an audit entry.
#' @export
nma_missing_sensitivity <- function(x, scenarios) {
  x <- .nma_unwrap(x)
  if (!is.list(scenarios) || !length(scenarios) || !all(vapply(scenarios, is.function, logical(1))))
    stop("`scenarios` must be a named list of functions.", call. = FALSE)
  if (is.null(names(scenarios))) names(scenarios) <- paste0("scenario", seq_along(scenarios))
  lapply(names(scenarios), function(nm) {
    z <- scenarios[[nm]](x)
    z$audit <- .nma_audit_bind(z$audit, .nma_audit_row("sensitivity", "missing_data_scenario",
                                                        reason = nm, automatic = FALSE))
    z
  }) |> stats::setNames(names(scenarios))
}

#' Zero-event continuity-correction sensitivity
#' @param x Arm-level data.
#' @param corrections Numeric corrections to compare.
#' @param event,total Event and denominator columns.
#' @param measure `RR` or `OR`.
#' @return Named list of contrast datasets.
#' @export
nma_zero_event_sensitivity <- function(x, corrections = c(0, 0.25, 0.5, 1), event, total,
                                       measure = c("RR", "OR")) {
  measure <- match.arg(measure)
  ans <- lapply(corrections, function(cc) nma_effects(x, measure = measure, event = event, total = total,
                                                       continuity = cc))
  stats::setNames(ans, paste0("cc_", corrections))
}
